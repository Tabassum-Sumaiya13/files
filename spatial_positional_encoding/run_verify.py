#!/usr/bin/env python3
"""
run_verify.py — Dataset-agnostic feature verification (survival NOT required).

Answers "do the pipeline's spatial features actually work on this cohort?" using
only the cohort itself — no survival, no outcome labels needed. Survival, when a
dataset has it, is a SEPARATE per-dataset downstream check; this tool is the
common one that runs for every dataset regardless.

THE BASELINE
    Composition proportions (fraction of each lineage / cell type per sample).
    This is the abundance-only null of the whole field: "spatial arrangement
    adds nothing over how much of each cell type is present." Every spatial
    feature must be shown to carry something this baseline does not.

THREE UNIVERSAL CHECKS (label-free)
  1. Null z          — feature vs its own within-sample spatial-shuffle null
                       (Delaunay graph fixed, cell labels permuted). |z| large
                       => the value reflects spatial organisation, not chance.
  2. Spatial-specific— 1 - R^2 of regressing the feature on the composition
                       baseline (5-fold). High => the feature is NOT just a
                       re-encoding of abundance; it adds beyond the baseline.
  3. Stability r     — split each sample's cells into two random halves, rebuild
                       the graph on each, recompute; Pearson r across samples.
                       High => reproducible, not a small-sample artifact.

OPTIONAL (only if a categorical label is available)
  4. Separability    — GroupKFold-by-patient classification AUC/accuracy of a
                       provided label, comparing baseline vs baseline+block vs a
                       width-matched noise block.

Usage
    python run_verify.py --dataset UPMC
    python run_verify.py --dataset UPMC --taxonomy native
    python run_verify.py --dataset CRC --taxonomy lineage --n-perm 20
    python run_verify.py --dataset UPMC --label-col <col in manifest>   # +separability
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
sys.path.insert(0, str(Path(__file__).parent / "src"))

from src import spatial_features as sf
from src.cohort import LINEAGES, load_cohort


def _banner(t: str):
    print(f"\n{'=' * 74}\n  {t}\n{'=' * 74}")


# ---------------------------------------------------------------------------
# Category vocabulary + fixed lineage grouping (shared across the whole cohort)
# ---------------------------------------------------------------------------
def build_vocab(cohort, label_col: str):
    """Return (categories, immune_idx, tumour_idx, stromal_idx).

    categories : fixed column order so every sample's vector aligns.
    *_idx      : category indices per lineage, decided once at cohort level by
                 the majority lineage of each category's cells (identity map when
                 label_col == 'lineage').
    """
    if label_col == "lineage":
        cats = list(LINEAGES)
        return cats, [0], [1], [2]

    lin_counter = {}
    cats_set = set()
    for _, df in cohort.iter_samples(columns=[label_col, "lineage"]):
        lab = df[label_col].astype(str).values
        lin = df["lineage"].astype(str).values
        for c, l in zip(lab, lin):
            cats_set.add(c)
            lin_counter.setdefault(c, {}).setdefault(l, 0)
            lin_counter[c][l] += 1
    cats = sorted(cats_set)
    buckets = {"immune": [], "tumour": [], "stromal": []}
    for j, c in enumerate(cats):
        counts = lin_counter.get(c, {})
        if counts:
            maj = max(counts, key=counts.get)
            if maj in buckets:
                buckets[maj].append(j)
    return cats, buckets["immune"], buckets["tumour"], buckets["stromal"]


def _idx_of(labels, cat_index):
    return np.fromiter((cat_index.get(l, -1) for l in labels), dtype=int, count=len(labels))


def _enrich(edges, idx, K, immune, tumour, stromal):
    counts = np.bincount(idx[idx >= 0], minlength=K).astype(float)
    M = sf._count_matrix(edges, idx, K)
    return sf.enrichment_scalars(M, counts, immune, tumour, stromal)


# ---------------------------------------------------------------------------
# Per-sample computation: real, permutation-null, split-half
# ---------------------------------------------------------------------------
def compute_cohort(cohort, label_col, cats, immune, tumour, stromal, n_perm, seed, limit):
    cat_index = {c: i for i, c in enumerate(cats)}
    K = len(cats)
    rng = np.random.RandomState(seed)

    sids, props, real, nullz, half1, half2 = [], [], [], [], [], []
    n_deg = 0
    read_cols = list(dict.fromkeys(["X", "Y", label_col, "lineage"]))  # dedupe (lineage==label_col)
    for sid, df in cohort.iter_samples(columns=read_cols, limit=limit):
        coords = df[["X", "Y"]].values.astype(float)
        labels = df[label_col].astype(str).values
        idx = _idx_of(labels, cat_index)
        edges = sf.delaunay_edges(coords)
        if edges.shape[0] and (idx < 0).any():
            keep = (idx[edges[:, 0]] >= 0) & (idx[edges[:, 1]] >= 0)
            edges = edges[keep]
        if edges.shape[0] == 0:
            n_deg += 1
            continue

        counts = np.bincount(idx[idx >= 0], minlength=K).astype(float)
        p = counts / counts.sum() if counts.sum() > 0 else counts
        r = _enrich(edges, idx, K, immune, tumour, stromal)

        # permutation null: fix graph, shuffle labels
        null = np.empty((n_perm, len(r)))
        for s in range(n_perm):
            null[s] = _enrich(edges, rng.permutation(idx), K, immune, tumour, stromal)
        nz = (r - null.mean(0)) / (null.std(0) + 1e-9)

        # split-half stability: two disjoint random halves, each its own graph
        order = rng.permutation(len(df))
        h = len(order) // 2
        hs = []
        for part in (order[:h], order[h:]):
            e2 = sf.delaunay_edges(coords[part])
            i2 = idx[part]
            if e2.shape[0] and (i2 < 0).any():
                k2 = (i2[e2[:, 0]] >= 0) & (i2[e2[:, 1]] >= 0)
                e2 = e2[k2]
            hs.append(_enrich(e2, i2, K, immune, tumour, stromal) if e2.shape[0] else np.full(len(r), np.nan))

        sids.append(sid)
        props.append(p)
        real.append(r)
        nullz.append(nz)
        half1.append(hs[0])
        half2.append(hs[1])

    return {
        "sids": sids,
        "props": np.array(props),
        "real": np.array(real),
        "nullz": np.array(nullz),
        "half1": np.array(half1),
        "half2": np.array(half2),
        "n_degenerate": n_deg,
    }


# ---------------------------------------------------------------------------
# Metric assembly
# ---------------------------------------------------------------------------
def spatial_specific(real: np.ndarray, props: np.ndarray) -> np.ndarray:
    """1 - CV R^2 of feature ~ composition. High = adds beyond abundance baseline."""
    from sklearn.linear_model import LinearRegression
    from sklearn.model_selection import KFold
    from sklearn.metrics import r2_score

    out = []
    kf = KFold(n_splits=5, shuffle=True, random_state=0)
    for f in range(real.shape[1]):
        y = real[:, f]
        if np.std(y) < 1e-9:
            out.append(0.0)
            continue
        preds = np.zeros_like(y)
        for tr, te in kf.split(props):
            preds[te] = LinearRegression().fit(props[tr], y[tr]).predict(props[te])
        # clamp composition-R^2 to [0,1]: R^2<0 means composition is anti-predictive,
        # i.e. explains none of the feature -> spatial_specific = 1.
        comp_r2 = float(np.clip(r2_score(y, preds), 0.0, 1.0))
        out.append(1.0 - comp_r2)
    return np.array(out)


def stability_r(half1: np.ndarray, half2: np.ndarray) -> np.ndarray:
    out = []
    for f in range(half1.shape[1]):
        a, b = half1[:, f], half2[:, f]
        ok = np.isfinite(a) & np.isfinite(b)
        if ok.sum() < 3 or np.std(a[ok]) < 1e-9 or np.std(b[ok]) < 1e-9:
            out.append(np.nan)
        else:
            out.append(float(np.corrcoef(a[ok], b[ok])[0, 1]))
    return np.array(out)


def build_matrix(res, feature_names) -> pd.DataFrame:
    real, nullz, props = res["real"], res["nullz"], res["props"]
    spec = spatial_specific(real, props)
    stab = stability_r(res["half1"], res["half2"])
    med_z = np.median(nullz, axis=0)
    frac_sig = np.mean(np.abs(nullz) > 2, axis=0)
    zero_var = np.std(real, axis=0) < 1e-9

    rows = []
    for f, name in enumerate(feature_names):
        rows.append({
            "feature": name,
            "real_mean": float(np.mean(real[:, f])),
            "null_z_median": float(med_z[f]),
            "frac_|z|>2": float(frac_sig[f]),
            "spatial_specific": float(spec[f]),
            "stability_r": float(stab[f]),
            "zero_variance": bool(zero_var[f]),
            "verdict": _verdict(med_z[f], spec[f], stab[f], zero_var[f]),
        })
    return pd.DataFrame(rows)


def _verdict(z, spec, stab, zero_var):
    if zero_var:
        return "DEGENERATE (constant)"
    real_signal = abs(z) >= 2
    adds = spec >= 0.5
    stable = (stab is not None) and np.isfinite(stab) and stab >= 0.5
    if real_signal and adds and stable:
        return "STRONG (real + adds + stable)"
    if real_signal and adds:
        return "REAL + adds beyond baseline"
    if real_signal:
        return "real spatial signal"
    return "weak / composition-like"


# ---------------------------------------------------------------------------
# Optional: separability against a provided categorical label
# ---------------------------------------------------------------------------
def separability(res, cohort, label_col, feature_names):
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.model_selection import StratifiedGroupKFold
    from sklearn.metrics import roc_auc_score, accuracy_score

    man = cohort.manifest.set_index("acquisition_id")
    if label_col not in man.columns:
        print(f"  [separability] label '{label_col}' not in manifest columns "
              f"{list(man.columns)} — skipped")
        return None
    sids = res["sids"]
    y = man.loc[sids, label_col]
    patients = man.loc[sids, "patient_id"].values
    ok = y.notna().values
    if ok.sum() < 20 or y[ok].nunique() < 2:
        print(f"  [separability] label '{label_col}' unusable (n={ok.sum()}, "
              f"classes={y[ok].nunique()}) — skipped")
        return None

    y = pd.Categorical(y[ok]).codes
    patients = patients[ok]
    props = res["props"][ok]
    block = res["real"][ok]
    rng = np.random.RandomState(42)
    noise = rng.randn(*block.shape)
    multiclass = len(np.unique(y)) > 2

    def score(X):
        n_splits = min(5, len(np.unique(patients)))
        skf = StratifiedGroupKFold(n_splits=n_splits, shuffle=True, random_state=0)
        s = []
        for tr, te in skf.split(X, y, patients):
            clf = RandomForestClassifier(n_estimators=300, random_state=0, n_jobs=-1)
            clf.fit(X[tr], y[tr])
            if multiclass:
                s.append(accuracy_score(y[te], clf.predict(X[te])))
            else:
                s.append(roc_auc_score(y[te], clf.predict_proba(X[te])[:, 1]))
        return float(np.mean(s))

    metric = "accuracy" if multiclass else "AUC"
    sets = {
        "Composition baseline": props,
        "Baseline + spatial block": np.hstack([props, block]),
        "Baseline + noise": np.hstack([props, noise]),
        "Spatial block alone": block,
    }
    rows = [{"feature_set": k, metric: score(v), "n_feats": v.shape[1]} for k, v in sets.items()]
    return pd.DataFrame(rows), metric


# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dataset", required=True, help="Ingested dataset folder name (UPMC, CRC, …)")
    ap.add_argument("--taxonomy", default="lineage", choices=["lineage", "native"],
                    help="lineage = portable 3-way (default, comparable across cohorts); "
                         "native = per-dataset cluster_label resolution")
    ap.add_argument("--n-perm", type=int, default=20, help="permutations for the null (default 20)")
    ap.add_argument("--label-col", default=None,
                    help="manifest column for the optional separability check")
    ap.add_argument("--limit", type=int, default=None, help="cap samples (fast iteration)")
    ap.add_argument("--seed", type=int, default=1029)
    args = ap.parse_args()

    _banner(f"FEATURE VERIFICATION — {args.dataset}  (taxonomy={args.taxonomy})")
    cohort = load_cohort(args.dataset)
    label_col = "lineage" if args.taxonomy == "lineage" else "cluster_label"
    print(f"  samples={len(cohort.sample_ids)}  markers={len(cohort.marker_columns)}  "
          f"survival_present={cohort.has_survival()}")

    cats, immune, tumour, stromal = build_vocab(cohort, label_col)
    print(f"  categories ({label_col}): {len(cats)}  "
          f"[immune={len(immune)}, tumour={len(tumour)}, stromal={len(stromal)}]")

    res = compute_cohort(cohort, label_col, cats, immune, tumour, stromal,
                         args.n_perm, args.seed, args.limit)
    print(f"  scored {len(res['sids'])} samples  ({res['n_degenerate']} skipped as degenerate)")

    matrix = build_matrix(res, sf.ENRICH_FEATURE_NAMES)

    _banner("VERIFICATION MATRIX (spatial enrichment block)")
    print(f"  baseline = composition proportions ({len(cats)} categories)\n")
    with pd.option_context("display.width", 200, "display.max_columns", None):
        print(matrix.to_string(index=False,
              formatters={c: (lambda v: f"{v:+.3f}") for c in
                          ["real_mean", "null_z_median", "spatial_specific", "stability_r"]}))
    print("\n  Read: null_z_median = spatial signal (|z|>=2 real); "
          "spatial_specific = 1-R^2 vs baseline (>=0.5 adds beyond abundance); "
          "stability_r = split-half reproducibility (>=0.5 good).")

    out_dir = cohort.processed_dir / "verification"
    out_dir.mkdir(parents=True, exist_ok=True)
    mpath = out_dir / f"verify_{args.taxonomy}.csv"
    matrix.to_csv(mpath, index=False)
    print(f"\n  Matrix saved: {mpath}")

    if args.label_col:
        _banner(f"OPTIONAL SEPARABILITY — predicting '{args.label_col}' (GroupKFold by patient)")
        sep = separability(res, cohort, args.label_col, sf.ENRICH_FEATURE_NAMES)
        if sep is not None:
            sep_df, metric = sep
            print()
            print(sep_df.to_string(index=False, formatters={metric: lambda v: f"{v:.3f}"}))
            sep_df.to_csv(out_dir / f"separability_{args.taxonomy}_{args.label_col}.csv", index=False)
            print(f"\n  A real block beats both 'Composition baseline' and 'Baseline + noise'.")


if __name__ == "__main__":
    main()
