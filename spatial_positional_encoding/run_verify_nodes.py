#!/usr/bin/env python3
"""
run_verify_nodes.py — Verification matrix for the NODE FEATURE block
(celltype-conditioned functional markers). Survival NOT required.

The enrichment block has had three label-free tests on two cohorts since
run_verify.py; the node-marker block has never had any. Its entire evidential
basis is one survival C-index on one cohort. This tool gives it the same
treatment, using only the cohort itself.

THE BASELINE
    Composition proportions — "you only need to know the cell-type mix, not what
    those cells are expressing." Every node feature must carry something this
    does not.

THREE LABEL-FREE CHECKS
  1. cond_z          — feature vs its own within-sample CELLTYPE-SHUFFLE null.
                       Permute the cell-type labels, keeping each cell's marker
                       values attached to it, and recompute. The conditioning set
                       becomes a random size-matched draw from the same tissue, so
                       the null value is ~the bulk mean of that marker.
                       => cond_z asks: IS THIS MARKER ACTUALLY ENRICHED IN THIS
                       CELL TYPE, or would reading it anywhere give the same
                       number? That is the premise celltype-conditioning rests on,
                       and it has never been tested. NOTE this is NOT a spatial
                       claim — node features are per-sample means and carry no
                       spatial information; the graph is not involved anywhere.
  2. composition_specific — 1 - R^2 of regressing the feature on the composition
                       baseline (5-fold). High => the marker state is not just a
                       re-encoding of which cell types are present.
  3. stability_r     — split each sample's cells into two random halves, recompute
                       the conditioned mean on each; Pearson r across samples.
                       High => the mean is estimated from enough cells to be
                       reproducible.

SUPPORT
    Reported, never imputed. A sample with no cells of the conditioning type is
    NaN and excluded, not zero-filled — "no such cells" is not "no signal".
    n_samples_supported / median_cells qualify every row.

Usage
    python run_verify_nodes.py --dataset UPMC
    python run_verify_nodes.py --dataset CRC --taxonomy native
    python run_verify_nodes.py --dataset CRC --n-perm 20
"""
from __future__ import annotations

import argparse
import importlib
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
sys.path.insert(0, str(Path(__file__).parent / "src"))

from src import node_features as nf
from src.cohort import LINEAGES, load_cohort
from run_verify import stability_r  # NaN-aware split-half correlation, shared


def _banner(t: str):
    print(f"\n{'=' * 74}\n  {t}\n{'=' * 74}")


def _load_dp_module(name: str):
    """data_preprocessing/ is a sibling package, not on the path by default."""
    dp = Path(__file__).resolve().parents[1] / "data_preprocessing"
    if str(dp) not in sys.path:
        sys.path.insert(0, str(dp))
    return importlib.import_module(name)


# ---------------------------------------------------------------------------
# Per-sample computation
# ---------------------------------------------------------------------------
def compute_cohort(cohort, conditioning, marker_map, taxonomy, n_perm, seed, limit):
    """real / null-z / split-half / composition for every sample."""
    rng = np.random.RandomState(seed)
    used_cols = sorted(set(marker_map.values()))
    read_cols = list(dict.fromkeys(["cluster_label", "lineage"] + used_cols))

    # Fixed baseline vocabulary so every sample's composition vector aligns.
    if taxonomy == "lineage":
        base_cats = list(LINEAGES)
    else:
        seen = set()
        for _, df in cohort.iter_samples(columns=["cluster_label"], limit=limit):
            seen.update(df["cluster_label"].astype(str).unique())
        base_cats = sorted(seen)
    base_col = "lineage" if taxonomy == "lineage" else "cluster_label"
    base_index = {c: i for i, c in enumerate(base_cats)}

    sids, props, real, sup, nullz, half1, half2 = [], [], [], [], [], [], []
    for sid, df in cohort.iter_samples(columns=read_cols, limit=limit):
        if len(df) < 2:
            continue
        labels = df["cluster_label"].astype(str).values
        values = {c: df[c].values.astype(float) for c in used_cols}

        # Integer-code the labels once; each feature becomes a boolean lookup over
        # codes, so a permutation is one O(n) gather instead of an isin per shuffle.
        cats, codes = np.unique(labels, return_inverse=True)
        cat_index = {c: i for i, c in enumerate(cats)}
        sel = np.zeros((len(nf.FEATURES), len(cats)), dtype=bool)
        for i, (name, _, _) in enumerate(nf.FEATURES):
            for lab in conditioning.get(name, []):
                j = cat_index.get(lab)
                if j is not None:
                    sel[i, j] = True

        def _feats(cd, pos=None):
            """Conditioned means over the cells addressed by `cd`.

            `pos` selects the same cells out of the marker arrays — required for
            the split-half, where `cd` covers only half the cells and the marker
            values must be subset to match.
            """
            vals = values if pos is None else {c: v[pos] for c, v in values.items()}
            out = np.full(len(nf.FEATURES), np.nan)
            n = np.zeros(len(nf.FEATURES), dtype=int)
            for i, (name, marker, _) in enumerate(nf.FEATURES):
                col = marker_map.get(marker)
                if col is None:
                    continue
                mask = sel[i][cd]
                k = int(mask.sum())
                n[i] = k
                if k:
                    out[i] = float(vals[col][mask].mean())
            return out, n

        r, n_sup = _feats(codes)

        # celltype-shuffle null: marker values stay on their cells, labels move
        null = np.full((n_perm, len(nf.FEATURES)), np.nan)
        for s in range(n_perm):
            null[s] = _feats(rng.permutation(codes))[0]
        with np.errstate(invalid="ignore"):
            nz = (r - np.nanmean(null, axis=0)) / (np.nanstd(null, axis=0) + 1e-9)

        # split-half: two disjoint random halves of the cells
        order = rng.permutation(len(df))
        h = len(order) // 2
        hs = [_feats(codes[part], part)[0] for part in (order[:h], order[h:])]

        bl = df[base_col].astype(str).values
        bidx = np.fromiter((base_index.get(v, -1) for v in bl), dtype=int, count=len(bl))
        cnt = np.bincount(bidx[bidx >= 0], minlength=len(base_cats)).astype(float)
        p = cnt / cnt.sum() if cnt.sum() > 0 else cnt

        sids.append(sid)
        props.append(p)
        real.append(r)
        sup.append(n_sup)
        nullz.append(nz)
        half1.append(hs[0])
        half2.append(hs[1])

    return {
        "sids": sids,
        "props": np.array(props),
        "real": np.array(real),
        "support": np.array(sup),
        "nullz": np.array(nullz),
        "half1": np.array(half1),
        "half2": np.array(half2),
        "base_cats": base_cats,
    }


# ---------------------------------------------------------------------------
# Metrics
# ---------------------------------------------------------------------------
def composition_specific(real: np.ndarray, props: np.ndarray) -> np.ndarray:
    """1 - CV R^2 of feature ~ composition, computed on SUPPORTED samples only.

    run_verify.spatial_specific assumes a fully populated column; node features
    legitimately have NaN where a sample lacks the conditioning cell type, so the
    mask is per feature rather than global.
    """
    from sklearn.linear_model import LinearRegression
    from sklearn.model_selection import KFold
    from sklearn.metrics import r2_score

    out = []
    for f in range(real.shape[1]):
        y = real[:, f]
        ok = np.isfinite(y)
        X = props[ok]
        y = y[ok]
        if len(y) < 10 or np.std(y) < 1e-9:
            out.append(np.nan)
            continue
        kf = KFold(n_splits=min(5, len(y)), shuffle=True, random_state=0)
        preds = np.zeros_like(y)
        for tr, te in kf.split(X):
            preds[te] = LinearRegression().fit(X[tr], y[tr]).predict(X[te])
        out.append(1.0 - float(np.clip(r2_score(y, preds), 0.0, 1.0)))
    return np.array(out)


def _verdict(z, spec, stab, frac_support, marker_present):
    if not marker_present:
        return "UNAVAILABLE (marker not in this panel)"
    if not np.isfinite(z):
        return "NO SUPPORT (no sample has this cell type)"
    if frac_support < 0.5:
        return f"LOW SUPPORT ({frac_support:.0%} of samples)"
    adds = np.isfinite(spec) and spec >= 0.5
    stable = np.isfinite(stab) and stab >= 0.5
    # A node feature asserts "marker M is characteristically expressed in celltype
    # T". z <= -2 says M is DEPLETED in T relative to a size-matched random draw
    # from the same tissue — the premise the feature is built on is false for this
    # cohort. That is a real, reproducible signal, but reporting it as STRONG
    # would invert its meaning, so it gets its own verdict.
    if z <= -2:
        return "CONTRADICTED (marker depleted in its own celltype)"
    if z >= 2 and adds and stable:
        return "STRONG (enriched + adds + stable)"
    if z >= 2 and adds:
        return "enriched + adds beyond baseline"
    if z >= 2:
        return "enriched in its celltype"
    return "weak / composition-like"


def build_matrix(res, marker_map) -> pd.DataFrame:
    real, nullz, sup = res["real"], res["nullz"], res["support"]
    spec = composition_specific(real, res["props"])
    stab = stability_r(res["half1"], res["half2"])
    n_samples = real.shape[0]

    rows = []
    for f, (name, marker, _) in enumerate(nf.FEATURES):
        supported = sup[:, f] > 0
        frac_sup = float(supported.mean()) if n_samples else 0.0
        z = nullz[:, f]
        z_ok = z[np.isfinite(z)]
        med_z = float(np.median(z_ok)) if z_ok.size else np.nan
        rows.append({
            "feature": name,
            "marker": marker,
            "n_samples_supported": int(supported.sum()),
            "frac_supported": frac_sup,
            "median_cells": float(np.median(sup[supported, f])) if supported.any() else 0.0,
            "real_mean": float(np.nanmean(real[:, f])) if np.isfinite(real[:, f]).any() else np.nan,
            "cond_z_median": med_z,
            "frac_|z|>2": float(np.mean(np.abs(z_ok) > 2)) if z_ok.size else np.nan,
            "composition_specific": float(spec[f]),
            "stability_r": float(stab[f]),
            "verdict": _verdict(med_z, spec[f], stab[f], frac_sup, marker in marker_map),
        })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dataset", required=True, help="Ingested dataset folder name (UPMC, CRC, …)")
    ap.add_argument("--taxonomy", default="lineage", choices=["lineage", "native"],
                    help="vocabulary for the COMPOSITION BASELINE only; conditioning is "
                         "always by Cell Ontology term (default lineage)")
    ap.add_argument("--n-perm", type=int, default=20, help="celltype-shuffle permutations (default 20)")
    ap.add_argument("--limit", type=int, default=None, help="cap samples (fast iteration)")
    ap.add_argument("--seed", type=int, default=1029)
    args = ap.parse_args()

    _banner(f"NODE FEATURE VERIFICATION — {args.dataset}  (baseline={args.taxonomy})")
    cohort = load_cohort(args.dataset)
    registry = _load_dp_module("registry")
    markers = _load_dp_module("markers")

    conditioning = nf.resolve_conditioning(registry.dataset_rows(args.dataset))
    marker_map, missing = nf.resolve_markers(cohort.marker_columns, markers)
    print(f"  samples={len(cohort.sample_ids)}  markers={len(cohort.marker_columns)}  "
          f"survival_present={cohort.has_survival()}")
    print(f"  markers resolved: {len(marker_map)}/{len(nf.CANONICAL_MARKERS)}"
          + (f"  MISSING: {missing}" if missing else ""))

    _banner("CONDITIONING — which native cell types each feature reads its marker in")
    summary = nf.conditioning_summary(conditioning)
    with pd.option_context("display.width", 200, "display.max_colwidth", 90):
        print(summary.to_string(index=False))
    print("\n  Resolved from celltype_registry.csv via Cell Ontology terms — not hand-mapped.")

    res = compute_cohort(cohort, conditioning, marker_map, args.taxonomy,
                         args.n_perm, args.seed, args.limit)
    if not res["sids"]:
        print("\n  [ERROR] no scoreable samples.")
        sys.exit(1)
    print(f"\n  scored {len(res['sids'])} samples")

    matrix = build_matrix(res, marker_map)

    _banner("NODE FEATURE VERIFICATION MATRIX")
    print(f"  baseline = composition proportions ({len(res['base_cats'])} categories)\n")
    with pd.option_context("display.width", 250, "display.max_columns", None):
        print(matrix.to_string(index=False, formatters={
            c: (lambda v: f"{v:+.3f}" if np.isfinite(v) else "     nan")
            for c in ["real_mean", "cond_z_median", "composition_specific", "stability_r"]}))
    print("\n  Read: cond_z_median = marker enrichment in its own celltype vs a size-matched "
          "random draw (|z|>=2 real);\n        composition_specific = 1-R^2 vs the celltype-mix "
          "baseline (>=0.5 adds); stability_r = split-half (>=0.5 good).")
    print("  cond_z is NOT a spatial claim — node features are per-sample means and use no graph.")

    out_dir = cohort.processed_dir / "verification"
    out_dir.mkdir(parents=True, exist_ok=True)
    mpath = out_dir / f"verify_nodes_{args.taxonomy}.csv"
    matrix.to_csv(mpath, index=False)
    summary.to_csv(out_dir / "node_conditioning.csv", index=False)
    print(f"\n  Matrix saved: {mpath}")
    print(f"  Conditioning saved: {out_dir / 'node_conditioning.csv'}")


if __name__ == "__main__":
    main()
