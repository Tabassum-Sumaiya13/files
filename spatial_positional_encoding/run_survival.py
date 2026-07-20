#!/usr/bin/env python3
"""
run_survival.py — Optional per-dataset survival downstream check.

This is NOT how the feature blocks are validated. Validation is survival-free and
lives in run_verify.py (enrichment) and run_verify_nodes.py (node markers), so a
cohort with no survival — CRC — is still fully verified. This tool answers a
separate, narrower question for the cohorts that do have outcome data:

    do the already-validated features also predict patient survival?

PROTOCOL
    RandomSurvivalForest, 100 trees, random_state=1029, GroupKFold by patient_id
    (identical hyperparameters to the retired legacy runner and to the spatsurv
    baseline it was matched against). Grouping by patient is what stops two
    regions from the same patient landing in train and test.

FEATURE SETS
    Composition baseline        cell-type proportions — the abundance-only null
    + Enrichment                the 5 verified spatial scalars
    + Nodes                     the 8 verified celltype-conditioned markers
    + Enrichment + Nodes        both blocks
    + Noise                     WIDTH-MATCHED random block — the control that
                                catches a block "winning" purely by widening the
                                feature matrix. A block that does not beat this
                                has not earned its columns.
    Enrichment alone / Nodes alone

Comparisons are PAIRED: every feature set is scored on the same folds, so the
per-fold delta against the baseline is a matched difference and its spread across
folds is reported rather than implied.

Usage
    python run_survival.py --dataset UPMC
    python run_survival.py --dataset UPMC --taxonomy native --n-splits 10
"""
from __future__ import annotations

# --- sklearn version pin (MUST run before any sklearn/sksurv import) ---------
# scikit-survival 0.27 (system-site) needs scikit-learn <= 1.8, but the user-site
# install carries 1.9.0 and shadows it:
#     ImportError: cannot import name 'DTYPE' from 'sklearn.tree._tree'
# Forcing system-site to the front resolves sklearn 1.8.0 and imports cleanly.
import sys
from pathlib import Path as _Path

_SYS_SITE = r"C:\Python312\Lib\site-packages"
if _Path(_SYS_SITE).exists():
    sys.path.insert(0, _SYS_SITE)

import argparse
import warnings

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

sys.path.insert(0, str(_Path(__file__).parent))
sys.path.insert(0, str(_Path(__file__).parent / "src"))

import importlib

from src import node_features as nf
from src import spatial_features as sf
from src.cohort import load_cohort
from run_verify import build_vocab

RSF_N_ESTIMATORS = 100
RSF_RANDOM_STATE = 1029
DEFAULT_N_SPLITS = 10


def _banner(t: str):
    print(f"\n{'=' * 74}\n  {t}\n{'=' * 74}")


def _load_dp_module(name: str):
    dp = _Path(__file__).resolve().parents[1] / "data_preprocessing"
    if str(dp) not in sys.path:
        sys.path.insert(0, str(dp))
    return importlib.import_module(name)


# ---------------------------------------------------------------------------
# Feature construction — one row per sample, both blocks at once
# ---------------------------------------------------------------------------
def build_features(cohort, dataset, taxonomy):
    """Return (sids, props, enrichment, nodes, node_support)."""
    registry = _load_dp_module("registry")
    markers = _load_dp_module("markers")

    label_col = "lineage" if taxonomy == "lineage" else "cluster_label"
    cats, immune, tumour, stromal = build_vocab(cohort, label_col)
    cat_index = {c: i for i, c in enumerate(cats)}
    K = len(cats)

    conditioning = nf.resolve_conditioning(registry.dataset_rows(dataset))
    marker_map, missing = nf.resolve_markers(cohort.marker_columns, markers)
    if missing:
        print(f"  [markers] not in this panel, their features will be NaN: {missing}")

    used_markers = sorted(set(marker_map.values()))
    read_cols = list(dict.fromkeys(["X", "Y", label_col, "lineage", "cluster_label"] + used_markers))

    sids, props, enr, nodes, sup = [], [], [], [], []
    for sid, df in cohort.iter_samples(columns=read_cols):
        if len(df) < 3:
            continue
        labels = df[label_col].astype(str).values
        idx = np.fromiter((cat_index.get(l, -1) for l in labels), dtype=int, count=len(labels))
        edges = sf.delaunay_edges(df[["X", "Y"]].values.astype(float))
        if edges.shape[0] and (idx < 0).any():
            keep = (idx[edges[:, 0]] >= 0) & (idx[edges[:, 1]] >= 0)
            edges = edges[keep]
        if edges.shape[0] == 0:
            continue
        counts = np.bincount(idx[idx >= 0], minlength=K).astype(float)
        M = sf._count_matrix(edges, idx, K)

        native = df["cluster_label"].astype(str).values
        values = {c: df[c].values.astype(float) for c in used_markers}
        nv, ns = nf.sample_node_features(native, values, conditioning, marker_map)

        sids.append(sid)
        props.append(counts / counts.sum() if counts.sum() > 0 else counts)
        enr.append(sf.enrichment_scalars(M, counts, immune, tumour, stromal))
        nodes.append(nv)
        sup.append(ns)

    return sids, np.array(props), np.array(enr), np.array(nodes), np.array(sup), cats


def _impute(block: np.ndarray, name: str) -> np.ndarray:
    """Column-median impute, reporting how much was imputed.

    Node features are legitimately NaN where a sample has no cells of the
    conditioning type. The forest cannot take NaN, so imputation is unavoidable
    here — but it is reported, not silent, because an imputed column is weaker
    evidence than a measured one.
    """
    out = block.copy()
    n_bad = int(np.isnan(out).sum())
    if n_bad:
        med = np.nanmedian(out, axis=0)
        med = np.where(np.isfinite(med), med, 0.0)
        inds = np.where(np.isnan(out))
        out[inds] = np.take(med, inds[1])
        print(f"  [impute] {name}: {n_bad} NaN cell(s) filled with the column median "
              f"({n_bad / out.size:.2%} of the block)")
    return out


# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dataset", required=True)
    ap.add_argument("--taxonomy", default="lineage", choices=["lineage", "native"])
    ap.add_argument("--n-splits", type=int, default=DEFAULT_N_SPLITS)
    ap.add_argument("--seed", type=int, default=RSF_RANDOM_STATE)
    args = ap.parse_args()

    from sklearn.model_selection import GroupKFold
    from sksurv.ensemble import RandomSurvivalForest
    from sksurv.metrics import concordance_index_censored

    _banner(f"SURVIVAL DOWNSTREAM — {args.dataset}  (taxonomy={args.taxonomy})")
    cohort = load_cohort(args.dataset)

    if not cohort.has_survival():
        print(f"  '{args.dataset}' has no usable survival in its manifest.\n"
              f"  This is NOT a failure: the feature blocks are validated without survival by\n"
              f"      python run_verify.py       --dataset {args.dataset}\n"
              f"      python run_verify_nodes.py --dataset {args.dataset}\n"
              f"  Survival is an optional extra for cohorts that have outcome data.")
        return

    man = cohort.manifest.set_index("acquisition_id")
    sids, props, enr, nodes, sup, cats = build_features(cohort, args.dataset, args.taxonomy)
    print(f"  scored {len(sids)} samples  |  composition baseline = {len(cats)} categories")

    nodes = _impute(nodes, "node block")
    t = pd.to_numeric(man.loc[sids, "survival_day"], errors="coerce").values
    e = pd.to_numeric(man.loc[sids, "survival_status"], errors="coerce").values
    patients = man.loc[sids, "patient_id"].astype(str).values

    ok = np.isfinite(t) & np.isfinite(e)
    if ok.sum() < 20 or e[ok].sum() < 5:
        print(f"  Too few usable samples/events ({ok.sum()} samples, {int(np.nansum(e[ok]))} events).")
        return
    t, e, patients = t[ok], e[ok].astype(bool), patients[ok]
    props, enr, nodes = props[ok], enr[ok], nodes[ok]
    y = np.array(list(zip(e, t)), dtype=[("event", bool), ("time", float)])
    print(f"  {len(t)} samples, {int(e.sum())} events ({e.mean():.1%}), "
          f"{len(np.unique(patients))} patients")

    rng = np.random.RandomState(args.seed)
    noise = rng.randn(len(t), enr.shape[1] + nodes.shape[1])

    sets = {
        "Composition baseline": props,
        "Baseline + Enrichment": np.hstack([props, enr]),
        "Baseline + Nodes": np.hstack([props, nodes]),
        "Baseline + Enrichment + Nodes": np.hstack([props, enr, nodes]),
        "Baseline + Noise (control)": np.hstack([props, noise]),
        "Enrichment alone": enr,
        "Nodes alone": nodes,
    }

    n_splits = min(args.n_splits, len(np.unique(patients)))
    folds = list(GroupKFold(n_splits=n_splits).split(props, y, patients))
    print(f"  GroupKFold by patient: {n_splits} folds\n")

    per_fold = {}
    for name, X in sets.items():
        cs = []
        for tr, te in folds:
            if y["event"][tr].sum() == 0 or y["event"][te].sum() == 0:
                cs.append(np.nan)
                continue
            m = RandomSurvivalForest(n_estimators=RSF_N_ESTIMATORS,
                                     random_state=RSF_RANDOM_STATE, n_jobs=-1)
            m.fit(X[tr], y[tr])
            cs.append(concordance_index_censored(y["event"][te], y["time"][te],
                                                 m.predict(X[te]))[0])
        per_fold[name] = np.array(cs, dtype=float)
        print(f"  {name:32s} C = {np.nanmean(per_fold[name]):.3f}")

    base = per_fold["Composition baseline"]
    ctrl = per_fold["Baseline + Noise (control)"]
    rows = []
    for name, cs in per_fold.items():
        d_base = cs - base
        d_ctrl = cs - ctrl
        rows.append({
            "feature_set": name,
            "n_feats": sets[name].shape[1],
            "C_index": float(np.nanmean(cs)),
            "C_sd_across_folds": float(np.nanstd(cs)),
            "delta_vs_baseline": float(np.nanmean(d_base)),
            "delta_vs_baseline_sd": float(np.nanstd(d_base)),
            "delta_vs_noise": float(np.nanmean(d_ctrl)),
            "beats_baseline_and_noise": bool(np.nanmean(d_base) > 0 and np.nanmean(d_ctrl) > 0),
        })
    out = pd.DataFrame(rows)

    _banner("SURVIVAL RESULT (paired over the same folds)")
    with pd.option_context("display.width", 220, "display.max_columns", None):
        print(out.to_string(index=False, formatters={
            c: (lambda v: f"{v:+.3f}") for c in
            ["C_index", "C_sd_across_folds", "delta_vs_baseline",
             "delta_vs_baseline_sd", "delta_vs_noise"]}))
    print("\n  A block has earned its columns only if it beats BOTH the composition\n"
          "  baseline AND the width-matched noise control.")

    out_dir = cohort.processed_dir / "verification"
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / f"survival_{args.taxonomy}.csv"
    out.to_csv(path, index=False)
    print(f"\n  Saved: {path}")


if __name__ == "__main__":
    main()
