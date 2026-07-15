#!/usr/bin/env python3
"""
run_survival.py — Survival prediction comparison: PE features vs baselines.

Matches the spatsurv-main baseline approach:
  - RandomSurvivalForest (scikit-survival)
  - GroupKFold by patient_id (no data leakage)
  - concordance_index_censored for evaluation
  - Same RSF hyperparameters (n_estimators=100, random_state=1029)

Feature sets compared:
  1. Celltype proportions (baseline from spatsurv)
  2. PE features (our contribution)
  3. PE + Celltype proportions (combined)
  4. Raw X,Y pooled stats (sanity check)

Usage:
    python run_survival.py
    python run_survival.py --pooling mean_std
"""

import sys
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from typing import Dict, List, Tuple
from sklearn.model_selection import GroupKFold
from sklearn.preprocessing import StandardScaler
import warnings
warnings.filterwarnings('ignore')

sys.path.insert(0, str(Path(__file__).parent))

from config import get_config


def _banner(title: str):
    print(f"\n{'=' * 70}")
    print(f"  {title}")
    print(f"{'=' * 70}")


# ======================================================================
# Data Loading
# ======================================================================

def load_survival_labels(cfg) -> pd.DataFrame:
    """Load survival data from sample_metadata.csv."""
    _banner("LOADING SURVIVAL LABELS")

    meta = pd.read_csv(cfg.metadata_path)
    print(f"  Raw metadata: {meta.shape}")
    print(f"  Columns: {list(meta.columns)}")

    # Keep needed columns
    surv = meta[[
        cfg.META_ACQ_COL, cfg.META_PATIENT_COL, cfg.META_COVERSLIP_COL,
        cfg.META_SURVIVAL_TIME_COL, cfg.META_SURVIVAL_STATUS_COL,
    ]].copy()

    surv.columns = ['acquisition_id', 'patient_id', 'coverslip_label',
                     'survival_time', 'event']

    # Drop missing
    before = len(surv)
    surv = surv.dropna(subset=['survival_time', 'event'])
    surv = surv[surv['survival_time'] > 0]
    surv['event'] = surv['event'].astype(int)
    surv = surv[surv['event'].isin([0, 1])]
    after = len(surv)

    print(f"\n  After cleaning: {after} samples (dropped {before - after})")
    print(f"  Events (deaths): {surv['event'].sum()} / {after} "
          f"({surv['event'].mean() * 100:.1f}%)")
    print(f"  Survival time: [{surv['survival_time'].min():.0f}, "
          f"{surv['survival_time'].max():.0f}] days")
    print(f"  Unique patients: {surv['patient_id'].nunique()}")

    return surv


def load_celltype_proportions(cfg) -> pd.DataFrame:
    """
    Compute celltype proportions per sample from expression data —
    matching spatsurv's baseline exactly.
    """
    _banner("COMPUTING CELLTYPE PROPORTIONS (spatsurv baseline)")

    expr_path = cfg.expression_path
    print(f"  Reading expression data: {expr_path}")
    print(f"  (This may take a moment for large files...)")

    # Only need sample_id, cluster columns
    expr = pd.read_csv(expr_path, usecols=[cfg.EXPR_SAMPLE_COL, cfg.EXPR_CLUSTER_COL])
    expr = expr.rename(columns={
        cfg.EXPR_SAMPLE_COL: 'acquisition_id',
        cfg.EXPR_CLUSTER_COL: 'cluster',
    })
    expr['acquisition_id'] = expr['acquisition_id'].astype(str)

    print(f"  Expression rows: {len(expr):,}")
    print(f"  Unique samples: {expr['acquisition_id'].nunique()}")
    print(f"  Unique clusters: {expr['cluster'].nunique()}")

    # Compute proportions: (samples × clusters) normalised by row sum
    prop = (
        expr.groupby(['acquisition_id', 'cluster']).size()
        .unstack(fill_value=0)
    )
    prop = prop.div(prop.sum(axis=1), axis=0)

    print(f"  Proportion matrix: {prop.shape}")
    print(f"  Columns (cluster IDs): {list(prop.columns)}")

    # Rename columns
    prop.columns = [f'celltype_{c}' for c in prop.columns]

    print(f"  Head:\n{prop.head(3).to_string()}")

    return prop


def pool_pe_features(
    cfg,
    sample_ids: List[str],
    pooling: str = "mean_std",
) -> pd.DataFrame:
    """
    Pool per-cell PE features into per-sample features.

    pooling methods:
      - "mean": mean of each PE dim
      - "mean_std": concatenate mean + std of each PE dim
    """
    _banner("POOLING PE FEATURES -> PER-SAMPLE")

    encoding_dir = cfg.encoding_dir
    features = {}

    for sid in sample_ids:
        enc_file = encoding_dir / f"encoding_{sid}.parquet"
        if not enc_file.exists():
            continue

        df = pd.read_parquet(enc_file)
        pe_cols = [c for c in df.columns if c.startswith('pe_')]

        # Use only valid PE cells
        if 'is_valid_pe' in df.columns:
            df = df[df['is_valid_pe']]

        if len(df) == 0 or len(pe_cols) == 0:
            continue

        X = df[pe_cols].values
        if pooling == 'mean':
            pooled = X.mean(axis=0)
        elif pooling == 'mean_std':
            pooled = np.concatenate([X.mean(axis=0), X.std(axis=0)])
        else:
            raise ValueError(f"Unknown pooling: {pooling}")

        features[sid] = pooled

    if not features:
        print("  [WARNING] No PE features loaded!")
        return pd.DataFrame()

    result = pd.DataFrame.from_dict(features, orient='index')
    result.index.name = 'acquisition_id'
    n_dims = result.shape[1]

    if pooling == 'mean_std':
        k_pe = n_dims // 2
        col_names = [f'pe_mean_{i}' for i in range(k_pe)] + [f'pe_std_{i}' for i in range(k_pe)]
    else:
        col_names = [f'pe_mean_{i}' for i in range(n_dims)]
    result.columns = col_names

    print(f"  Pooled PE: {result.shape[0]} samples × {result.shape[1]} features")
    print(f"  Pooling method: {pooling}")

    return result


def pool_raw_coord_features(
    cfg,
    sample_ids: List[str],
) -> pd.DataFrame:
    """Pool raw X, Y coordinates as a sanity-check baseline."""
    _banner("POOLING RAW COORDINATE FEATURES")

    features = {}
    for sid in sample_ids:
        sample_file = cfg.samples_dir / f"sample_{sid}.parquet"
        if not sample_file.exists():
            continue
        df = pd.read_parquet(sample_file)
        X = df[['X', 'Y']].values
        pooled = np.concatenate([X.mean(axis=0), X.std(axis=0)])
        features[sid] = pooled

    if not features:
        return pd.DataFrame()

    result = pd.DataFrame.from_dict(features, orient='index')
    result.index.name = 'acquisition_id'
    result.columns = ['raw_X_mean', 'raw_Y_mean', 'raw_X_std', 'raw_Y_std']

    print(f"  Raw coord features: {result.shape}")
    return result


# ======================================================================
# RSF Training & Evaluation
# ======================================================================

def run_rsf_comparison(
    survival: pd.DataFrame,
    feature_sets: Dict[str, pd.DataFrame],
    cfg,
    output_dir: Path,
) -> Dict:
    """
    Train RandomSurvivalForest for each feature set, evaluate with C-index.

    Matches spatsurv: GroupKFold by patient_id, RSF(n_estimators=100),
    concordance_index_censored.
    """
    try:
        from sksurv.ensemble import RandomSurvivalForest
        from sksurv.metrics import concordance_index_censored
        print("\n  [OK] scikit-survival loaded")
    except ImportError:
        print("\n  [ERROR] scikit-survival not installed!")
        print("  Install it: pip install scikit-survival")
        return {}

    _banner("RANDOM SURVIVAL FOREST COMPARISON")

    results = {}

    for name, feat_df in feature_sets.items():
        if feat_df.empty:
            print(f"\n  [SKIP] {name}: no features")
            continue

        # Align: samples present in both survival + features
        common = sorted(set(feat_df.index) & set(survival['acquisition_id']))
        if len(common) < 20:
            print(f"\n  [SKIP] {name}: only {len(common)} common samples (need >=20)")
            continue

        # Build arrays
        feat_aligned = feat_df.loc[common]
        surv_aligned = survival.set_index('acquisition_id').loc[common]

        X = feat_aligned.values
        y_time = surv_aligned['survival_time'].values
        y_event = surv_aligned['event'].values.astype(bool)
        patients = surv_aligned['patient_id'].values

        # Structured survival array
        y_surv = np.array(
            list(zip(y_event, y_time)),
            dtype=[('event', bool), ('time', '<f8')]
        )

        # Scale features
        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(X)

        # GroupKFold by patient
        n_unique_patients = len(np.unique(patients))
        n_splits = min(cfg.rsf_n_splits, n_unique_patients)
        gkf = GroupKFold(n_splits=n_splits)

        print(f"\n  {name}: {len(common)} samples, {X.shape[1]} features, "
              f"{n_splits} folds (by patient)")

        cv_scores = []
        for fold_i, (train_idx, test_idx) in enumerate(gkf.split(X_scaled, y_surv, patients)):
            try:
                rsf = RandomSurvivalForest(
                    n_estimators=cfg.rsf_n_estimators,
                    n_jobs=-1,
                    random_state=cfg.rsf_random_state,
                )
                rsf.fit(X_scaled[train_idx], y_surv[train_idx])

                # Predict risk scores
                risk = rsf.predict(X_scaled[test_idx])

                # Concordance index
                c_idx = concordance_index_censored(
                    y_event[test_idx], y_time[test_idx], risk
                )[0]

                cv_scores.append(c_idx)
                print(f"    Fold {fold_i}: C-index = {c_idx:.4f} "
                      f"(train={len(train_idx)}, test={len(test_idx)})")

            except Exception as e:
                print(f"    Fold {fold_i}: FAILED — {e}")
                cv_scores.append(np.nan)

        mean_c = np.nanmean(cv_scores)
        std_c = np.nanstd(cv_scores)

        results[name] = {
            'c_index_mean': mean_c,
            'c_index_std': std_c,
            'c_index_per_fold': cv_scores,
            'n_samples': len(common),
            'n_features': X.shape[1],
            'n_folds': n_splits,
        }

        print(f"    >> {name}: C-index = {mean_c:.4f} +/- {std_c:.4f}")

    return results


def print_final_summary(results: Dict, output_dir: Path):
    """Print and save the final comparison table."""
    _banner("FINAL RESULTS")

    if not results:
        print("  No results to report!")
        return

    # Build summary table
    rows = []
    for name, r in sorted(results.items(), key=lambda x: x[1]['c_index_mean'], reverse=True):
        rows.append({
            'Feature Set': name,
            'C-index': f"{r['c_index_mean']:.4f} ± {r['c_index_std']:.4f}",
            'C-index (mean)': r['c_index_mean'],
            'C-index (std)': r['c_index_std'],
            'Samples': r['n_samples'],
            'Features': r['n_features'],
        })

    summary_df = pd.DataFrame(rows)

    print(f"\n  {'Feature Set':<35s} {'C-index':<20s} {'Samples':<10s} {'Features':<10s}")
    print(f"  {'-' * 75}")
    for _, row in summary_df.iterrows():
        print(f"  {row['Feature Set']:<35s} {row['C-index']:<20s} "
              f"{row['Samples']:<10} {row['Features']:<10}")

    # Compare to spatsurv baseline
    baseline_c = 0.704
    best = summary_df.iloc[0]
    print(f"\n  Spatsurv baseline C-index: {baseline_c}")
    print(f"  Best model: {best['Feature Set']} = {best['C-index']}")
    if best['C-index (mean)'] > baseline_c:
        print(f"  [OK] BEATS baseline by {best['C-index (mean)'] - baseline_c:.4f}")
    else:
        print(f"  [!!] Below baseline by {baseline_c - best['C-index (mean)']:.4f}")

    # Save
    save_path = output_dir / "survival_results.csv"
    summary_df.to_csv(save_path, index=False)
    print(f"\n  Results saved: {save_path}")


# ======================================================================
# Main
# ======================================================================

def main():
    parser = argparse.ArgumentParser(description="Survival RSF Evaluation")
    parser.add_argument('--pooling', default='mean_std',
                        choices=['mean', 'mean_std'],
                        help="How to pool cell->sample features")
    args = parser.parse_args()

    cfg = get_config(pooling_method=args.pooling)
    cfg.print_summary()

    output_dir = cfg.processed_dir / "outputs" / "survival"
    output_dir.mkdir(parents=True, exist_ok=True)

    # --- Load survival labels ---
    survival = load_survival_labels(cfg)

    # --- Build feature sets ---
    sample_ids = survival['acquisition_id'].tolist()

    feature_sets = {}

    # 1. Celltype proportions (spatsurv baseline)
    try:
        celltype_prop = load_celltype_proportions(cfg)
        feature_sets['Celltype Proportions'] = celltype_prop
    except Exception as e:
        print(f"  [WARNING] Could not load celltype props: {e}")

    # 2. PE features
    pe_features = pool_pe_features(cfg, sample_ids, pooling=args.pooling)
    if not pe_features.empty:
        feature_sets['PE (Laplacian)'] = pe_features

    # 3. Raw coordinate stats (sanity check)
    raw_features = pool_raw_coord_features(cfg, sample_ids)
    if not raw_features.empty:
        feature_sets['Raw Coordinates'] = raw_features

    # 4. Combined: PE + Celltype
    if 'Celltype Proportions' in feature_sets and 'PE (Laplacian)' in feature_sets:
        combined = pd.concat(
            [feature_sets['Celltype Proportions'], feature_sets['PE (Laplacian)']],
            axis=1, join='inner',
        )
        if not combined.empty:
            feature_sets['PE + Celltype'] = combined

    # 5. Combined: PE + Raw
    if 'PE (Laplacian)' in feature_sets and 'Raw Coordinates' in feature_sets:
        combined = pd.concat(
            [feature_sets['PE (Laplacian)'], feature_sets['Raw Coordinates']],
            axis=1, join='inner',
        )
        if not combined.empty:
            feature_sets['PE + Raw Coords'] = combined

    print(f"\n  Feature sets to evaluate: {list(feature_sets.keys())}")

    if len(feature_sets) < 1:
        print("  [ERROR] No feature sets available — run encoding pipeline first!")
        return

    # --- Run RSF ---
    results = run_rsf_comparison(survival, feature_sets, cfg, output_dir)

    # --- Summary ---
    print_final_summary(results, output_dir)


if __name__ == "__main__":
    main()
