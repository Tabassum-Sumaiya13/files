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
from pathlib import Path as _Path

# --- sklearn version fix (must run before ANY sklearn/sksurv import) ---
# The user-site sklearn (1.9) is missing symbols scikit-survival needs; the
# system-site install (1.8) is compatible. Force it to the front of the path.
_SYS_SITE = r'C:\Python312\Lib\site-packages'
if _Path(_SYS_SITE).exists():
    sys.path.insert(0, _SYS_SITE)

import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from typing import Dict, List, Tuple
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

def load_survival_labels(cfg, exclude_normal: bool = False,
                         censor_doc: bool = False) -> pd.DataFrame:
    """
    Load survival data from sample_metadata.csv.

    exclude_normal : drop Normal mucosa acquisitions. They contain no tumour, so
        their spatial architecture is unrelated to the outcome, yet they inherit
        the patient's survival label -> pure label noise. The base paper filters
        them (rsf_prediction.py --exclude-normal); this pipeline never did.

    censor_doc : treat "dead of other causes" as CENSORED rather than as an event
        (cause-specific hazard). `survival_status == 1` lumps DOD (dead of
        disease) together with DOC (heart attack, etc). No tumour feature can
        predict a patient dying of something unrelated, so counting DOC as an
        event puts a hard ceiling on any achievable C-index.
    """
    _banner("LOADING SURVIVAL LABELS")

    meta = pd.read_csv(cfg.metadata_path)
    print(f"  Raw metadata: {meta.shape}")

    keep = [cfg.META_ACQ_COL, cfg.META_PATIENT_COL, cfg.META_COVERSLIP_COL,
            cfg.META_SURVIVAL_TIME_COL, cfg.META_SURVIVAL_STATUS_COL]
    extra = [c for c in ['tissuetype', 'status'] if c in meta.columns]
    surv = meta[keep + extra].copy()
    surv.columns = ['acquisition_id', 'patient_id', 'coverslip_label',
                    'survival_time', 'event'] + extra

    # --- cohort definition ---
    if exclude_normal and 'tissuetype' in surv.columns:
        before = len(surv)
        surv = surv[surv['tissuetype'] != 'Normal mucosa']
        print(f"  [exclude-normal] dropped {before - len(surv)} Normal mucosa "
              f"-> {len(surv)} rows")

    # --- label definition ---
    if censor_doc:
        if 'status' not in surv.columns:
            print("  [WARNING] --censor-doc requested but no 'status' column; skipping")
        else:
            n_doc = (surv['status'] == 'DOC').sum()
            # cause-specific: ONLY dead-of-disease counts as an event
            surv['event'] = (surv['status'] == 'DOD').astype(float)
            surv.loc[surv['status'].isna(), 'event'] = np.nan
            print(f"  [censor-doc] {n_doc} 'dead of other causes' re-coded as censored")

    # --- drop missing ---
    before = len(surv)
    surv = surv.dropna(subset=['survival_time', 'event'])
    surv = surv[surv['survival_time'] > 0]
    surv['event'] = surv['event'].astype(int)
    surv = surv[surv['event'].isin([0, 1])]

    print(f"\n  After cleaning: {len(surv)} samples (dropped {before - len(surv)})")
    print(f"  Survival time: [{surv['survival_time'].min():.0f}, "
          f"{surv['survival_time'].max():.0f}] days")

    # Report at BOTH levels. Acquisition-level counts are pseudo-replication:
    # 3 acquisitions from one patient are not 3 independent events. The
    # patient-level event count is what actually limits the model.
    pat = surv.groupby('patient_id')['event'].max()
    print(f"\n  Acquisition level: {len(surv)} samples, {surv['event'].sum()} events")
    print(f"  PATIENT level    : {len(pat)} patients, {int(pat.sum())} events "
          f"<- this is the real sample size")
    print(f"  Events per variable at 272 feats: {pat.sum() / 272:.3f} (rule of thumb: >=10)")

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


def load_cell_count(
    cfg,
    sample_ids: List[str],
) -> pd.DataFrame:
    """
    Trivial 1-feature control: number of cells per sample.

    This is the number every degenerate feature block secretly collapses into
    (see the critical review). Any real feature MUST beat this to mean anything.
    """
    _banner("LOADING CELL COUNT (trivial control)")

    features = {}
    for sid in sample_ids:
        sample_file = cfg.samples_dir / f"sample_{sid}.parquet"
        if not sample_file.exists():
            continue
        df = pd.read_parquet(sample_file, columns=['X'])
        features[sid] = [len(df)]

    if not features:
        print("  [WARNING] No sample files found for cell count!")
        return pd.DataFrame()

    result = pd.DataFrame.from_dict(features, orient='index')
    result.index.name = 'acquisition_id'
    result.columns = ['cell_count']
    print(f"  Cell count: {result.shape[0]} samples, "
          f"range [{int(result['cell_count'].min())}, {int(result['cell_count'].max())}]")
    return result


def make_noise_block(index, n_features: int, seed: int = 0) -> pd.DataFrame:
    """
    A block of pure random numbers, same width as a block under test.

    This is THE control for "does my block add?". On a wide baseline the forest
    samples ~sqrt(p) features per split, so ANY added block is mostly ignored and
    scores slightly worse than the baseline through dilution alone. Comparing a
    real block against dC > 0 therefore condemns good blocks. The right reference
    is a noise block of IDENTICAL width: real must beat noise, not beat zero.
    """
    rng = np.random.RandomState(seed)
    return pd.DataFrame(
        rng.randn(len(index), n_features),
        index=index,
        columns=[f'noise_{i}' for i in range(n_features)],
    )


def _is_noise(name: str) -> bool:
    n = name.lower()
    return n.startswith('noise') or '+ noise' in n


def _add_width_matched_noise(feature_sets: Dict[str, pd.DataFrame],
                             baseline_name: str,
                             common: List[str],
                             seed: int = 42) -> Dict[str, pd.DataFrame]:
    """Ensure every `baseline + block` set has a noise control of IDENTICAL width.

    A 272-wide block judged against a 28-wide noise control measures the width
    difference, not the signal. Each real block gets a control built the same way
    it was: baseline + (W - baseline_width) random columns.
    """
    base_df = feature_sets.get(baseline_name)
    if base_df is None:
        return feature_sets

    base_w = base_df.loc[common].shape[1]
    have = {feature_sets[nm].loc[common].shape[1]
            for nm in feature_sets if _is_noise(nm)}
    need = {feature_sets[nm].loc[common].shape[1]
            for nm in feature_sets
            if nm != baseline_name and not _is_noise(nm)}

    out = dict(feature_sets)
    for w in sorted(need - have):
        k = w - base_w
        if k <= 0:
            continue  # standalone block narrower than the baseline; no control
        noise = make_noise_block(base_df.index, k, seed=seed)
        out[f'Celltype + Noise({k})'] = pd.concat([base_df, noise], axis=1, join='inner')
        print(f"  [auto-noise] width-matched control added: Celltype + Noise({k}) "
              f"-> {w} feats")
    return out


def load_neighbor_features(cfg, graph_name: str) -> pd.DataFrame:
    """Load a precomputed row-normalised neighbour matrix (see src/neighbor_features.py)."""
    path = cfg.processed_dir / "neighbor_features" / f"neighbor_{graph_name}.parquet"
    if not path.exists():
        print(f"  [WARNING] {path.name} not found — run: python src/neighbor_features.py")
        return pd.DataFrame()
    df = pd.read_parquet(path)
    print(f"  [{graph_name:9s}] {df.shape[0]} samples x {df.shape[1]} features")
    return df


def load_marker_states(cfg) -> pd.DataFrame:
    """Load the celltype-conditioned functional-marker node features.

    See src/marker_states.py — 8 validated functional markers (Ki67, PDL1,
    GranzymeB, PD1, ICOS, FoxP3, CD45RO, HLA-DR), each read in the celltype it
    matters in. This is the only place per-cell marker expression (the node
    feature) enters the survival model.
    """
    path = cfg.processed_dir / "neighbor_features" / "marker_states.parquet"
    if not path.exists():
        print(f"  [WARNING] {path.name} not found — run: python src/marker_states.py")
        return pd.DataFrame()
    df = pd.read_parquet(path)
    print(f"  [markers  ] {df.shape[0]} samples x {df.shape[1]} features: {list(df.columns)}")
    return df


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

def _load_sksurv():
    """Import scikit-survival, working around the known sklearn version clash."""
    try:
        from sksurv.ensemble import RandomSurvivalForest
        from sksurv.metrics import concordance_index_censored
    except ImportError:
        # Known env issue (see critical review): user-site sklearn conflicts with
        # sksurv; the system-site install works. Retry after putting it first.
        sys.path.insert(0, r'C:\Python312\Lib\site-packages')
        from sksurv.ensemble import RandomSurvivalForest
        from sksurv.metrics import concordance_index_censored
    return RandomSurvivalForest, concordance_index_censored


def _surv_array(event: np.ndarray, time: np.ndarray) -> np.ndarray:
    """Build the structured (event, time) array scikit-survival expects."""
    return np.array(list(zip(event, time)), dtype=[('event', bool), ('time', '<f8')])


def evaluate_seeded(
    X: np.ndarray,
    y_event: np.ndarray,
    y_time: np.ndarray,
    patients: np.ndarray,
    cfg,
    n_seeds: int,
    permute: bool = False,
) -> np.ndarray:
    """
    Repeat StratifiedGroupKFold CV over `n_seeds` different shuffles and return
    the per-seed mean C-index (one number per seed).

    - Grouped by patient (no leakage), stratified on event (kills empty folds).
    - permute=True shuffles the labels first -> gives the NOISE FLOOR (null).
    - No feature scaling: RSF splits on ranks, so it is scale-invariant, and
      fitting a scaler on the full set would be the leakage pattern anyway.
    """
    from sklearn.model_selection import StratifiedGroupKFold

    RSF, cidx = _load_sksurv()
    n_splits = min(cfg.rsf_n_splits, len(np.unique(patients)))

    seed_means = []
    for seed in range(n_seeds):
        ev, tm = y_event, y_time
        if permute:
            rng = np.random.RandomState(9000 + seed)
            perm = rng.permutation(len(y_event))
            ev, tm = y_event[perm], y_time[perm]

        y_surv = _surv_array(ev, tm)
        skf = StratifiedGroupKFold(n_splits=n_splits, shuffle=True, random_state=seed)

        fold_scores = []
        for train_idx, test_idx in skf.split(X, ev.astype(int), patients):
            if ev[test_idx].sum() == 0:
                continue  # C-index undefined with no events in the test fold
            rsf = RSF(
                n_estimators=cfg.rsf_n_estimators,
                n_jobs=-1,
                random_state=cfg.rsf_random_state,
            )
            rsf.fit(X[train_idx], y_surv[train_idx])
            risk = rsf.predict(X[test_idx])
            fold_scores.append(cidx(ev[test_idx], tm[test_idx], risk)[0])

        if fold_scores:
            seed_means.append(float(np.mean(fold_scores)))

    return np.array(seed_means)


def _ci(scores: np.ndarray) -> Tuple[float, float, float]:
    """Return (mean, 2.5th pct, 97.5th pct) of a seed-score array."""
    if len(scores) == 0:
        return (np.nan, np.nan, np.nan)
    return (float(scores.mean()),
            float(np.percentile(scores, 2.5)),
            float(np.percentile(scores, 97.5)))


def run_validation(
    survival: pd.DataFrame,
    feature_sets: Dict[str, pd.DataFrame],
    cfg,
    output_dir: Path,
    n_seeds: int,
    n_null_seeds: int,
    baseline_name: str = 'Celltype Proportions',
    tag: str = '',
) -> pd.DataFrame:
    """
    Honest validation: for every feature set report
      - real C-index with a 95% CI over `n_seeds` seeds
      - its own permutation null (noise floor)
      - paired ΔC vs the baseline, over the SAME seeds/folds

    All feature sets are aligned to one shared cohort so the CV folds are
    identical across sets -> ΔC per seed is a true paired difference.
    """
    _banner("VALIDATION: real vs null vs baseline (seeded, with CIs)")

    # --- shared cohort so folds match across every feature set ---
    surv_idx = set(survival['acquisition_id'])
    common = surv_idx.copy()
    for feat_df in feature_sets.values():
        if not feat_df.empty:
            common &= set(feat_df.index)
    common = sorted(common)

    if len(common) < 20:
        print(f"  [ERROR] Only {len(common)} samples shared across all feature sets.")
        return pd.DataFrame()

    surv_aligned = survival.set_index('acquisition_id').loc[common]
    y_time = surv_aligned['survival_time'].values.astype(float)
    y_event = surv_aligned['event'].values.astype(bool)
    patients = surv_aligned['patient_id'].values

    print(f"  Shared cohort: {len(common)} samples, "
          f"{len(np.unique(patients))} patients, {int(y_event.sum())} events")
    print(f"  Seeds: {n_seeds} (real), {n_null_seeds} (null)")

    # --- width-matched dilution controls ---
    # Every real block needs a noise block of its OWN width. Without this the
    # verdict depended on which noise row happened to exist: `Celltype + delaunay`
    # (272 feats, C=0.6899) read 'REAL SIGNAL' against Noise(256) in Experiment A
    # and 'indistinguishable from noise' against Noise(12) in Experiment B —
    # identical score, opposite verdict.
    feature_sets = _add_width_matched_noise(feature_sets, baseline_name, common)

    # --- real scores first (needed as baseline for ΔC) ---
    real_scores: Dict[str, np.ndarray] = {}
    null_scores: Dict[str, np.ndarray] = {}
    for name, feat_df in feature_sets.items():
        if feat_df.empty:
            continue
        X = feat_df.loc[common].values
        print(f"\n  [{name}] {X.shape[1]} features -> running {n_seeds}+{n_null_seeds} fits/seed...")
        real_scores[name] = evaluate_seeded(X, y_event, y_time, patients, cfg, n_seeds)
        null_scores[name] = evaluate_seeded(X, y_event, y_time, patients, cfg,
                                            n_null_seeds, permute=True)

    base = real_scores.get(baseline_name)

    # The dilution reference, keyed BY WIDTH: a block of width W is judged only
    # against a noise block of width W. Any real block must beat THIS, not
    # merely beat zero.
    noise_by_width: Dict[int, Tuple[str, np.ndarray]] = {}
    for nm, sc in real_scores.items():
        if _is_noise(nm):
            noise_by_width[feature_sets[nm].loc[common].shape[1]] = (nm, sc)

    # --- build table ---
    rows = []
    for name in real_scores:
        rmean, rlo, rhi = _ci(real_scores[name])
        nmean, _, nhi = _ci(null_scores[name])
        n_feat = feature_sets[name].loc[common].shape[1]

        # paired ΔC vs baseline (same seeds -> same folds)
        if base is not None and name != baseline_name and len(base) == len(real_scores[name]):
            delta = real_scores[name] - base
            dmean, dlo, dhi = _ci(delta)
        else:
            dmean = dlo = dhi = np.nan

        # paired ΔC vs the WIDTH-MATCHED noise block (the dilution-corrected test)
        is_noise_row = _is_noise(name)
        ref = noise_by_width.get(n_feat)
        if (ref is not None and not is_noise_row and name != baseline_name
                and len(ref[1]) == len(real_scores[name])):
            vn = real_scores[name] - ref[1]
            vn_mean, vn_lo, vn_hi = _ci(vn)
            vn_ref = ref[0]
        else:
            vn_mean = vn_lo = vn_hi = np.nan
            vn_ref = ''

        # verdict — judged against noise when available, else against baseline
        if np.isnan(rmean):
            verdict = "no result"
        elif rmean <= nhi:
            verdict = "= NOISE (below its own null)"
        elif name == baseline_name:
            verdict = "baseline"
        elif is_noise_row:
            verdict = "dilution control"
        elif not np.isnan(vn_mean):
            if vn_lo > 0:
                verdict = "REAL SIGNAL (beats noise)"
            elif vn_hi < 0:
                verdict = "worse than noise"
            else:
                verdict = "indistinguishable from noise"
        elif np.isnan(dmean):
            verdict = "above null"
        elif dlo > 0:
            verdict = "ADDS over baseline"
        elif dhi < 0:
            verdict = "HURTS vs baseline"
        else:
            verdict = "no gain over baseline"

        rows.append({
            'Feature Set': name,
            'Feats': n_feat,
            'C_mean': rmean, 'C_lo': rlo, 'C_hi': rhi,
            'Null_mean': nmean, 'Null_hi': nhi,
            'dC_mean': dmean, 'dC_lo': dlo, 'dC_hi': dhi,
            'vsNoise_mean': vn_mean, 'vsNoise_lo': vn_lo, 'vsNoise_hi': vn_hi,
            'NoiseRef': vn_ref,
            'Verdict': verdict,
        })

    summary = pd.DataFrame(rows).sort_values('C_mean', ascending=False).reset_index(drop=True)

    _banner("FINAL RESULTS")
    print(f"\n  {'Feature Set':<26s} {'Feats':>5s}  {'C-index (95% CI)':<22s} "
          f"{'Null':>6s}  {'dC vs baseline':<22s} {'dC vs NOISE':<22s} Verdict")
    print(f"  {'-' * 132}")
    for _, r in summary.iterrows():
        c_str = f"{r['C_mean']:.3f} [{r['C_lo']:.3f},{r['C_hi']:.3f}]"
        d_str = ("n/a" if np.isnan(r['dC_mean'])
                 else f"{r['dC_mean']:+.3f} [{r['dC_lo']:+.3f},{r['dC_hi']:+.3f}]")
        v_str = ("n/a" if np.isnan(r['vsNoise_mean'])
                 else f"{r['vsNoise_mean']:+.3f} [{r['vsNoise_lo']:+.3f},{r['vsNoise_hi']:+.3f}]")
        print(f"  {r['Feature Set']:<26s} {int(r['Feats']):>5d}  {c_str:<22s} "
              f"{r['Null_mean']:>6.3f}  {d_str:<22s} {v_str:<22s} {r['Verdict']}")

    print(f"\n  How to read this:")
    print(f"    - C-index (95% CI): real performance across seeds. Wide CI = unstable.")
    print(f"    - Null: this block's score with shuffled labels. Real must clear it.")
    print(f"    - dC vs baseline: paired gain over '{baseline_name}'.")
    print(f"    - dC vs NOISE: THE test. Same-width random block absorbs the dilution")
    print(f"      cost, so this isolates real signal. Signal only if its CI > 0.")

    # Experiment-tagged so a later run cannot silently clobber an earlier
    # experiment's results.
    save_path = output_dir / (f"survival_validation_{tag}.csv" if tag
                              else "survival_validation.csv")
    summary.to_csv(save_path, index=False)
    print(f"\n  Results saved: {save_path}")
    return summary


# ======================================================================
# Main
# ======================================================================

def main():
    parser = argparse.ArgumentParser(description="Survival RSF Evaluation")
    parser.add_argument('--pooling', default='mean_std',
                        choices=['mean', 'mean_std'],
                        help="How to pool cell->sample features")
    parser.add_argument('--n-seeds', type=int, default=20,
                        help="Number of CV seeds for the real score + CI")
    parser.add_argument('--n-null-seeds', type=int, default=10,
                        help="Number of CV seeds for the permutation null")
    parser.add_argument('--experiment', default='pe', choices=['pe', 'A', 'B', 'C'],
                        help="'pe' = PE vs celltype; 'A' = graph comparison "
                             "(readout fixed = neighbour matrix, graph is the variable); "
                             "'B' = two-view niches vs the best graph readout; "
                             "'C' = 5 abundance-corrected enrichment scalars vs "
                             "the 256-wide flattened readout")
    parser.add_argument('--exclude-normal', action='store_true',
                        help="Drop Normal mucosa acquisitions (the base paper does this)")
    parser.add_argument('--censor-doc', action='store_true',
                        help="Treat 'dead of other causes' as censored (cause-specific hazard)")
    parser.add_argument('--with-markers', action='store_true',
                        help="Add celltype-conditioned functional-marker NODE features "
                             "(8 feats: Ki67/PDL1/GranzymeB/PD1/ICOS/FoxP3/CD45RO/HLA-DR). "
                             "Appends a '+ Markers' variant to every Celltype-based set "
                             "and a 'Markers alone' set. Run: python src/marker_states.py")
    parser.add_argument('--marker-keep', default=None,
                        help="Comma-separated marker column names to KEEP (drops the rest). "
                             "Ranked strongest->weakest: cd4_foxp3, tcell_cd45ro, "
                             "tumor_mac_pdl1, tumor_ki67, cd8_granzymeb, cd4_icos, cd4_pd1, "
                             "apc_mac_hladr. Fewer, stronger markers = less dilution = "
                             "tighter CI. Only used with --with-markers.")
    parser.add_argument('--n-splits', type=int, default=None,
                        help="CV folds (default from config=10). Lower it when events are "
                             "scarce: --censor-doc leaves 16 patient events, so 10 folds "
                             "gives ~1.6 events/fold and the C-index becomes a coin flip.")
    args = parser.parse_args()

    cfg = get_config(pooling_method=args.pooling)
    if args.n_splits is not None:
        cfg.rsf_n_splits = args.n_splits
    cfg.print_summary()

    output_dir = cfg.processed_dir / "outputs" / "survival"
    output_dir.mkdir(parents=True, exist_ok=True)

    # --- Load survival labels ---
    survival = load_survival_labels(cfg, exclude_normal=args.exclude_normal,
                                    censor_doc=args.censor_doc)

    # --- Build feature sets ---
    sample_ids = survival['acquisition_id'].tolist()

    feature_sets = {}

    # 1. Celltype proportions (spatsurv baseline) — used by both experiments
    try:
        celltype_prop = load_celltype_proportions(cfg)
        feature_sets['Celltype Proportions'] = celltype_prop
    except Exception as e:
        print(f"  [WARNING] Could not load celltype props: {e}")

    if args.experiment == 'A':
        # ==============================================================
        # EXPERIMENT A — does the GRAPH matter?
        #   readout FIXED = row-normalised neighbour matrix (validated
        #     against the base paper's data/raw/k10/*.npy)
        #   cohort  FIXED
        #   variable = the graph only
        # Every config is the same width (16 + 256 = 272), so dilution is
        # identical across rows and the comparison is fair by construction.
        # ==============================================================
        _banner("EXPERIMENT A: graph comparison (readout fixed)")
        ct = feature_sets.get('Celltype Proportions')
        if ct is None:
            print("  [ERROR] Celltype proportions required for Experiment A")
            return

        for gname in ['knn10', 'radius20', 'radius50', 'delaunay']:
            nbr = load_neighbor_features(cfg, gname)
            if nbr.empty:
                continue
            combo = pd.concat([ct, nbr], axis=1, join='inner')
            if not combo.empty:
                feature_sets[f'Celltype + {gname}'] = combo

        # The dilution control: same width as a neighbour block (256 random cols).
        # Any graph must beat THIS row to count as real signal.
        n_noise = 256
        noise = make_noise_block(ct.index, n_noise, seed=42)
        feature_sets[f'Celltype + Noise({n_noise})'] = pd.concat(
            [ct, noise], axis=1, join='inner')

    elif args.experiment == 'B':
        # ==============================================================
        # EXPERIMENT B — is two-view propagation a better READOUT?
        #   graph   FIXED = Delaunay (the winner of Experiment A)
        #   variable = the readout: neighbour matrix vs two-view niches
        #
        # Two-view is only 12 features, so it is NOT diluted the way a 256-wide
        # block is. Its noise control is matched to ITS width (12), not 256 —
        # comparing a 12-wide block against a 256-wide noise block would be
        # rigged in its favour.
        # ==============================================================
        _banner("EXPERIMENT B: two-view niches vs neighbour matrix")
        ct = feature_sets.get('Celltype Proportions')
        if ct is None:
            print("  [ERROR] Celltype proportions required for Experiment B")
            return

        tv_path = cfg.processed_dir / "neighbor_features" / "two_view_niches.parquet"
        if not tv_path.exists():
            print(f"  [ERROR] {tv_path.name} missing — run: python src/two_view.py")
            return
        tv = pd.read_parquet(tv_path)
        print(f"  [two-view ] {tv.shape[0]} samples x {tv.shape[1]} niche proportions")

        feature_sets['Two-view alone'] = tv
        feature_sets['Celltype + TwoView'] = pd.concat([ct, tv], axis=1, join='inner')

        # width-matched dilution control for the 12-feature block
        noise12 = make_noise_block(ct.index, tv.shape[1], seed=42)
        feature_sets[f'Celltype + Noise({tv.shape[1]})'] = pd.concat(
            [ct, noise12], axis=1, join='inner')

        # reference: the best neighbour-matrix readout on the same graph
        nbr = load_neighbor_features(cfg, 'delaunay')
        if not nbr.empty:
            feature_sets['Celltype + delaunay'] = pd.concat([ct, nbr], axis=1, join='inner')
            feature_sets['Celltype + delaunay + TwoView'] = pd.concat(
                [ct, nbr, tv], axis=1, join='inner')

    elif args.experiment == 'C':
        # ==============================================================
        # EXPERIMENT C — 5 abundance-corrected scalars vs 256 flattened.
        #   graph  FIXED = Delaunay (Experiment A winner)
        #   cohort FIXED
        #   variable = the readout width/correction
        #
        # The 256 P[i][j] columns are confounded with global composition, which
        # the celltype baseline ALREADY supplies -> the block is largely a noisy
        # copy of the baseline (dC = +0.012). Dividing by p_j removes the shared
        # part, leaving only what abundance cannot explain. If the spatial signal
        # is real but diluted, 5 corrected features should recover more of it
        # than 256 raw ones. At ~27 patient events, EPV goes 0.099 -> 5.4.
        #
        # Width-matched noise controls — Noise(5) and Noise(256) — are generated
        # automatically in run_validation().
        # ==============================================================
        _banner("EXPERIMENT C: 5 enrichment scalars vs 256 flattened")
        ct = feature_sets.get('Celltype Proportions')
        if ct is None:
            print("  [ERROR] Celltype proportions required for Experiment C")
            return

        enr_path = cfg.processed_dir / "neighbor_features" / "enrichment.parquet"
        if not enr_path.exists():
            print(f"  [ERROR] {enr_path.name} missing — "
                  f"run: python src/enrichment_features.py")
            return
        enr = pd.read_parquet(enr_path)
        print(f"  [enrichment] {enr.shape[0]} samples x {enr.shape[1]} feats: "
              f"{list(enr.columns)}")

        feature_sets['Enrichment alone'] = enr
        feature_sets['Celltype + Enrichment'] = pd.concat(
            [ct, enr], axis=1, join='inner')

        # reference: the base paper's flattened readout on the same graph
        nbr = load_neighbor_features(cfg, 'delaunay')
        if not nbr.empty:
            feature_sets['Celltype + delaunay(256)'] = pd.concat(
                [ct, nbr], axis=1, join='inner')

    else:
        # ==============================================================
        # Default: the PE experiment
        # ==============================================================
        pe_features = pool_pe_features(cfg, sample_ids, pooling=args.pooling)
        if not pe_features.empty:
            feature_sets['PE (Laplacian)'] = pe_features

        cell_count = load_cell_count(cfg, sample_ids)
        if not cell_count.empty:
            feature_sets['Cell Count'] = cell_count

        raw_features = pool_raw_coord_features(cfg, sample_ids)
        if not raw_features.empty:
            feature_sets['Raw Coordinates'] = raw_features

        if 'Celltype Proportions' in feature_sets and 'PE (Laplacian)' in feature_sets:
            combined = pd.concat(
                [feature_sets['Celltype Proportions'], feature_sets['PE (Laplacian)']],
                axis=1, join='inner',
            )
            if not combined.empty:
                feature_sets['PE + Celltype'] = combined

        if 'PE (Laplacian)' in feature_sets and 'Raw Coordinates' in feature_sets:
            combined = pd.concat(
                [feature_sets['PE (Laplacian)'], feature_sets['Raw Coordinates']],
                axis=1, join='inner',
            )
            if not combined.empty:
                feature_sets['PE + Raw Coords'] = combined

    # --- node features: celltype-conditioned functional markers ---
    # Appends a '+ Markers' variant to every Celltype-based set (baseline,
    # delaunay, delaunay(256), Enrichment, ...) so each gets a paired ΔC vs its
    # own no-marker version. The width-matched noise control is added
    # automatically in run_validation(), so the honest verdict still applies.
    if args.with_markers:
        _banner("ADDING NODE FEATURES (celltype-conditioned functional markers)")
        markers = load_marker_states(cfg)
        if not markers.empty and args.marker_keep:
            keep = [c.strip() for c in args.marker_keep.split(',')]
            missing = [c for c in keep if c not in markers.columns]
            if missing:
                print(f"  [WARNING] --marker-keep names not found (ignored): {missing}")
            markers = markers[[c for c in keep if c in markers.columns]]
            print(f"  [marker-keep] reduced to {markers.shape[1]} markers: {list(markers.columns)}")
        if not markers.empty:
            base_names = [n for n in list(feature_sets)
                          if n.startswith('Celltype')
                          and 'Noise' not in n and 'Markers' not in n]
            for name in base_names:
                combo = pd.concat([feature_sets[name], markers], axis=1, join='inner')
                if not combo.empty:
                    feature_sets[f'{name} + Markers'] = combo
                    print(f"  [+markers] {name} + Markers -> {combo.shape[1]} feats")
            feature_sets['Markers alone'] = markers

    print(f"\n  Feature sets to evaluate: {list(feature_sets.keys())}")

    if len(feature_sets) < 1:
        print("  [ERROR] No feature sets available — run encoding pipeline first!")
        return

    # --- Run seeded validation (real vs null vs baseline) ---
    tag = args.experiment
    if args.exclude_normal:
        tag += '_exnorm'
    if args.censor_doc:
        tag += '_dod'
    if args.with_markers:
        n_mk = len(args.marker_keep.split(',')) if args.marker_keep else 8
        tag += f'_mk{n_mk}'

    run_validation(
        survival, feature_sets, cfg, output_dir,
        n_seeds=args.n_seeds, n_null_seeds=args.n_null_seeds,
        tag=tag,
    )


if __name__ == "__main__":
    main()
