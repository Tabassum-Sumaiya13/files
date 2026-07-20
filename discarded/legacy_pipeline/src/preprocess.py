"""
src/preprocess.py — Data loading, validation, merge, QC, normalisation, export.

Reads the three raw CSV files, standardises column names, merges locations
with expression, applies per-sample coordinate normalisation, and exports
one parquet per sample.

Every step prints verbose debug output so you always know what happened.
"""

import pandas as pd
import numpy as np
from pathlib import Path
from typing import List, Tuple, Optional
from tqdm import tqdm
import warnings


def _banner(title: str):
    """Print a clearly visible section header."""
    print(f"\n{'=' * 70}")
    print(f"  {title}")
    print(f"{'=' * 70}")


def load_marker_names(cfg) -> List[str]:
    """Load the canonical marker name list from marker_names.csv."""
    path = cfg.marker_names_path
    print(f"\n[load_marker_names] Reading: {path}")
    marker_df = pd.read_csv(path, header=None)
    markers = marker_df[0].tolist()
    print(f"  => {len(markers)} markers: {markers[:5]} ... {markers[-3:]}")
    return markers


def load_qc_sample_ids(cfg) -> List[str]:
    """Load the list of 307 QC-passing acquisition IDs."""
    path = cfg.qc_ids_path
    print(f"\n[load_qc_sample_ids] Reading: {path}")
    df = pd.read_csv(path, header=None)
    ids = df[0].astype(str).tolist()
    print(f"  => {len(ids)} QC-passing sample IDs loaded")
    print(f"  => First 3: {ids[:3]}")
    return ids


def load_raw_data(cfg) -> Tuple[pd.DataFrame, pd.DataFrame, List[str]]:
    """
    Load and standardise the two raw CSV files + marker names.

    Returns
    -------
    locations : DataFrame  [acquisition_id, cell_id, X, Y, cluster_id, cluster_label]
    expression : DataFrame [acquisition_id, cell_id, cluster_id, cluster_label, <marker_cols>]
    marker_cols : list[str]
    """
    marker_cols = load_marker_names(cfg)

    # --- Locations ---
    _banner("LOADING: cell_locations_and_labels.csv")
    print(f"  Path: {cfg.locations_path}")
    locations = pd.read_csv(cfg.locations_path)
    print(f"  Raw shape: {locations.shape}")
    print(f"  Raw columns: {list(locations.columns)}")
    print(f"  Sample dtypes:\n{locations.dtypes.to_string()}")

    # Standardise column names
    locations = locations.rename(columns={
        cfg.LOC_ACQ_COL:  cfg.ACQ_COL,
        cfg.LOC_CELL_COL: cfg.CELL_COL,
        cfg.LOC_X_COL:    cfg.X_COL,
        cfg.LOC_Y_COL:    cfg.Y_COL,
        cfg.LOC_CLUSTER_ID_COL:    cfg.CLUSTER_ID_COL,
        cfg.LOC_CLUSTER_LABEL_COL: cfg.CLUSTER_LABEL_COL,
    })
    locations[cfg.ACQ_COL] = locations[cfg.ACQ_COL].astype(str)
    locations[cfg.CELL_COL] = locations[cfg.CELL_COL].astype(str)
    locations[cfg.CLUSTER_LABEL_COL] = locations[cfg.CLUSTER_LABEL_COL].astype(str)

    # Keep only needed columns
    keep = [cfg.ACQ_COL, cfg.CELL_COL, cfg.X_COL, cfg.Y_COL,
            cfg.CLUSTER_ID_COL, cfg.CLUSTER_LABEL_COL]
    locations = locations[[c for c in keep if c in locations.columns]]

    print(f"\n  After rename: {locations.shape}")
    print(f"  Columns: {list(locations.columns)}")
    print(f"  Unique samples: {locations[cfg.ACQ_COL].nunique()}")
    print(f"  Head:\n{locations.head(3).to_string()}")

    # --- Expression ---
    _banner("LOADING: labeled_arcsinh_norm_data.csv")
    print(f"  Path: {cfg.expression_path}")
    expression = pd.read_csv(cfg.expression_path)
    print(f"  Raw shape: {expression.shape}")
    print(f"  Raw columns: {list(expression.columns)}")

    # Standardise
    expression = expression.rename(columns={
        cfg.EXPR_SAMPLE_COL:       cfg.ACQ_COL,
        cfg.EXPR_CELL_COL:         cfg.CELL_COL,
        cfg.EXPR_CLUSTER_COL:      cfg.CLUSTER_ID_COL,
        cfg.EXPR_CLUSTER_LABEL_COL: cfg.CLUSTER_LABEL_COL,
    })
    expression[cfg.ACQ_COL] = expression[cfg.ACQ_COL].astype(str)
    expression[cfg.CELL_COL] = expression[cfg.CELL_COL].astype(str)

    # Drop junk columns
    for col in ['Unnamed: 0']:
        if col in expression.columns:
            expression = expression.drop(columns=[col])

    # Keep only: acq, cell, cluster_id, cluster_label, markers
    expr_keep = [cfg.ACQ_COL, cfg.CELL_COL, cfg.CLUSTER_ID_COL, cfg.CLUSTER_LABEL_COL] + marker_cols
    expr_keep = [c for c in expr_keep if c in expression.columns]
    expression = expression[expr_keep]

    print(f"\n  After rename: {expression.shape}")
    print(f"  Marker columns present: {len([c for c in marker_cols if c in expression.columns])}/{len(marker_cols)}")
    missing = [c for c in marker_cols if c not in expression.columns]
    if missing:
        print(f"  [WARNING] Missing markers: {missing}")
        marker_cols = [c for c in marker_cols if c in expression.columns]
    print(f"  Head:\n{expression.head(3).to_string()}")

    return locations, expression, marker_cols


def validate_and_clean(
    locations: pd.DataFrame,
    expression: pd.DataFrame,
    marker_cols: List[str],
    cfg,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Validate both DataFrames, drop duplicates, handle missing values."""
    _banner("STEP 1: VALIDATION & CLEANING")

    # 1a. Duplicates in locations
    n_dup = locations.duplicated().sum()
    if n_dup > 0:
        locations = locations.drop_duplicates(keep='first')
        print(f"  [CLEANED] Dropped {n_dup} duplicate rows from locations")
    else:
        print(f"  [OK] No duplicate rows in locations")

    # 1b. Duplicates in expression
    n_dup = expression.duplicated().sum()
    if n_dup > 0:
        expression = expression.drop_duplicates(keep='first')
        print(f"  [CLEANED] Dropped {n_dup} duplicate rows from expression")
    else:
        print(f"  [OK] No duplicate rows in expression")

    # 1c. Duplicate (acquisition_id, cell_id) pairs
    dup_pairs = locations.duplicated(subset=[cfg.ACQ_COL, cfg.CELL_COL]).sum()
    if dup_pairs > 0:
        locations = locations.drop_duplicates(subset=[cfg.ACQ_COL, cfg.CELL_COL], keep='first')
        print(f"  [CLEANED] Dropped {dup_pairs} duplicate (sample, cell) pairs from locations")
    else:
        print(f"  [OK] No duplicate (sample, cell) pairs in locations")

    dup_pairs = expression.duplicated(subset=[cfg.ACQ_COL, cfg.CELL_COL]).sum()
    if dup_pairs > 0:
        expression = expression.drop_duplicates(subset=[cfg.ACQ_COL, cfg.CELL_COL], keep='first')
        print(f"  [CLEANED] Dropped {dup_pairs} duplicate (sample, cell) pairs from expression")
    else:
        print(f"  [OK] No duplicate (sample, cell) pairs in expression")

    # 1d. Missing coordinates
    missing_xy = locations[[cfg.X_COL, cfg.Y_COL]].isnull().sum().sum()
    if missing_xy > 0:
        before = len(locations)
        locations = locations.dropna(subset=[cfg.X_COL, cfg.Y_COL])
        print(f"  [CLEANED] Dropped {before - len(locations)} rows with missing X/Y")
    else:
        print(f"  [OK] No missing coordinate values")

    # 1e. Missing marker values
    missing_markers = expression[marker_cols].isnull().sum().sum()
    if missing_markers > 0:
        before = len(expression)
        expression = expression.dropna(subset=marker_cols)
        print(f"  [CLEANED] Dropped {before - len(expression)} rows with missing marker values")
    else:
        print(f"  [OK] No missing marker values")

    # 1f. Type conversion
    locations[[cfg.X_COL, cfg.Y_COL]] = locations[[cfg.X_COL, cfg.Y_COL]].astype(float)

    print(f"\n  After cleaning:")
    print(f"    Locations:  {len(locations):,} rows × {len(locations.columns)} cols")
    print(f"    Expression: {len(expression):,} rows × {len(expression.columns)} cols")
    print(f"    Locations unique samples:  {locations[cfg.ACQ_COL].nunique()}")
    print(f"    Expression unique samples: {expression[cfg.ACQ_COL].nunique()}")

    return locations, expression


def filter_qc_samples(
    locations: pd.DataFrame,
    expression: pd.DataFrame,
    cfg,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Keep only the 307 QC-passing samples from the baseline paper."""
    _banner("STEP 2: QC SAMPLE FILTERING")

    if not cfg.use_qc_sample_list:
        print("  [SKIP] use_qc_sample_list is False — keeping all samples")
        return locations, expression

    qc_ids = load_qc_sample_ids(cfg)
    qc_set = set(qc_ids)

    before_loc = locations[cfg.ACQ_COL].nunique()
    before_expr = expression[cfg.ACQ_COL].nunique()

    locations = locations[locations[cfg.ACQ_COL].isin(qc_set)]
    expression = expression[expression[cfg.ACQ_COL].isin(qc_set)]

    after_loc = locations[cfg.ACQ_COL].nunique()
    after_expr = expression[cfg.ACQ_COL].nunique()

    print(f"  Locations:  {before_loc} -> {after_loc} samples")
    print(f"  Expression: {before_expr} -> {after_expr} samples")
    print(f"  Locations rows:  {len(locations):,}")
    print(f"  Expression rows: {len(expression):,}")

    return locations, expression


def merge_data(
    locations: pd.DataFrame,
    expression: pd.DataFrame,
    marker_cols: List[str],
    cfg,
) -> pd.DataFrame:
    """Inner-merge locations + expression on (acquisition_id, cell_id)."""
    _banner("STEP 3: MERGE LOCATIONS + EXPRESSION")

    print(f"  Locations:  {len(locations):,} rows")
    print(f"  Expression: {len(expression):,} rows")

    # The expression file has cluster_id and cluster_label too — drop them
    # before merge to avoid duplicates, we keep the ones from locations.
    expr_for_merge = expression.drop(
        columns=[c for c in [cfg.CLUSTER_ID_COL, cfg.CLUSTER_LABEL_COL]
                 if c in expression.columns],
        errors='ignore',
    )

    merged = locations.merge(expr_for_merge, on=[cfg.ACQ_COL, cfg.CELL_COL], how='inner')

    dropped = len(locations) - len(merged)
    print(f"\n  After merge: {len(merged):,} rows")
    print(f"  Dropped:     {dropped:,} ({dropped / max(len(locations), 1):.1%})")
    print(f"  Columns:     {list(merged.columns)}")
    print(f"  Unique samples: {merged[cfg.ACQ_COL].nunique()}")

    # Check for any sample with high drop rate
    sample_sizes_before = locations.groupby(cfg.ACQ_COL).size()
    sample_sizes_after = merged.groupby(cfg.ACQ_COL).size()
    for sid in sample_sizes_before.index:
        before = sample_sizes_before.get(sid, 0)
        after = sample_sizes_after.get(sid, 0)
        if before > 0 and (before - after) / before > 0.05:
            print(f"  [WARNING] {sid}: {before - after}/{before} cells lost ({(before - after) / before:.1%})")

    return merged


def filter_small_samples(merged: pd.DataFrame, cfg) -> pd.DataFrame:
    """Drop samples with fewer than min_cells_per_sample cells."""
    _banner("STEP 4: FILTER SMALL SAMPLES")

    sample_sizes = merged.groupby(cfg.ACQ_COL).size()
    before = len(sample_sizes)
    small = sample_sizes[sample_sizes < cfg.min_cells_per_sample]

    if len(small) > 0:
        print(f"  Dropping {len(small)} samples with < {cfg.min_cells_per_sample} cells:")
        for sid, n in small.items():
            print(f"    {sid}: {n} cells")
        merged = merged[merged[cfg.ACQ_COL].isin(sample_sizes[sample_sizes >= cfg.min_cells_per_sample].index)]
    else:
        print(f"  [OK] All {before} samples have >= {cfg.min_cells_per_sample} cells")

    print(f"\n  After filter: {merged[cfg.ACQ_COL].nunique()} samples, {len(merged):,} cells")
    return merged


def normalize_coordinates(merged: pd.DataFrame, cfg) -> pd.DataFrame:
    """Per-sample z-score normalisation of X, Y coordinates."""
    _banner("STEP 5: COORDINATE NORMALISATION")

    if not cfg.normalize_coords:
        print("  [SKIP] normalize_coords is False")
        return merged

    print("  Applying per-sample z-score to X, Y...")
    normalized = merged.copy()

    issues = []
    for sample_id, grp in normalized.groupby(cfg.ACQ_COL):
        mask = normalized[cfg.ACQ_COL] == sample_id

        x_mean, y_mean = grp[cfg.X_COL].mean(), grp[cfg.Y_COL].mean()
        x_std, y_std = grp[cfg.X_COL].std(), grp[cfg.Y_COL].std()

        if x_std == 0 or y_std == 0:
            issues.append(sample_id)

        x_std = max(x_std, 1e-8)
        y_std = max(y_std, 1e-8)

        normalized.loc[mask, cfg.X_COL] = (grp[cfg.X_COL] - x_mean) / x_std
        normalized.loc[mask, cfg.Y_COL] = (grp[cfg.Y_COL] - y_mean) / y_std

    if issues:
        print(f"  [WARNING] {len(issues)} samples had zero variance in X or Y: {issues[:5]}")
    else:
        print(f"  [OK] All samples normalised successfully")

    # Verify
    check = normalized.groupby(cfg.ACQ_COL)[[cfg.X_COL, cfg.Y_COL]].agg(['mean', 'std'])
    print(f"\n  Post-normalisation stats (first 3 samples):")
    print(f"{check.head(3).to_string()}")

    return normalized


def quality_control(merged: pd.DataFrame, marker_cols: List[str], cfg) -> dict:
    """Run quality control checks and print diagnostics."""
    _banner("STEP 6: QUALITY CONTROL")

    qc = {}

    # Sample size distribution
    sample_sizes = merged.groupby(cfg.ACQ_COL).size()
    qc['n_samples'] = len(sample_sizes)
    qc['n_cells'] = len(merged)
    qc['cells_per_sample_mean'] = sample_sizes.mean()
    qc['cells_per_sample_median'] = sample_sizes.median()
    qc['cells_per_sample_min'] = sample_sizes.min()
    qc['cells_per_sample_max'] = sample_sizes.max()

    print(f"  Samples: {qc['n_samples']}")
    print(f"  Total cells: {qc['n_cells']:,}")
    print(f"  Cells/sample: mean={qc['cells_per_sample_mean']:.0f}, "
          f"median={qc['cells_per_sample_median']:.0f}, "
          f"min={qc['cells_per_sample_min']}, max={qc['cells_per_sample_max']}")

    # Marker statistics
    marker_data = merged[marker_cols].values
    zero_rates = (marker_data == 0).mean(axis=0)
    high_zero = [(marker_cols[i], zero_rates[i]) for i in range(len(marker_cols)) if zero_rates[i] > 0.9]
    if high_zero:
        print(f"\n  [WARNING] {len(high_zero)} markers have >90% zeros:")
        for name, rate in high_zero:
            print(f"    {name}: {rate:.1%}")
    else:
        print(f"\n  [OK] No markers have >90% zeros")

    # Cell type distribution
    n_types = merged[cfg.CLUSTER_LABEL_COL].nunique()
    print(f"\n  Unique cell types: {n_types}")
    print(f"  Top 5 cell types:")
    top5 = merged[cfg.CLUSTER_LABEL_COL].value_counts().head(5)
    for ct, count in top5.items():
        print(f"    {ct}: {count:,} ({count / len(merged):.1%})")

    # NaN/inf check
    nan_count = merged[marker_cols].isnull().sum().sum()
    inf_count = (~np.isfinite(merged[marker_cols].values)).sum()
    print(f"\n  NaN in markers: {nan_count}")
    print(f"  Inf in markers: {inf_count}")
    if nan_count > 0 or inf_count > 0:
        print(f"  [ERROR] Data has NaN or Inf — fix before proceeding!")

    return qc


def export_samples(normalized: pd.DataFrame, marker_cols: List[str], cfg) -> pd.DataFrame:
    """Write one parquet per sample and a manifest."""
    _banner("STEP 7: EXPORT PER-SAMPLE PARQUET FILES")

    cfg.ensure_dirs()

    sample_info = []
    for sample_id, grp in tqdm(normalized.groupby(cfg.ACQ_COL), desc="Exporting"):
        safe_id = str(sample_id).replace('/', '_').replace('\\', '_')
        out_path = cfg.samples_dir / f"sample_{safe_id}.parquet"
        grp.to_parquet(out_path, compression='snappy', index=False)
        sample_info.append({
            cfg.ACQ_COL: sample_id,
            'n_cells': len(grp),
            'n_cell_types': grp[cfg.CLUSTER_LABEL_COL].nunique(),
            'file_path': str(out_path),
        })

    manifest = pd.DataFrame(sample_info)
    manifest.to_parquet(cfg.manifest_path, compression='snappy', index=False)

    print(f"\n  Exported {len(manifest)} samples to {cfg.samples_dir}")
    print(f"  Manifest saved to {cfg.manifest_path}")
    print(f"  Total cells: {manifest['n_cells'].sum():,}")
    print(f"  Head:\n{manifest.head(5).to_string()}")

    return manifest


def verify_output(manifest: pd.DataFrame, marker_cols: List[str], cfg):
    """Re-read exported files and verify no NaN/inf slipped through."""
    _banner("STEP 8: INDEPENDENT VERIFICATION")

    bad_files = []
    total_cells = 0
    for _, row in manifest.iterrows():
        df = pd.read_parquet(row['file_path'])
        total_cells += len(df)
        if df[marker_cols].isnull().any().any():
            bad_files.append(row[cfg.ACQ_COL])
        if not np.isfinite(df[[cfg.X_COL, cfg.Y_COL]].values).all():
            bad_files.append(row[cfg.ACQ_COL])

    print(f"  Checked {len(manifest)} files")
    print(f"  Total cells across all files: {total_cells:,}")
    if bad_files:
        print(f"  [ERROR] {len(set(bad_files))} files have NaN/inf: {sorted(set(bad_files))[:10]}")
    else:
        print(f"  [OK] All output files are clean (no NaN/inf)")


# ======================================================================
# Main entry point
# ======================================================================

def run_preprocessing(cfg) -> Tuple[pd.DataFrame, pd.DataFrame, List[str]]:
    """Run the full preprocessing pipeline.

    Returns
    -------
    normalized : DataFrame with all cleaned, normalised data
    manifest   : DataFrame mapping sample IDs to output files
    marker_cols: list of marker column names
    """
    _banner("PREPROCESSING PIPELINE START")
    cfg.ensure_dirs()

    # Load
    locations, expression, marker_cols = load_raw_data(cfg)

    # Debug subset
    if cfg.debug_n_samples is not None:
        _banner(f"DEBUG MODE: limiting to {cfg.debug_n_samples} samples")
        keep = locations[cfg.ACQ_COL].unique()[:cfg.debug_n_samples]
        locations = locations[locations[cfg.ACQ_COL].isin(keep)]
        expression = expression[expression[cfg.ACQ_COL].isin(keep)]
        print(f"  Kept {len(keep)} samples")

    # Step 1: Validate
    locations, expression = validate_and_clean(locations, expression, marker_cols, cfg)

    # Step 2: QC filter
    locations, expression = filter_qc_samples(locations, expression, cfg)

    # Step 3: Merge
    merged = merge_data(locations, expression, marker_cols, cfg)

    # Step 4: Filter small
    merged = filter_small_samples(merged, cfg)

    # Step 5: Normalise coords
    normalized = normalize_coordinates(merged, cfg)

    # Step 6: QC
    quality_control(normalized, marker_cols, cfg)

    # Step 7: Export
    manifest = export_samples(normalized, marker_cols, cfg)

    # Step 8: Verify
    verify_output(manifest, marker_cols, cfg)

    _banner("PREPROCESSING COMPLETE")
    print(f"  Final: {len(manifest)} samples, {manifest['n_cells'].sum():,} cells")
    print(f"  {len(marker_cols)} markers: {marker_cols[:5]}...")
    print(f"  Output: {cfg.samples_dir}")

    return normalized, manifest, marker_cols