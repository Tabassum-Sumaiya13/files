"""
processor.py — Transforms a VALIDATED external dataset into this project's
canonical processed format: one parquet per sample (same shape as
spatial_positional_encoding/data/processed/samples/) plus a manifest, and a
full change-log of every step performed so the resulting cohort is
auditable later.

Usage (normally via run_ingest.py, not directly):

    from processor import process_dataset
    from report import ChangeLog
    manifest = process_dataset(adapter_config_module, out_dir, ChangeLog("my_dataset"))
"""
from pathlib import Path

import numpy as np
import pandas as pd

import schema
from report import ChangeLog


def _read(path) -> pd.DataFrame:
    path = Path(path)
    if path.suffix == ".csv":
        return pd.read_csv(path)
    if path.suffix == ".parquet":
        return pd.read_parquet(path)
    if path.suffix in (".pkl", ".pickle"):
        return pd.read_pickle(path)
    raise ValueError(f"Unsupported file type: {path.suffix} ({path})")


def process_dataset(cfg, out_dir, log: ChangeLog = None) -> pd.DataFrame:
    log = log or ChangeLog(cfg.DATASET_NAME)
    out_dir = Path(out_dir)
    samples_dir = out_dir / "samples"
    samples_dir.mkdir(parents=True, exist_ok=True)

    min_cells = getattr(cfg, "MIN_CELLS_PER_SAMPLE", schema.MIN_CELLS_PER_SAMPLE_DEFAULT)
    celltype_map = getattr(cfg, "CELLTYPE_MAP", {})

    # --- Load + rename to canonical columns -----------------------------
    locations = _read(cfg.LOCATIONS_PATH).rename(columns=cfg.LOCATIONS_COLUMN_MAP)
    expression = _read(cfg.EXPRESSION_PATH).rename(columns=cfg.EXPRESSION_COLUMN_MAP)
    metadata = _read(cfg.METADATA_PATH).rename(columns=cfg.METADATA_COLUMN_MAP)
    marker_cols = cfg.marker_columns(expression)
    log.step("load", f"locations={len(locations):,} rows, expression={len(expression):,} rows, "
                      f"metadata={len(metadata):,} rows, {len(marker_cols)} marker columns detected")

    locations[schema.ACQ_COL] = locations[schema.ACQ_COL].astype(str)
    locations[schema.CELL_COL] = locations[schema.CELL_COL].astype(str)
    expression[schema.ACQ_COL] = expression[schema.ACQ_COL].astype(str)
    expression[schema.CELL_COL] = expression[schema.CELL_COL].astype(str)
    metadata[schema.ACQ_COL] = metadata[schema.ACQ_COL].astype(str)

    # --- Dedup ------------------------------------------------------------
    n0 = len(locations)
    locations = locations.drop_duplicates(subset=[schema.ACQ_COL, schema.CELL_COL])
    log.step("dedup:locations", f"{n0:,} -> {len(locations):,} rows ({n0 - len(locations):,} duplicate (sample, cell) pairs dropped)")

    n0 = len(expression)
    expression = expression.drop_duplicates(subset=[schema.ACQ_COL, schema.CELL_COL])
    log.step("dedup:expression", f"{n0:,} -> {len(expression):,} rows ({n0 - len(expression):,} duplicate (sample, cell) pairs dropped)")

    # --- Drop missing coords / marker values -------------------------------
    n0 = len(locations)
    locations = locations.dropna(subset=[schema.X_COL, schema.Y_COL])
    log.step("drop_missing:coords", f"{n0:,} -> {len(locations):,} rows ({n0 - len(locations):,} dropped for missing X/Y)")

    n0 = len(expression)
    expression = expression.dropna(subset=marker_cols)
    log.step("drop_missing:markers", f"{n0:,} -> {len(expression):,} rows ({n0 - len(expression):,} dropped for missing marker values)")

    locations[schema.X_COL] = locations[schema.X_COL].astype(float)
    locations[schema.Y_COL] = locations[schema.Y_COL].astype(float)

    # --- arcsinh normalise marker intensities (if not already normalised) --
    if getattr(cfg, "APPLY_ARCSINH", True):
        cofactor = getattr(cfg, "ARCSINH_COFACTOR", 5.0)
        expression[marker_cols] = np.arcsinh(expression[marker_cols].astype(float) / cofactor)
        log.step("normalise_markers", f"arcsinh(x / {cofactor}) applied to {len(marker_cols)} marker columns")
    else:
        log.step("normalise_markers", "skipped — adapter_config.APPLY_ARCSINH=False, data already normalised")

    # --- Merge locations + expression on (acquisition_id, cell_id) --------
    expr_for_merge = expression[[schema.ACQ_COL, schema.CELL_COL] + marker_cols]
    n_loc = len(locations)
    merged = locations.merge(expr_for_merge, on=[schema.ACQ_COL, schema.CELL_COL], how="inner")
    log.step("merge_locations_expression",
             f"{n_loc:,} location rows -> {len(merged):,} merged rows "
             f"({n_loc - len(merged):,} location rows had no matching expression row and were dropped)")

    # --- Filter samples below the cell-count floor -------------------------
    sizes = merged.groupby(schema.ACQ_COL).size()
    small = sizes[sizes < min_cells].index.tolist()
    n_before = merged[schema.ACQ_COL].nunique()
    merged = merged[~merged[schema.ACQ_COL].isin(small)]
    log.step("filter_small_samples",
             f"{n_before} -> {merged[schema.ACQ_COL].nunique()} samples "
             f"(dropped {len(small)} sample(s) below MIN_CELLS_PER_SAMPLE={min_cells}: {small[:10]})")

    # --- Cell-type -> lineage harmonisation --------------------------------
    merged["lineage"] = merged[schema.CLUSTER_LABEL_COL].astype(str).map(celltype_map)
    n_unmapped = int(merged["lineage"].isnull().sum())
    if n_unmapped:
        dropped_types = sorted(merged.loc[merged["lineage"].isnull(), schema.CLUSTER_LABEL_COL].astype(str).unique())
        n_before = len(merged)
        merged = merged.dropna(subset=["lineage"])
        log.step("celltype_to_lineage",
                 f"{n_unmapped:,} cells ({n_unmapped / n_before:.1%}) had a cell type with no lineage in the "
                 f"cell-type registry and were dropped: {dropped_types[:10]} — to keep them, give them a "
                 f"lineage in celltype_registry.csv (see registry.py); the reason each is currently excluded "
                 f"is recorded in that file's `notes` column",
                 level="warn")
    else:
        log.step("celltype_to_lineage",
                 f"all {len(merged):,} cells mapped to a lineage: {merged['lineage'].value_counts().to_dict()}")

    # --- Per-sample z-score of coordinates ----------------------------------
    merged = merged.reset_index(drop=True)
    for sid, grp in merged.groupby(schema.ACQ_COL):
        mask = merged[schema.ACQ_COL] == sid
        for c in (schema.X_COL, schema.Y_COL):
            std = grp[c].std()
            std = std if std > 1e-8 else 1e-8
            merged.loc[mask, c] = (grp[c] - grp[c].mean()) / std
    log.step("normalise_coords", "per-sample z-score applied to X, Y (matches spatial_positional_encoding/src/preprocess.py Step 5)")

    # --- Attach survival metadata at the sample level -----------------------
    meta_cols = [c for c in [schema.ACQ_COL, schema.PATIENT_COL, schema.SURVIVAL_TIME_COL, schema.SURVIVAL_STATUS_COL]
                 if c in metadata.columns]
    metadata_small = metadata[meta_cols].drop_duplicates(subset=[schema.ACQ_COL])
    n_before = merged[schema.ACQ_COL].nunique()
    merged = merged.merge(metadata_small, on=schema.ACQ_COL, how="inner")
    n_after = merged[schema.ACQ_COL].nunique()
    # Survival is optional: it's carried through only if the metadata had it.
    has_surv = all(c in merged.columns for c in schema.SURVIVAL_COLS)
    log.step("merge_metadata",
             f"{n_before} -> {n_after} samples kept a metadata row "
             f"({n_before - n_after} sample(s) had no metadata row and were dropped); "
             f"survival {'present' if has_surv else 'ABSENT (survival-less cohort)'}")

    # --- Export one parquet per sample + a manifest -------------------------
    manifest_rows = []
    for sid, g in merged.groupby(schema.ACQ_COL):
        safe = str(sid).replace("/", "_").replace("\\", "_")
        path = samples_dir / f"sample_{safe}.parquet"
        g.to_parquet(path, index=False)
        manifest_rows.append({
            schema.ACQ_COL: sid,
            schema.PATIENT_COL: g[schema.PATIENT_COL].iloc[0],
            "n_cells": len(g),
            "n_cell_types": g[schema.CLUSTER_LABEL_COL].nunique(),
            schema.SURVIVAL_TIME_COL: g[schema.SURVIVAL_TIME_COL].iloc[0] if has_surv else np.nan,
            schema.SURVIVAL_STATUS_COL: g[schema.SURVIVAL_STATUS_COL].iloc[0] if has_surv else np.nan,
            "file_path": str(path),
        })
    manifest = pd.DataFrame(manifest_rows)
    manifest.to_parquet(out_dir / "manifest.parquet", index=False)
    (out_dir / "marker_columns.txt").write_text("\n".join(marker_cols), encoding="utf-8")

    log.step("export", f"{len(manifest)} sample parquet files -> {samples_dir}, manifest -> {out_dir / 'manifest.parquet'}")
    if has_surv:
        surv_note = (f"{int(manifest[schema.SURVIVAL_STATUS_COL].sum())} events "
                     f"({manifest[schema.SURVIVAL_STATUS_COL].mean():.1%} event rate)")
    else:
        surv_note = "no survival (feature-verification cohort only)"
    log.step(
        "final_cohort",
        f"{int(manifest['n_cells'].sum()):,} cells, {len(manifest)} samples, "
        f"{manifest[schema.PATIENT_COL].nunique()} patients, {surv_note}",
    )

    return manifest
