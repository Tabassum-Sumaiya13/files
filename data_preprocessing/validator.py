"""
validator.py — Checks a raw external dataset against schema.py BEFORE any
processing happens. Nothing in `processed/` is written until this passes
(or the user explicitly forces processing despite FAILs).

Usage (normally via run_ingest.py, not directly):

    from validator import validate_dataset
    report = validate_dataset(adapter_config_module)
    print(report.to_markdown())
    report.ready   # True if there are zero FAIL checks
"""
from pathlib import Path

import pandas as pd

import schema
from report import ValidationReport


def _read(path) -> pd.DataFrame:
    path = Path(path)
    if path.suffix == ".csv":
        return pd.read_csv(path)
    if path.suffix == ".parquet":
        return pd.read_parquet(path)
    if path.suffix in (".pkl", ".pickle"):
        return pd.read_pickle(path)
    raise ValueError(f"Unsupported file type: {path.suffix} ({path})")


def validate_dataset(cfg) -> ValidationReport:
    """cfg is an adapter_config module — see datasets/_template/adapter_config.py."""
    report = ValidationReport(dataset_name=cfg.DATASET_NAME)
    print(f"\n{'=' * 70}\n  VALIDATING: {cfg.DATASET_NAME}\n{'=' * 70}")

    # 1. Files exist -----------------------------------------------------
    paths = {
        "locations": cfg.LOCATIONS_PATH,
        "expression": cfg.EXPRESSION_PATH,
        "metadata": cfg.METADATA_PATH,
    }
    for label, path in paths.items():
        p = Path(path)
        if p.exists():
            report.add(f"file:{label}", "PASS", str(p))
        else:
            report.add(f"file:{label}", "FAIL", f"not found: {p}")

    if report.n_fail:
        report.add("stopped_early", "FAIL", "cannot check columns — one or more raw files are missing")
        return report

    locations = _read(paths["locations"]).rename(columns=cfg.LOCATIONS_COLUMN_MAP)
    expression = _read(paths["expression"]).rename(columns=cfg.EXPRESSION_COLUMN_MAP)
    metadata = _read(paths["metadata"]).rename(columns=cfg.METADATA_COLUMN_MAP)

    # 2. Required columns present after mapping --------------------------
    for label, df, required in [
        ("locations", locations, schema.LOCATIONS_REQUIRED),
        ("expression", expression, schema.EXPRESSION_REQUIRED),
        ("metadata", metadata, schema.METADATA_REQUIRED),
    ]:
        missing = [c for c in required if c not in df.columns]
        if missing:
            report.add(f"columns:{label}", "FAIL",
                       f"missing after column-map rename: {missing} "
                       f"(fix {label.upper()}_COLUMN_MAP in adapter_config.py)")
        else:
            report.add(f"columns:{label}", "PASS", f"all {len(required)} canonical columns present")

    if report.n_fail:
        report.add("stopped_early", "FAIL", "cannot run further checks until required columns exist")
        return report

    locations[schema.ACQ_COL] = locations[schema.ACQ_COL].astype(str)
    locations[schema.CELL_COL] = locations[schema.CELL_COL].astype(str)
    expression[schema.ACQ_COL] = expression[schema.ACQ_COL].astype(str)
    expression[schema.CELL_COL] = expression[schema.CELL_COL].astype(str)
    metadata[schema.ACQ_COL] = metadata[schema.ACQ_COL].astype(str)

    # 3. dtypes ------------------------------------------------------------
    for col in (schema.X_COL, schema.Y_COL):
        if pd.api.types.is_numeric_dtype(locations[col]):
            report.add(f"dtype:{col}", "PASS", "numeric")
        else:
            report.add(f"dtype:{col}", "FAIL", f"{col} is not numeric — check for stray text/units in the raw column")

    # Survival is OPTIONAL. If present, sanity-check its dtype; if absent, WARN
    # (the cohort still ingests and can be feature-verified — only the separate,
    # per-dataset survival check is unavailable).
    has_survival = all(c in metadata.columns for c in schema.SURVIVAL_COLS)
    if not has_survival:
        missing_surv = [c for c in schema.SURVIVAL_COLS if c not in metadata.columns]
        report.add("survival:present", "WARN",
                   f"no survival columns {missing_surv} — cohort ingests for feature "
                   f"work, but the survival downstream check won't run "
                   f"(fill METADATA_COLUMN_MAP if this cohort has survival data)")
    else:
        report.add("survival:present", "PASS", "survival_day + survival_status present")
        if pd.api.types.is_numeric_dtype(metadata[schema.SURVIVAL_TIME_COL]):
            report.add("dtype:survival_day", "PASS", "numeric")
        else:
            report.add("dtype:survival_day", "FAIL", "survival_day is not numeric")

        status_vals = set(pd.to_numeric(metadata[schema.SURVIVAL_STATUS_COL], errors="coerce").dropna().unique().tolist())
        if status_vals <= {0.0, 1.0}:
            report.add("dtype:survival_status", "PASS", "binary 0/1")
        else:
            report.add("dtype:survival_status", "WARN",
                       f"unexpected values {sorted(status_vals)} (expected exactly {{0, 1}} — "
                       f"0=censored/alive, 1=event/dead)")

    # 4. Marker columns present --------------------------------------------
    marker_cols = cfg.marker_columns(expression)
    if len(marker_cols) == 0:
        report.add("markers:present", "FAIL", "no numeric marker columns detected — check marker_columns() in adapter_config.py")
        return report
    report.add("markers:present", "PASS", f"{len(marker_cols)} marker columns: {marker_cols[:8]}{'...' if len(marker_cols) > 8 else ''}")

    # 5. ID consistency across the 3 tables ---------------------------------
    loc_ids = set(locations[schema.ACQ_COL])
    expr_ids = set(expression[schema.ACQ_COL])
    meta_ids = set(metadata[schema.ACQ_COL])
    common = loc_ids & expr_ids & meta_ids
    report.add(
        "ids:overlap",
        "PASS" if common else "FAIL",
        f"{len(common)} acquisition_ids common to all 3 tables "
        f"(locations={len(loc_ids)}, expression={len(expr_ids)}, metadata={len(meta_ids)})",
    )
    orphans = (loc_ids | expr_ids) - meta_ids
    if orphans:
        report.add("ids:no_metadata", "WARN",
                   f"{len(orphans)} acquisition_ids have cell data but no metadata row "
                   f"— they will be dropped at the merge_metadata step: {sorted(orphans)[:5]}")

    # 6. Patient grouping present & multi-patient ---------------------------
    if schema.PATIENT_COL in metadata.columns:
        n_patients = metadata[schema.PATIENT_COL].nunique()
        n_samples = metadata[schema.ACQ_COL].nunique()
        report.add(
            "grouping:patient_id",
            "PASS" if n_patients > 1 else "FAIL",
            f"{n_patients} unique patients across {n_samples} samples "
            f"(needed for patient-grouped CV — a single patient_id would leak across folds)",
        )
    else:
        report.add("grouping:patient_id", "FAIL", "patient_id column missing from metadata")

    # 7. Sample size distribution -------------------------------------------
    min_cells = getattr(cfg, "MIN_CELLS_PER_SAMPLE", schema.MIN_CELLS_PER_SAMPLE_DEFAULT)
    sizes = locations.groupby(schema.ACQ_COL).size()
    n_small = int((sizes < min_cells).sum())
    report.add(
        "sample_size",
        "WARN" if n_small > 0 else "PASS",
        f"{n_small}/{len(sizes)} samples below MIN_CELLS_PER_SAMPLE={min_cells} "
        f"(median cells/sample={sizes.median():.0f}) — small samples are dropped during processing",
    )

    # 8. Cell-type taxonomy mapping coverage --------------------------------
    native_types = set(locations[schema.CLUSTER_LABEL_COL].astype(str).unique())
    celltype_map = getattr(cfg, "CELLTYPE_MAP", {})
    mapped = {t for t in native_types if t in celltype_map}
    unmapped = native_types - mapped
    bad_lineage = {celltype_map[t] for t in mapped} - set(schema.LINEAGES)
    if bad_lineage:
        report.add("celltype:mapping", "FAIL",
                   f"CELLTYPE_MAP contains lineage value(s) not in {schema.LINEAGES}: {bad_lineage}")
    elif not mapped:
        report.add("celltype:mapping", "FAIL", "CELLTYPE_MAP is empty — no native cell type maps to a lineage")
    elif unmapped:
        n_unmapped_cells = int(locations[schema.CLUSTER_LABEL_COL].astype(str).isin(unmapped).sum())
        report.add("celltype:mapping", "WARN",
                   f"{len(mapped)}/{len(native_types)} native types mapped "
                   f"({n_unmapped_cells:,}/{len(locations):,} cells unmapped) — "
                   f"unmapped types dropped at processing: {sorted(unmapped)[:10]}")
    else:
        report.add("celltype:mapping", "PASS", f"all {len(native_types)} native types mapped to a lineage")

    # 9. Marker panel overlap with feature-critical markers ------------------
    have_lineage = [m for m in schema.LINEAGE_VALIDATION_MARKERS if m in marker_cols]
    report.add(
        "markers:lineage_validation",
        "PASS" if len(have_lineage) >= 4 else "WARN",
        f"{len(have_lineage)}/{len(schema.LINEAGE_VALIDATION_MARKERS)} canonical lineage markers present: {have_lineage} "
        f"(used to sanity-check CELLTYPE_MAP against expression, not required to proceed)",
    )
    have_node = {feat: mk for feat, (mk, _) in schema.NODE_MARKER_FEATURES.items() if mk in marker_cols}
    report.add(
        "markers:node_features",
        "PASS" if len(have_node) > 0 else "WARN",
        f"{len(have_node)}/{len(schema.NODE_MARKER_FEATURES)} node-marker features reproducible on this cohort: "
        f"{list(have_node) if have_node else '(none)'}",
    )
    have_recommended = [m for m in schema.RECOMMENDED_MARKER_SET if m in marker_cols]
    report.add(
        "markers:recommended_pair",
        "PASS" if len(have_recommended) == len(schema.RECOMMENDED_MARKER_SET) else "WARN",
        f"this project's best config needs {schema.RECOMMENDED_MARKER_SET} (C=0.733, RESULT_REPORT.md Table 5); "
        f"present here: {have_recommended if have_recommended else '(none)'}",
    )

    # 10. Missingness ---------------------------------------------------------
    n_null_xy = int(locations[[schema.X_COL, schema.Y_COL]].isnull().sum().sum())
    report.add("missing:coords", "PASS" if n_null_xy == 0 else "WARN", f"{n_null_xy} null X/Y values (rows dropped at processing)")
    n_null_markers = int(expression[marker_cols].isnull().sum().sum())
    report.add("missing:markers", "PASS" if n_null_markers == 0 else "WARN", f"{n_null_markers} null marker values (rows dropped at processing)")

    return report
