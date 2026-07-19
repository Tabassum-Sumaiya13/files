"""
adapter_config.py — UPMC (this project's baseline head-and-neck CODEX cohort).

This re-ingests the ORIGINAL baseline dataset through the same
data_preprocessing path every external cohort uses, so UPMC and any new
cohort (CRC, …) land in one identical canonical format
(datasets/UPMC/processed/) that spatial_positional_encoding then reads.

Source files (already the project's cleaned UPMC tables):
    ../../../data/raw/dataset_info/cell_locations_and_labels.csv
    ../../../data/raw/dataset_info/labeled_arcsinh_norm_data.csv   (already arcsinh)
    ../../../data/raw/dataset_info/sample_metadata.csv             (has survival)

Run:
    python run_ingest.py --dataset UPMC
"""
from pathlib import Path

DATASET_NAME = "UPMC"

# ---------------------------------------------------------------------------
# 1. Raw files — point back at the project's existing UPMC tables under
#    <repo>/data/raw/dataset_info/. (datasets/UPMC/ -> repo root is ../../../..)
# ---------------------------------------------------------------------------
_REPO_ROOT = Path(__file__).resolve().parents[3]
_RAW = _REPO_ROOT / "data" / "raw" / "dataset_info"

LOCATIONS_PATH = _RAW / "cell_locations_and_labels.csv"
EXPRESSION_PATH = _RAW / "labeled_arcsinh_norm_data.csv"
METADATA_PATH = _RAW / "sample_metadata.csv"

# ---------------------------------------------------------------------------
# 2. Column renames -> canonical schema.
#    locations : ACQUISITION_ID/CELL_ID/CLUSTER_LABEL renamed; X/Y already canonical.
#    expression: sample_id->acquisition_id, cluster->cluster_id (so it is
#                excluded from markers); cell_id already canonical.
#    metadata  : acquisition_id / patient_id / survival_day / survival_status
#                are ALREADY the canonical names -> empty map.
# ---------------------------------------------------------------------------
LOCATIONS_COLUMN_MAP = {
    "ACQUISITION_ID": "acquisition_id",
    "CELL_ID": "cell_id",
    "CLUSTER_LABEL": "cluster_label",
    "CLUSTER_ID": "cluster_id",
}
EXPRESSION_COLUMN_MAP = {
    "sample_id": "acquisition_id",
    "cluster": "cluster_id",
}
METADATA_COLUMN_MAP = {
    # all four canonical metadata columns already present under their canonical
    # names (acquisition_id, patient_id, survival_day, survival_status)
}


# ---------------------------------------------------------------------------
# 3. Marker columns. The expression table carries an unnamed pandas index
#    ("Unnamed: 0") and the integer cluster/id columns, all numeric, which the
#    default "every numeric column" rule would misread as markers. Exclude the
#    id/label/index columns explicitly -> the 39 real protein markers.
# ---------------------------------------------------------------------------
_NON_MARKER_COLS = {
    "acquisition_id", "cell_id", "cluster_id", "cluster_label",
    "cluster", "sample_id",
}


def marker_columns(expression_df) -> list:
    import pandas as pd
    return [
        c for c in expression_df.columns
        if c not in _NON_MARKER_COLS
        and not str(c).startswith("Unnamed")
        and pd.api.types.is_numeric_dtype(expression_df[c])
    ]


# ---------------------------------------------------------------------------
# 4. Cell-type -> lineage taxonomy. The canonical 16-type UPMC mapping
#    (schema.UPMC_CANONICAL_CELLTYPES). Verified to match the CLUSTER_LABEL
#    values on disk exactly.
# ---------------------------------------------------------------------------
CELLTYPE_MAP = {
    "APC": "immune",
    "B cell": "immune",
    "CD4 T cell": "immune",
    "CD8 T cell": "immune",
    "Granulocyte": "immune",
    "Macrophage": "immune",
    "Naive immune cell": "immune",
    "Tumor": "tumour",
    "Tumor (CD15+)": "tumour",
    "Tumor (CD20+)": "tumour",
    "Tumor (CD21+)": "tumour",
    "Tumor (Ki67+)": "tumour",
    "Tumor (Podo+)": "tumour",
    "Lymph vessel": "stromal",
    "Stromal / Fibroblast": "stromal",
    "Vessel": "stromal",
}

# ---------------------------------------------------------------------------
# 5. Normalisation — labeled_arcsinh_norm_data.csv is ALREADY arcsinh-normalised
#    (per its name and the project's preprocessing), so do not re-apply it.
# ---------------------------------------------------------------------------
APPLY_ARCSINH = False
ARCSINH_COFACTOR = 5.0

# ---------------------------------------------------------------------------
# 6. QC threshold.
# ---------------------------------------------------------------------------
MIN_CELLS_PER_SAMPLE = 50
