"""
adapter_config.py — TEMPLATE.

Copy this whole `_template` folder to
data_preprocessing/datasets/<your_dataset_name>/ and fill in every value
below for the new cohort. validator.py and processor.py only import what
this module defines — they never touch the raw files directly.

Then run:
    python run_ingest.py --dataset <your_dataset_name>
"""
from pathlib import Path

import registry

# A short, filesystem-safe name — must match the folder name under datasets/
DATASET_NAME = "CRC"

# ---------------------------------------------------------------------------
# 1. Where the raw files live. Put them in datasets/<name>/raw/ (gitignored —
#    never commit patient-level data) and point these at that folder.
#    Accepted formats: .csv, .parquet, .pkl
# ---------------------------------------------------------------------------
_RAW = Path(__file__).parent / "raw"

LOCATIONS_PATH = _RAW / "cell_locations.csv"    # per-cell X, Y, cell type
EXPRESSION_PATH = _RAW / "cell_expression.csv"  # per-cell marker intensities
METADATA_PATH = _RAW / "sample_metadata.csv"    # per-sample survival/clinical fields

# ---------------------------------------------------------------------------
# 2. Column renames: {native_name_in_raw_file: canonical_name}.
#
#    Canonical names required (schema.py):
#      locations : acquisition_id, cell_id, X, Y, cluster_label
#      expression: acquisition_id, cell_id   (+ marker columns — left as-is)
#      metadata  : acquisition_id, patient_id, survival_day, survival_status
#
#    Only list columns that need RENAMING — if a raw file already uses the
#    canonical name, leave it out of the map.
# ---------------------------------------------------------------------------
LOCATIONS_COLUMN_MAP = {
    # "SampleID": "acquisition_id",
    # "CellID": "cell_id",
    # "centroid_x": "X",
    # "centroid_y": "Y",
    # "cell_type": "cluster_label",
}
EXPRESSION_COLUMN_MAP = {
    # "SampleID": "acquisition_id",
    # "CellID": "cell_id",
}
METADATA_COLUMN_MAP = {
    # "sample": "acquisition_id",
    # "patient": "patient_id",
    # "os_days": "survival_day",
    # "os_event": "survival_status",   # must end up 0 = censored/alive, 1 = event/dead
}


# ---------------------------------------------------------------------------
# 3. Marker columns — how to find them in the expression table.
#    Default rule: every numeric column that isn't an ID/label column.
#    Override the body of this function if the raw file needs a different
#    rule (e.g. an explicit marker_names.csv, or a column-name prefix/suffix).
# ---------------------------------------------------------------------------
_NON_MARKER_COLS = {"acquisition_id", "cell_id", "cluster_id", "cluster_label"}


def marker_columns(expression_df) -> list:
    import pandas as pd
    return [
        c for c in expression_df.columns
        if c not in _NON_MARKER_COLS and pd.api.types.is_numeric_dtype(expression_df[c])
    ]


# ---------------------------------------------------------------------------
# 4. Cell-type taxonomy — DO NOT write the mapping here.
#
#    It lives in ../../celltype_registry.csv, one row per native cell-type
#    label, so that it is cited, versioned, ontology-grounded and falsifiable
#    rather than a bare dict. Add one row per native type of this cohort:
#
#      dataset,native_label,lineage,cl_term_id,cl_label,source,verified,notes,evidence_override
#      MYCOHORT,CD8+ T cell,immune,CL:0000625,"CD8-positive, alpha-beta T cell",<paper>,no,,
#
#    Rules the validator enforces (see registry.py):
#      * EVERY native label in the raw data needs a row, or ingest FAILs — an
#        incomplete map can no longer silently shrink the cohort.
#      * `lineage` must agree with what `cl_term_id` anchors to in
#        schema.CL_LINEAGE_ANCHOR. Add the term there if it is new.
#      * A blank `lineage` means "deliberately excluded"; `notes` must then say
#        why (artifact, unassigned cluster, ambiguous doublet, ...).
#      * validator.py then confronts every row with THIS cohort's own marker
#        data (check `celltype:marker_evidence`) and blocks the ingest if a
#        material contradiction has no recorded justification. To keep a
#        contradicted row, put the reason in `evidence_override` — never reach
#        for --force, which suppresses every check at once.
# ---------------------------------------------------------------------------
CELLTYPE_MAP = registry.celltype_map(DATASET_NAME)

# ---------------------------------------------------------------------------
# 5. Normalisation — does the marker data still need arcsinh, or is it
#    already normalised (e.g. you exported an already-processed table)?
#    The UPMC baseline cohort uses arcsinh(raw_intensity) with no cofactor
#    scaling beyond the divide; CyTOF/IMC panels commonly use cofactor=5.
#    Check your source paper's methods section for the value it used.
# ---------------------------------------------------------------------------
APPLY_ARCSINH = True
ARCSINH_COFACTOR = 5.0

# ---------------------------------------------------------------------------
# 6. QC threshold — samples with fewer cells than this are dropped.
# ---------------------------------------------------------------------------
MIN_CELLS_PER_SAMPLE = 50
