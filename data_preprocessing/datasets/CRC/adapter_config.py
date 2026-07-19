"""
adapter_config.py — CRC (Schürch et al. 2020 colorectal-cancer CODEX cohort).

This cohort ships as a SINGLE table
(raw/CRC_clusters_neighborhoods_markers.csv, ~258k cells, 140 tissue
regions, 35 patients) that already contains cell locations, marker
intensities AND the per-cell grouping columns. So all three source paths
below point at that one file — the processor reads it three times and each
COLUMN_MAP pulls out the columns that view needs.

Run:
    python run_ingest.py --dataset CRC

NOTE — survival labels are NOT in this raw file (see METADATA_COLUMN_MAP
below). The only outcome variable present is `groups` (1 = CLR /
2 = DII). survival_day / survival_status must be supplied before
run_ingest can complete.
"""
from pathlib import Path

# A short, filesystem-safe name — must match the folder name under datasets/
DATASET_NAME = "CRC"

# ---------------------------------------------------------------------------
# 1. Where the raw files live. This cohort is one combined table, so all
#    three logical views resolve to the same CSV.
# ---------------------------------------------------------------------------
_RAW = Path(__file__).parent / "raw"

_COMBINED = _RAW / "CRC_clusters_neighborhoods_markers.csv"
LOCATIONS_PATH = _COMBINED     # per-cell X, Y, cell type
EXPRESSION_PATH = _COMBINED    # per-cell marker intensities
METADATA_PATH = _COMBINED      # per-sample grouping / (survival) fields

# ---------------------------------------------------------------------------
# 2. Column renames: {native_name_in_raw_file: canonical_name}.
#
#    Canonical names required (schema.py):
#      locations : acquisition_id, cell_id, X, Y, cluster_label
#      expression: acquisition_id, cell_id   (+ marker columns — left as-is)
#      metadata  : acquisition_id, patient_id, survival_day, survival_status
#
#    Native columns used here:
#      "File Name"   -> acquisition_id  (reg001_A … 140 unique tissue regions,
#                                        each belongs to exactly one patient)
#      "CellID"      -> cell_id         (globally unique running index)
#      "X:X" / "Y:Y" -> X / Y           (global tile-stitched coordinates)
#      "ClusterName" -> cluster_label   (native 29-way cell-type labels)
#      "patients"    -> patient_id      (1..35, globally unique — group 1 uses
#                                        1..35, group 2 uses 2..31, no overlap)
# ---------------------------------------------------------------------------
LOCATIONS_COLUMN_MAP = {
    "File Name": "acquisition_id",
    "CellID": "cell_id",
    "X:X": "X",
    "Y:Y": "Y",
    "ClusterName": "cluster_label",
}
EXPRESSION_COLUMN_MAP = {
    "File Name": "acquisition_id",
    "CellID": "cell_id",
}
METADATA_COLUMN_MAP = {
    "File Name": "acquisition_id",
    "patients": "patient_id",
    # --- SURVIVAL IS MISSING FROM THE RAW FILE ------------------------------
    # There is no overall-survival day or event column in
    # CRC_clusters_neighborhoods_markers.csv. The only outcome proxy present
    # is `groups` (1 = CLR, good prognosis; 2 = DII, poor prognosis).
    #
    # A COLUMN_MAP can only RENAME a column — it cannot remap {1,2} -> {0,1}
    # nor invent survival_day. To make this cohort runnable you must add a
    # real per-patient survival table (Schürch et al. 2020, Supplementary
    # Table S1 has OS days + vital status) and merge it into raw/ as extra
    # columns, e.g. "OS_day" and "OS_event", then uncomment:
    #     "OS_day":   "survival_day",
    #     "OS_event": "survival_status",   # 0 = alive/censored, 1 = dead
    #
    # If instead you only want the CLR-vs-DII binary classification label,
    # derive survival_status from `groups` (2 -> 1, 1 -> 0) in a small
    # preprocessing step and add a placeholder survival_day column.
}


# ---------------------------------------------------------------------------
# 3. Marker columns — how to find them in the expression table.
#    The combined CSV carries dozens of numeric NON-marker columns (CellID,
#    ClusterID, EventID, groups, patients, X/Y, size, neighborhood ids, the
#    binary phenotype flags like "CD4+ICOS+", …), so the default "every
#    numeric column" rule would wrongly pick those up. The real fluorescence
#    markers are exactly the columns whose name carries a ":Cyc_<n>_ch_<n>"
#    channel tag. HOECHST1 and DRAQ5 are nuclear segmentation stains, not
#    phenotypic markers, so they are excluded -> 56 markers (the published
#    CRC panel size).
# ---------------------------------------------------------------------------
_NUCLEAR_STAINS = ("HOECHST1:", "DRAQ5:")


def marker_columns(expression_df) -> list:
    return [
        c for c in expression_df.columns
        if ":Cyc_" in c and not c.startswith(_NUCLEAR_STAINS)
    ]


# ---------------------------------------------------------------------------
# 4. Cell-type taxonomy: map EVERY native cell-type label to
#    "immune" | "tumour" | "stromal". Native types intentionally left out
#    (dirt, undefined, and the mixed/doublet clusters) are dropped by the
#    processor — forcing an artifact or an ambiguous doublet into a lineage
#    would be worse than losing those cells.
# ---------------------------------------------------------------------------
CELLTYPE_MAP = {
    # --- tumour -------------------------------------------------------------
    "tumor cells": "tumour",

    # --- immune -------------------------------------------------------------
    "CD68+CD163+ macrophages": "immune",
    "granulocytes": "immune",
    "CD8+ T cells": "immune",
    "CD4+ T cells CD45RO+": "immune",
    "B cells": "immune",
    "plasma cells": "immune",
    "immune cells": "immune",
    "Tregs": "immune",
    "CD4+ T cells": "immune",
    "CD68+ macrophages": "immune",
    "CD11b+CD68+ macrophages": "immune",
    "CD11b+ monocytes": "immune",
    "CD11c+ DCs": "immune",
    "NK cells": "immune",
    "CD3+ T cells": "immune",
    "CD68+ macrophages GzmB+": "immune",
    "CD4+ T cells GATA3+": "immune",
    "CD163+ macrophages": "immune",

    # --- stromal ------------------------------------------------------------
    "smooth muscle": "stromal",
    "stroma": "stromal",
    "vasculature": "stromal",
    "adipocytes": "stromal",
    "nerves": "stromal",
    "lymphatics": "stromal",

    # --- intentionally UNMAPPED (dropped during processing) -----------------
    #   "dirt"                        -> imaging artifact / debris
    #   "undefined"                   -> unassigned cluster
    #   "immune cells / vasculature"  -> ambiguous doublet cluster
    #   "tumor cells / immune cells"  -> ambiguous doublet cluster
}

# ---------------------------------------------------------------------------
# 5. Normalisation — the marker columns in this file are raw-scale
#    fluorescence intensities (0 … several thousand), so arcsinh compression
#    is applied, matching the project's baseline. Confirm the cofactor
#    against the source if you need an exact reproduction of the paper.
# ---------------------------------------------------------------------------
APPLY_ARCSINH = True
ARCSINH_COFACTOR = 5.0

# ---------------------------------------------------------------------------
# 6. QC threshold — samples with fewer cells than this are dropped.
# ---------------------------------------------------------------------------
MIN_CELLS_PER_SAMPLE = 50
