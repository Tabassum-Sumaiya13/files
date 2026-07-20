"""
adapter_config.py — CRC_doublets: the CRC cohort WITH its two ambiguous doublet
clusters included, instead of excluded.

WHY THIS VARIANT EXISTS
-----------------------
Production CRC drops 4 native types (6.9% of cells), two of which are doublets:
    "tumor cells / immune cells"   1,797 cells
    "immune cells / vasculature"   2,153 cells

Dropping them is defensible — forcing an ambiguous doublet into one lineage is a
guess — but the bias is NOT neutral with respect to what this project measures.
Those cells sit at the immune-tumour interface, and three of the five spatial
features (kl_tumor, immune_tumor, stroma_tumor) measure exactly that interface.
Systematically deleting interface cells thins the interface before measuring it.

An argument that the exclusion is harmless has to be demonstrated, not asserted.
This cohort is the demonstration: it reads the SAME raw file, mapping the two
doublets instead of dropping them, so

    run_verify.py --dataset CRC_doublets --perturb-map

produces a `drop:<doublet>` scenario that reproduces production CRC exactly and a
baseline that includes them. If no verdict differs between the two, the exclusion
does not carry the result.

This is a measurement variant, not a second cohort — do not report it as
independent evidence.
"""

from pathlib import Path

import registry

# A short, filesystem-safe name — must match the folder name under datasets/
DATASET_NAME = "CRC_doublets"

# ---------------------------------------------------------------------------
# 1. Where the raw files live. This cohort is one combined table, so all
#    three logical views resolve to the same CSV.
# ---------------------------------------------------------------------------
_RAW = Path(__file__).parent.parent / "CRC" / "raw"   # same raw file as CRC

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
# 4. Cell-type taxonomy: native cell-type label -> "immune" | "tumour" | "stromal".
#
#    NOT hardcoded here any more. The mapping lives in
#    ../../celltype_registry.csv — one row per native label, each carrying a
#    Cell Ontology term, a source citation (Schurch et al. 2020), and a
#    verification flag. Lineage is derived from the CL term via
#    schema.CL_LINEAGE_ANCHOR, so it is not a free choice per row.
#
#    The four types this cohort deliberately EXCLUDES (dirt, undefined, and the
#    two ambiguous doublet clusters) are registry rows with a blank lineage and
#    a recorded reason — the exclusion is now data with a rationale attached,
#    not a code comment. The doublet exclusion carries a directional-bias risk
#    that run_verify.py --perturb-map quantifies.
#
#    To change a mapping, edit the registry — not this file.
# ---------------------------------------------------------------------------
CELLTYPE_MAP = registry.celltype_map(DATASET_NAME)

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
