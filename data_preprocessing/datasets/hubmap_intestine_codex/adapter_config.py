"""
adapter_config.py — Hickey et al. / HuBMAP CODEX human intestine cohort
(Dryad doi.org/10.5061/dryad.pk0p2ngrf, "23_09_CODEX_HuBMAP_alldata_Dryad_merged.csv").

64 tissue regions (66 unique_region values — 2 renamed replicates per the
Dryad README) x 8 donors, 2.6M cells, 53 antibody markers.

*** IMPORTANT CAVEAT — NO REAL SURVIVAL DATA ***
This is healthy human intestine (all 8 donors have "History of cancer: no").
donor_metadata.csv has no outcome/survival field at all — there is nothing
to link a survival_day/survival_status to. survival_day/survival_status in
raw/sample_metadata.csv are PLACEHOLDER CONSTANTS (0, 0 = "censored at day
0" for every region), added only so this cohort satisfies schema.py's
METADATA_REQUIRED and can flow through validator.py/processor.py.

This adapter exists to check that the enrichment / cell-type-proportion /
node-marker feature extraction code (which only needs XY + cell type +
markers, not survival) runs correctly on a second cohort's shape. Do NOT
run run_survival.py's C-index evaluation against this cohort's output and
report it as a real result — there is no real endpoint here. See
../../README.md for real survival-endpoint candidates (Schurch, Jackson,
Danenberg, Keren).

Raw files here were built from the single merged Dryad CSV by a one-off
script (see conversation history / scratchpad split_hubmap_codex.py) —
column-subset + parquet export so validator.py/processor.py aren't reading
a 2.9GB / 75-column file twice.

See datasets/_template/adapter_config.py for the field-by-field contract.
"""
from pathlib import Path

DATASET_NAME = "hubmap_intestine_codex"

_RAW = Path(__file__).parent / "raw"
LOCATIONS_PATH = _RAW / "cell_locations.parquet"
EXPRESSION_PATH = _RAW / "cell_expression.parquet"
METADATA_PATH = _RAW / "sample_metadata.csv"

LOCATIONS_COLUMN_MAP = {
    "unique_region": "acquisition_id",
    "x": "X",
    "y": "Y",
    "Cell Type": "cluster_label",
}
EXPRESSION_COLUMN_MAP = {
    "unique_region": "acquisition_id",
    "HLADR": "HLA-DR",
    "Cytokeratin": "PanCK",
    "CD3": "CD3e",
}
METADATA_COLUMN_MAP = {}

# 3. Marker columns — explicit allow-list. 6 non-clustering markers
#    (OLFM4, FAP, CD25, CollIV, CK7, MUC6) have per-batch missingness that
#    would zero out the whole cohort via processor.py's all-columns dropna
#    (every row is null in at least one of them), so they're excluded here.
#    The 47 markers below are exactly the ones the Dryad paper used for
#    clustering and have zero nulls.
_NON_MARKER_COLS = {"acquisition_id", "cell_id", "cluster_id", "cluster_label"}
_CORE_PANEL_MARKERS = {
    "MUC2", "SOX9", "MUC1", "CD31", "Synapto", "CD49f", "CD15", "CHGA", "CDX2",
    "ITLN1", "CD4", "CD127", "Vimentin", "HLA-DR", "CD8", "CD11c", "CD44", "CD16",
    "BCL2", "CD3e", "CD123", "CD38", "CD90", "aSMA", "CD21", "NKG2D", "CD66",
    "CD57", "CD206", "CD68", "CD34", "aDef5", "CD7", "CD36", "CD138", "CD45RO",
    "PanCK", "CD117", "CD19", "Podoplanin", "CD45", "CD56", "CD69", "Ki67",
    "CD49a", "CD163", "CD161",
}


def marker_columns(expression_df) -> list:
    import pandas as pd
    return [
        c for c in expression_df.columns
        if c in _CORE_PANEL_MARKERS and c not in _NON_MARKER_COLS
        and pd.api.types.is_numeric_dtype(expression_df[c])
    ]


# 4. CELLTYPE_MAP — every one of the 25 native "Cell Type" values is mapped
#    (nothing dropped at the celltype_to_lineage step). schema.py's taxonomy
#    only has 3 buckets (immune / tumour / stromal); this cohort is healthy
#    tissue with no malignant compartment, so the 9 epithelial types
#    (Enterocyte and its CD57+/CD66+/MUC1+ subsets, Cycling TA, TA, Goblet,
#    Paneth, Neuroendocrine) are mapped to "tumour" as the closest
#    structural analogue: in the UPMC taxonomy "tumour" is simply the
#    non-immune, non-stromal PARENCHYMAL compartment (schema.py
#    UPMC_CANONICAL_CELLTYPES), which is exactly the role intestinal
#    epithelium plays here. This is a labeling convenience for
#    cross-cohort feature comparability (enrichment / proportion / node
#    features are computed per-lineage) — it does NOT imply these cells are
#    malignant, and lineage="tumour" should be read as "epithelial
#    parenchyma" for this cohort specifically.
#    immune: B, CD4+ T cell, CD7+ Immune, CD8+ T, DC, M1 Macrophage,
#            M2 Macrophage, NK, Neutrophil, Plasma
#    stromal: Endothelial, ICC, Lymphatic, Nerve, Smooth muscle, Stroma
#    tumour (= epithelial parenchyma here): CD57+ Enterocyte,
#            CD66+ Enterocyte, Cycling TA, Enterocyte, Goblet,
#            MUC1+ Enterocyte, Neuroendocrine, Paneth, TA
CELLTYPE_MAP = {
    "B": "immune", "CD4+ T cell": "immune", "CD7+ Immune": "immune",
    "CD8+ T": "immune", "DC": "immune", "M1 Macrophage": "immune",
    "M2 Macrophage": "immune", "NK": "immune", "Neutrophil": "immune",
    "Plasma": "immune",
    "Endothelial": "stromal", "ICC": "stromal", "Lymphatic": "stromal",
    "Nerve": "stromal", "Smooth muscle": "stromal", "Stroma": "stromal",
    "CD57+ Enterocyte": "tumour", "CD66+ Enterocyte": "tumour",
    "Cycling TA": "tumour", "Enterocyte": "tumour", "Goblet": "tumour",
    "MUC1+ Enterocyte": "tumour", "Neuroendocrine": "tumour",
    "Paneth": "tumour", "TA": "tumour",
}

# 5. Normalisation — already z-normalized per Dryad README, so:
APPLY_ARCSINH = False
ARCSINH_COFACTOR = None

# 6. QC threshold
MIN_CELLS_PER_SAMPLE = 50
