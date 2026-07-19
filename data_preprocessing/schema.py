"""
schema.py — Canonical data schema this project's pipeline expects, and the
reference marker / cell-type vocabularies a new external cohort is checked
against before it can be used.

`spatial_positional_encoding/` (preprocessing, graph construction, enrichment
features, node-marker features, survival evaluation) all assume these exact
column names and this cell-type taxonomy. Any external dataset must be
adapted to this schema — via a per-dataset `adapter_config.py`, see
datasets/_template/ — before validator.py / processor.py can touch it.
"""

# ---------------------------------------------------------------------------
# Canonical column names (must exist after adapter renaming)
# ---------------------------------------------------------------------------
ACQ_COL = "acquisition_id"          # one row per tissue sample/region
CELL_COL = "cell_id"                # unique within an acquisition
X_COL = "X"
Y_COL = "Y"
CLUSTER_LABEL_COL = "cluster_label"  # native cell-type name, pre lineage-mapping
PATIENT_COL = "patient_id"           # required for grouped cross-validation
SURVIVAL_TIME_COL = "survival_day"
SURVIVAL_STATUS_COL = "survival_status"  # 0 = censored/alive, 1 = event/dead

LOCATIONS_REQUIRED = [ACQ_COL, CELL_COL, X_COL, Y_COL, CLUSTER_LABEL_COL]
EXPRESSION_REQUIRED = [ACQ_COL, CELL_COL]  # marker columns checked separately
METADATA_REQUIRED = [ACQ_COL, PATIENT_COL, SURVIVAL_TIME_COL, SURVIVAL_STATUS_COL]

# ---------------------------------------------------------------------------
# Canonical 16-type / 3-lineage cell taxonomy this project's UPMC cohort uses
# (doc/RESULT_REPORT.md §1.1). A new cohort does not need these exact 16
# names, but every native cell type it has MUST map to one of the 3
# lineages below for the enrichment / celltype-proportion features to be
# comparable across cohorts.
# ---------------------------------------------------------------------------
LINEAGES = ("immune", "tumour", "stromal")

UPMC_CANONICAL_CELLTYPES = {
    "APC": "immune", "B cell": "immune", "CD4 T cell": "immune",
    "CD8 T cell": "immune", "Granulocyte": "immune", "Macrophage": "immune",
    "Naive immune cell": "immune",
    "Tumor": "tumour", "Tumor (CD15+)": "tumour", "Tumor (CD20+)": "tumour",
    "Tumor (CD21+)": "tumour", "Tumor (Ki67+)": "tumour", "Tumor (Podo+)": "tumour",
    "Lymph vessel": "stromal", "Stromal / Fibroblast": "stromal", "Vessel": "stromal",
}

# ---------------------------------------------------------------------------
# Markers used to sanity-check that lineage assignment matches expression
# (see spatial_positional_encoding/src/validate_groups.py). A cohort missing
# most of these can still be processed, but its cell-type -> lineage mapping
# can't be independently verified against the marker data.
# ---------------------------------------------------------------------------
LINEAGE_VALIDATION_MARKERS = [
    "CD45", "PanCK", "Vimentin", "aSMA", "CD31", "Podoplanin",
    "CD3e", "CD20", "CD21", "CD68",
]

# ---------------------------------------------------------------------------
# Markers required by the celltype-conditioned node features
# (spatial_positional_encoding/src/marker_states.py). Used to score how much
# of that feature block a new cohort can reproduce.
# ---------------------------------------------------------------------------
NODE_MARKER_FEATURES = {
    # feature_name: (marker, [lineages it's read in])
    "tumor_ki67":     ("Ki67", ["tumour"]),
    "tumor_mac_pdl1": ("PDL1", ["tumour", "immune"]),
    "cd8_granzymeb":  ("GranzymeB", ["immune"]),
    "cd4_pd1":        ("PD1", ["immune"]),
    "cd4_icos":       ("ICOS", ["immune"]),
    "cd4_foxp3":      ("FoxP3", ["immune"]),     # strongest single marker (RESULT_REPORT.md §8.2)
    "tcell_cd45ro":   ("CD45RO", ["immune"]),    # second strongest
    "apc_mac_hladr":  ("HLA-DR", ["immune"]),
}

# The 2-marker block that produced the best result in this project
# (C = 0.733, RESULT_REPORT.md Table 5) — check for these specifically.
RECOMMENDED_MARKER_SET = ["FoxP3", "CD45RO"]

MIN_CELLS_PER_SAMPLE_DEFAULT = 50
