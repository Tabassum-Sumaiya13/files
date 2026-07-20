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

# Metadata: patient_id is required (patient-grouped CV can't run without it).
# Survival is OPTIONAL — a cohort with no survival table still ingests and flows
# through feature extraction + the dataset-agnostic verification (run_verify.py).
# Survival is only needed for the separate, per-dataset survival downstream check.
METADATA_REQUIRED = [ACQ_COL, PATIENT_COL]
SURVIVAL_COLS = [SURVIVAL_TIME_COL, SURVIVAL_STATUS_COL]

# ---------------------------------------------------------------------------
# The 3-lineage vocabulary every cohort collapses to. A new cohort does not need
# any particular native cell-type names, but every native type it has MUST map
# to one of these 3 for the enrichment / celltype-proportion features to be
# comparable across cohorts.
#
# The per-dataset native -> lineage mapping itself is NOT here and is NOT
# hardcoded anywhere: it lives in celltype_registry.csv, one cited row per
# (dataset, native_label). See registry.py.
# ---------------------------------------------------------------------------
LINEAGES = ("immune", "tumour", "stromal")

# ---------------------------------------------------------------------------
# Cell Ontology anchor — the external grounding for the registry.
#
# Maps each CL term the registry uses to its lineage, so a registry row's
# lineage is a FUNCTION of its cl_term_id rather than a free per-row choice:
# two rows citing the same term can never be assigned different lineages, and
# registry.validate_registry() FAILs if a row's declared lineage disagrees.
#
# LIMITATION, stated plainly: this is a curated lookup of the terms this project
# actually uses, NOT a traversal of the Cell Ontology. A real traversal would
# derive lineage from `is_a` ancestry and needs the OBO file. This table catches
# inconsistency and gives every row a traceable term id; it does not prove the
# term itself is the correct one for that native label. Each registry row
# therefore carries verified=no until confirmed against the ontology
# (procedure: doc/CELLTYPE_MAPPING.md).
#
# Terms resolvable at https://www.ebi.ac.uk/ols4/ontologies/cl/classes?obo_id=<id>
# ---------------------------------------------------------------------------
CL_LINEAGE_ANCHOR = {
    # --- immune: leukocyte (CL:0000738) and its descendants ----------------
    "CL:0000738": "immune",   # leukocyte
    "CL:0000145": "immune",   # professional antigen presenting cell
    "CL:0000236": "immune",   # B cell
    "CL:0000084": "immune",   # T cell
    "CL:0000624": "immune",   # CD4-positive, alpha-beta T cell
    "CL:0000625": "immune",   # CD8-positive, alpha-beta T cell
    "CL:0000897": "immune",   # CD4-positive, alpha-beta memory T cell
    "CL:0000815": "immune",   # regulatory T cell
    "CL:0000094": "immune",   # granulocyte
    "CL:0000235": "immune",   # macrophage
    "CL:0000576": "immune",   # monocyte
    "CL:0000451": "immune",   # dendritic cell
    "CL:0000623": "immune",   # natural killer cell
    "CL:0000786": "immune",   # plasma cell
    # --- tumour: neoplastic cell -------------------------------------------
    "CL:0001063": "tumour",   # neoplastic cell
    # --- stromal: non-immune, non-neoplastic structural compartment --------
    "CL:0000057": "stromal",  # fibroblast
    "CL:0000115": "stromal",  # endothelial cell
    "CL:0002138": "stromal",  # endothelial cell of lymphatic vessel
    "CL:0000192": "stromal",  # smooth muscle cell
    "CL:0000136": "stromal",  # adipocyte
    "CL:0000125": "stromal",  # glial cell — neural-crest derived, NOT mesenchymal.
                              # Assigned stromal only as the closest of the three
                              # available lineages; low confidence by construction.
}

# ---------------------------------------------------------------------------
# Lineage marker panels — the EVIDENCE the registry is tested against.
#
# `core` markers are the ones the lineage score is computed from: chosen for
# specificity, so a high score means that lineage and not another.
# `supporting` markers are reported as context but never scored — they are
# informative yet cross-expressed (Vimentin is high in activated immune cells;
# CD68 appears on some tumour cells), so scoring them would blur the contrast
# the check exists to draw.
#
# Names here are CANONICAL. Cohorts name markers differently
# ('Cytokeratin - epithelia:Cyc_10_ch_2'); markers.resolve_panel() bridges that.
# ---------------------------------------------------------------------------
LINEAGE_MARKER_PANELS = {
    # Immune and stromal are HETEROGENEOUS lineages: a T cell is CD45+/CD20-/CD68-,
    # an endothelial cell is CD31+/aSMA-. Core panels therefore list the anchors of
    # each major subset, and the score aggregates them with MAX, not mean
    # (see lineage_evidence.py) — "bright for at least one marker of this lineage".
    # Core must carry an anchor for every leukocyte subset the cohorts' label
    # sets contain, decided a priori from the label vocabulary — not tuned to
    # make disagreements disappear. T/B/myeloid/DC (CD3e/CD20/CD68+CD11b/CD11c)
    # plus the two that were initially missing: NK (CD56) and plasma (CD38).
    "immune":  {"core": ["CD45", "CD3e", "CD20", "CD68", "CD11b", "CD11c", "CD56", "CD38"],
                "supporting": ["CD4", "CD8", "CD45RO", "FoxP3", "CD138", "CD57"]},
    # PanCK is the epithelial/carcinoma anchor; CDX2 the intestinal-epithelium one.
    # MUC1 is deliberately NOT core: it is a genuine PLASMA CELL marker as well as
    # an epithelial one, and scoring it made CRC's plasma-cell cluster read as
    # tumour. Kept as supporting so the evidence is still visible.
    "tumour":  {"core": ["PanCK", "CDX2"],
                "supporting": ["MUC1", "EGFR", "Ki67", "p53"]},
    "stromal": {"core": ["aSMA", "CD31", "Podoplanin", "CollagenIV", "CD34"],
                "supporting": ["Vimentin", "CD44"]},
}

# Flat list of every marker involved in lineage validation (any panel, any role).
LINEAGE_VALIDATION_MARKERS = sorted({
    m for panel in LINEAGE_MARKER_PANELS.values()
    for role in panel.values() for m in role
})

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
