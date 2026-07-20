"""
src/node_features.py — Portable celltype-conditioned marker features ("node features").

WHAT A NODE FEATURE IS
----------------------
The mean of one protein marker over ONLY the cells of the type it is biologically
read in, pooled to one number per sample. Ki67 in a tumour cell means
proliferation; Ki67 averaged over every cell is noise. The conditioning is the
feature.

WHY THIS MODULE EXISTS (what it replaces)
-----------------------------------------
src/marker_states.py builds the same 8 features but cannot run on any cohort but
UPMC, for three reasons:

  1. it conditions on hardcoded integer CLUSTER_IDs 0-15 — CRC has string
     cluster_labels and 29 of them;
  2. it reads bare marker names ("Ki67", "HLA-DR") — CRC's columns are named
     'Ki67 - proliferation:Cyc_12_ch_3';
  3. it reads the pre-canonical data/processed/ layout, whose paths no longer
     resolve.

schema.NODE_MARKER_FEATURES is the portable spec, but it collapses conditioning to
LINEAGE, so 5 of the 8 features condition on plain "immune" — cd4_foxp3 would be
FoxP3 averaged over B cells, macrophages and granulocytes too. That is the exact
bulk-average the celltype conditioning exists to avoid, and the feature name would
no longer describe what it measures.

HOW CONDITIONING IS RESOLVED HERE
---------------------------------
By CELL ONTOLOGY TERM, read from data_preprocessing/celltype_registry.csv. Each
feature names a set of CL terms; the registry maps native_label -> cl_term_id per
dataset; so the native labels a feature conditions on are DERIVED per cohort with
no hand mapping:

    cd4_pd1  conditions on {CL:0000624, CL:0000897, CL:0000815}
      UPMC -> ['CD4 T cell']
      CRC  -> ['CD4+ T cells', 'CD4+ T cells CD45RO+', 'CD4+ T cells GATA3+', 'Tregs']

This reuses the grounding the registry already carries and that validator.py
already checks, so a new cohort gets correct conditioning the moment its registry
rows exist. A cohort missing a cell type simply has no cells for that feature,
which is reported as SUPPORT rather than silently imputed (see below).

SUPPORT, NOT ZERO-FILL
----------------------
marker_states.py fills a sample with no cells of the conditioning type with 0.0.
On the arcsinh scale 0 means "no signal", but the truth is "no such cells" — a
different statement, and one that manufactures a fake data point. Here such a
sample is NaN and excluded, and every feature reports how many samples actually
supported it and with how many cells. A feature computed from 3 cells on half the
cohort is not the same measurement as one computed from 3,000 on all of it, and
the matrix has to say so.
"""
from __future__ import annotations

from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Conditioning sets, as Cell Ontology terms.
#
# Subsets are folded into their parent deliberately: UPMC has ONE 'CD4 T cell'
# cluster that contains its Tregs and memory cells, while CRC splits Tregs
# (CL:0000815) and CD45RO+ memory (CL:0000897) into their own clusters. Folding
# them in is what makes cd4_* comparable ACROSS the two cohorts rather than
# measuring a systematically different cell population in each.
# ---------------------------------------------------------------------------
CD4T = {"CL:0000624",   # CD4-positive, alpha-beta T cell
        "CL:0000897",   # CD4-positive, alpha-beta memory T cell
        "CL:0000815"}   # regulatory T cell
CD8T = {"CL:0000625"}   # CD8-positive, alpha-beta T cell
TCELL = CD4T | CD8T | {"CL:0000084"}      # + generic T cell
MACROPHAGE = {"CL:0000235"}
APC = {"CL:0000145",    # professional antigen presenting cell
       "CL:0000451",    # dendritic cell
       "CL:0000235"}    # macrophage
TUMOUR = {"CL:0001063"}  # neoplastic cell

# (feature name, canonical marker, conditioning CL terms)
# Names and marker choices are unchanged from doc/NODE_FEATURES_REPORT.md so the
# matrix is directly comparable to the survival work that motivated them.
FEATURES: List[Tuple[str, str, set]] = [
    ("tumor_ki67",     "Ki67",      TUMOUR),               # proliferation
    ("tumor_mac_pdl1", "PDL1",      TUMOUR | MACROPHAGE),  # checkpoint / evasion
    ("cd8_granzymeb",  "GranzymeB", CD8T),                 # cytotoxicity
    ("cd4_pd1",        "PD1",       CD4T),                 # Schurch top survival hit
    ("cd4_icos",       "ICOS",      CD4T),                 # T-cell activation
    ("cd4_foxp3",      "FoxP3",     CD4T),                 # Treg suppression
    ("tcell_cd45ro",   "CD45RO",    TCELL),                # memory / prior exposure
    ("apc_mac_hladr",  "HLA-DR",    APC),                  # antigen presentation
]

NODE_FEATURE_NAMES = [f[0] for f in FEATURES]
CANONICAL_MARKERS = sorted({f[1] for f in FEATURES})


def resolve_conditioning(registry_rows: pd.DataFrame) -> Dict[str, List[str]]:
    """{feature_name: [native labels in THIS cohort]} via the registry's CL terms.

    `registry_rows` is registry.dataset_rows(<dataset>). Rows with no lineage
    (deliberately excluded types) carry no CL term and so are never selected.
    """
    by_term: Dict[str, List[str]] = {}
    for r in registry_rows.itertuples():
        if r.lineage and r.cl_term_id:
            by_term.setdefault(r.cl_term_id, []).append(r.native_label)
    return {name: sorted({lab for term in terms for lab in by_term.get(term, [])})
            for name, _, terms in FEATURES}


def resolve_markers(marker_columns: Sequence[str], markers_module) -> Tuple[Dict[str, str], List[str]]:
    """{canonical marker: actual column name} for this cohort, plus what is missing.

    Delegates to data_preprocessing/markers.py so a cohort naming its markers
    'Ki67 - proliferation:Cyc_12_ch_3' resolves without a per-cohort alias table.
    """
    return markers_module.resolve_panel(CANONICAL_MARKERS, list(marker_columns))


def sample_node_features(
    labels: np.ndarray,
    marker_values: Dict[str, np.ndarray],
    conditioning: Dict[str, List[str]],
    marker_map: Dict[str, str],
) -> Tuple[np.ndarray, np.ndarray]:
    """One sample -> (values, support) arrays, both len(FEATURES).

    values[i]  is NaN when the sample has no cells of feature i's conditioning
               type, or when its marker is absent from this cohort's panel.
    support[i] is the number of cells the mean was taken over.
    """
    vals = np.full(len(FEATURES), np.nan)
    sup = np.zeros(len(FEATURES), dtype=int)
    for i, (name, marker, _) in enumerate(FEATURES):
        col = marker_map.get(marker)
        if col is None:
            continue
        mask = np.isin(labels, conditioning.get(name, []))
        n = int(mask.sum())
        sup[i] = n
        if n:
            vals[i] = float(np.mean(marker_values[col][mask]))
    return vals, sup


def conditioning_summary(conditioning: Dict[str, List[str]]) -> pd.DataFrame:
    """One row per feature: which native labels it resolved to in this cohort."""
    return pd.DataFrame([
        {"feature": name,
         "marker": marker,
         "n_labels": len(conditioning.get(name, [])),
         "conditioned_on": ", ".join(conditioning.get(name, [])) or "(none)"}
        for name, marker, _ in FEATURES
    ])
