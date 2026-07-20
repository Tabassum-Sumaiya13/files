"""
src/marker_states.py — Node features: celltype-conditioned functional-marker states.

WHY THIS EXISTS
---------------
Every cell carries 39 protein markers, but only its celltype LABEL ever reaches
the survival model (via celltype proportions + neighbour matrices). The markers
themselves — the node features — are dropped. This module puts a small, validated
subset back.

WHICH MARKERS, AND WHY THESE
----------------------------
Not all 39. With ~27 patient-level events, a 39x16 block (624 cols) is pure
dilution — the exact trap run_survival's `dC vs NOISE` control exists to catch.
Instead we take the functional/activation markers that carried prognostic signal
in the works this dataset comes from:

  - Schurch et al. 2020 (Cell) — the CODEX colorectal atlas this data IS.
    Their top survival hit was PD1+ CD4 T cells; ICOS, GranzymeB, CD45RO also key.
  - Jackson et al. 2020 (Nature) / Danenberg et al. 2022 (Nat Genetics) — IMC
    breast survival: functional-marker means (Ki67, PDL1, GranzymeB) prognostic.

A marker matters WHERE it is expressed, so each feature is conditioned on the
celltype it is biologically read in (Ki67 in tumour, GranzymeB in CD8 T, ...).
This is the celltype-conditioned functional state, not a bulk average.

READOUT
-------
Per sample, for each feature: MEAN of that marker over the cells of the named
celltype(s). Samples with none of that celltype get 0.0 (arcsinh scale: 0 = no
signal) — "no such cells, so no such functional signal in this tissue". Values
are already arcsinh-normalised in labeled_arcsinh_norm_data.csv.
"""

import numpy as np
import pandas as pd
from pathlib import Path
from typing import Dict, List, Tuple, Optional

# CLUSTER_ID -> celltype (verified against the data, same mapping as
# enrichment_features.py):
#   0 APC   1 B cell   2 CD4 T   3 CD8 T   4 Granulocyte   5 Lymph vessel
#   6 Macrophage   7 Naive immune   8 Stromal/Fibroblast   9 Tumor
#   10 Tumor(CD15+) 11 Tumor(CD20+) 12 Tumor(CD21+) 13 Tumor(Ki67+)
#   14 Tumor(Podo+) 15 Vessel
TUMOR = [9, 10, 11, 12, 13, 14]
MACROPHAGE = [6]
APC = [0]
CD4T = [2]
CD8T = [3]
TCELL = [2, 3]

# (feature name, marker column, celltype cluster ids it is read in)
FEATURES: List[Tuple[str, str, List[int]]] = [
    ("tumor_ki67",      "Ki67",      TUMOR),              # proliferation
    ("tumor_mac_pdl1",  "PDL1",      TUMOR + MACROPHAGE), # checkpoint / evasion
    ("cd8_granzymeb",   "GranzymeB", CD8T),               # cytotoxicity
    ("cd4_pd1",         "PD1",       CD4T),               # Schurch top survival hit
    ("cd4_icos",        "ICOS",      CD4T),               # T-cell activation
    ("cd4_foxp3",       "FoxP3",     CD4T),               # Treg suppression
    ("tcell_cd45ro",    "CD45RO",    TCELL),              # memory / prior exposure
    ("apc_mac_hladr",   "HLA-DR",    APC + MACROPHAGE),   # antigen presentation
]

FEATURE_NAMES = [f[0] for f in FEATURES]
_MARKERS_NEEDED = sorted({f[1] for f in FEATURES})


def build(expression_csv: Path, out_path: Path,
          sample_col: str = "sample_id",
          cluster_col: str = "cluster") -> pd.DataFrame:
    """Build the celltype-conditioned functional-marker table (one row per sample)."""
    usecols = [sample_col, cluster_col] + _MARKERS_NEEDED
    print(f"  Reading {len(usecols)} columns from {expression_csv.name} "
          f"(markers: {_MARKERS_NEEDED})")
    df = pd.read_csv(expression_csv, usecols=usecols)
    df[sample_col] = df[sample_col].astype(str)
    df[cluster_col] = df[cluster_col].astype(int)

    samples = sorted(df[sample_col].unique())
    print(f"  {len(df):,} cells, {len(samples)} samples")

    out = pd.DataFrame(index=pd.Index(samples, name="acquisition_id"),
                       columns=FEATURE_NAMES, dtype=float)

    for name, marker, clusters in FEATURES:
        sub = df[df[cluster_col].isin(clusters)]
        means = sub.groupby(sample_col)[marker].mean()
        out[name] = means.reindex(samples)
        n_missing = out[name].isna().sum()
        print(f"    [{name:16s}] {marker:10s} on {clusters} | "
              f"{n_missing} samples lack this celltype -> 0.0")

    # No cells of that celltype -> no such functional signal in this tissue.
    out = out.fillna(0.0)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out.to_parquet(out_path)
    print(f"\n  Built: {out.shape[0]} samples x {out.shape[1]} features -> {out_path.name}")
    print(f"\n  Distribution (arcsinh-normalised marker means):")
    print(out.describe().T[["mean", "std", "min", "max"]].to_string())
    return out


if __name__ == "__main__":
    import sys
    sys.path.insert(0, str(Path(__file__).parent.parent))
    from config import get_config

    cfg = get_config()
    build(
        expression_csv=cfg.expression_path,
        out_path=cfg.processed_dir / "neighbor_features" / "marker_states.parquet",
    )
