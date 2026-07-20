"""
adapter_config.py - Keren 2018 TNBC MIBI-TOF cohort.
Source: Keren et al. 2018, Cell - A Structured Tumor-Immune Microenvironment
DOI: 10.1016/j.cell.2018.08.039
"""
from pathlib import Path

DATASET_NAME = "keren_tnbc"

_RAW = Path(__file__).parent / "raw"

LOCATIONS_PATH = _RAW / "cell_locations.csv"
EXPRESSION_PATH = _RAW / "cell_expression.csv"
METADATA_PATH = _RAW / "sample_metadata.csv"

LOCATIONS_COLUMN_MAP = {
    "SampleID": "acquisition_id",
    "cellLabelInImage": "cell_id",
}
EXPRESSION_COLUMN_MAP = {
    "SampleID": "acquisition_id",
    "cellLabelInImage": "cell_id",
}
METADATA_COLUMN_MAP = {}

_NON_MARKER_COLS = {"acquisition_id", "cell_id"}

def marker_columns(expression_df):
    import pandas as pd
    return [c for c in expression_df.columns if c not in _NON_MARKER_COLS and pd.api.types.is_numeric_dtype(expression_df[c])]

CELLTYPE_MAP = {
    "Immune": "immune",
    "Tregs": "immune",
    "CD4_T": "immune",
    "CD8_T": "immune",
    "CD3_T": "immune",
    "NK": "immune",
    "B": "immune",
    "Neutrophils": "immune",
    "Macrophages": "immune",
    "DC": "immune",
    "DC_Mono": "immune",
    "Mono_Neu": "immune",
    "Other_immune": "immune",
    "Tumor": "tumour",
    "Keratin_positive_tumor": "tumour",
    "Endothelial": "stromal",
    "Mesenchymal_like": "stromal",
}

APPLY_ARCSINH = False
ARCSINH_COFACTOR = 5.0
MIN_CELLS_PER_SAMPLE = 50
