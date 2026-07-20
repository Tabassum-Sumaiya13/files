"""
src/feature_assembly.py — Combine positional encodings + protein markers into
the final per-cell feature vector.

Output schema per cell:
    [pe_0, pe_1, ..., pe_{k-1}, marker_0, marker_1, ..., marker_{m-1}]
"""

import numpy as np
import pandas as pd
from typing import Dict, Tuple, List


def assemble_features(
    pe_vectors: np.ndarray,
    markers: np.ndarray,
    marker_names: List[str],
) -> Tuple[np.ndarray, Dict]:
    """
    Concatenate PE vectors and marker values into a single feature matrix.

    Parameters
    ----------
    pe_vectors : (N, k_pe) positional encodings
    markers : (N, n_markers) protein expression values
    marker_names : list of marker column names

    Returns
    -------
    features : (N, k_pe + n_markers) combined feature matrix
    metadata : dict with dimension info and column names
    """
    n_cells = pe_vectors.shape[0]
    k_pe = pe_vectors.shape[1]
    n_markers = len(marker_names)

    if pe_vectors.shape[0] != markers.shape[0]:
        raise ValueError(f"Row mismatch: PE has {pe_vectors.shape[0]}, "
                         f"markers has {markers.shape[0]}")
    if markers.shape[1] != n_markers:
        raise ValueError(f"markers has {markers.shape[1]} columns but "
                         f"{n_markers} marker_names given")

    features = np.hstack([pe_vectors, markers])

    feature_names = [f'pe_{i}' for i in range(k_pe)] + list(marker_names)
    metadata = {
        'n_cells': n_cells,
        'n_pe_dims': k_pe,
        'n_marker_dims': n_markers,
        'total_dims': k_pe + n_markers,
        'feature_names': feature_names,
    }

    print(f"    [features] Assembled: ({n_cells}, {k_pe} PE + {n_markers} markers = {features.shape[1]} total)")

    return features, metadata


def create_output_dataframe(
    sample_id: str,
    cell_ids: np.ndarray,
    cluster_labels: np.ndarray,
    pe_vectors: np.ndarray,
    markers: np.ndarray,
    marker_names: List[str],
    valid_mask: np.ndarray,
) -> Tuple[pd.DataFrame, Dict]:
    """
    Create the output DataFrame for one sample with all features.

    Columns: acquisition_id, cell_id, cluster_label, is_valid_pe,
             pe_0..pe_{k-1}, <marker_cols>
    """
    features, metadata = assemble_features(pe_vectors, markers, marker_names)

    df = pd.DataFrame({
        'acquisition_id': [sample_id] * len(cell_ids),
        'cell_id': cell_ids,
        'cluster_label': cluster_labels,
        'is_valid_pe': valid_mask.astype(bool),
    })

    for i, name in enumerate(metadata['feature_names']):
        df[name] = features[:, i]

    return df, metadata