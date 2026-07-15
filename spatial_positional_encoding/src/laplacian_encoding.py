"""
src/laplacian_encoding.py — Compute positional encodings via Laplacian eigendecomposition.

Uses the symmetric normalised Laplacian: L_sym = I - D^{-1/2} A D^{-1/2}
Extracts the k_pe smallest non-trivial eigenvectors as positional encodings.
"""

import numpy as np
from scipy.sparse import csr_matrix, diags
from scipy.sparse.linalg import eigsh
from scipy.sparse.csgraph import connected_components
from typing import Tuple, Dict
import warnings

from src.graph_construction import get_largest_component


def compute_laplacian_encoding(
    adj: csr_matrix,
    k_pe: int = 8,
    sigma_offset: float = 1e-8,
    keep_largest_component: bool = True,
) -> Tuple[np.ndarray, np.ndarray, Dict]:
    """
    Compute positional encodings from graph Laplacian eigenvectors.

    Parameters
    ----------
    adj : (N, N) sparse adjacency matrix
    k_pe : number of eigenvector dimensions to keep
    sigma_offset : shift-invert stabiliser for eigsh
    keep_largest_component : if True, only encode the LCC; others get zeros

    Returns
    -------
    pe_vectors : (N, k_pe) array — zero rows for invalid cells
    valid_mask : (N,) boolean — True only for cells with real PE
    stats : dict with encoding diagnostics
    """
    n_cells = adj.shape[0]
    stats = {
        'n_cells_original': n_cells,
        'n_components': 1,
        'n_cells_in_lcc': n_cells,
        'n_dropped': 0,
        'error': None,
    }

    print(f"    [laplacian] Input: {n_cells} cells, k_pe={k_pe}")

    # --- Optionally restrict to largest connected component ---
    if keep_largest_component:
        adj_work, component_mask = get_largest_component(adj)
        n_dropped = n_cells - adj_work.shape[0]
        stats['n_cells_in_lcc'] = adj_work.shape[0]
        stats['n_dropped'] = n_dropped
        stats['drop_fraction'] = n_dropped / n_cells if n_cells > 0 else 0.0
        if n_dropped > 0:
            print(f"    [laplacian] Dropped {n_dropped} cells "
                  f"({stats['drop_fraction']:.1%}) from smaller components")
    else:
        n_comp, _ = connected_components(adj, directed=False)
        stats['n_components'] = n_comp
        adj_work = adj
        component_mask = np.ones(n_cells, dtype=bool)
        if n_comp > 1:
            print(f"    [laplacian] WARNING: {n_comp} connected components — "
                  f"multiple trivial eigenvalues expected")

    n_work = adj_work.shape[0]

    # --- Check we have enough cells ---
    # eigsh needs k < N strictly; we request k_pe+1 vectors (drop trivial)
    if n_work <= k_pe + 1:
        msg = (f"Too few cells in LCC ({n_work}) for k_pe+1={k_pe + 1} eigenvectors "
               f"(need at least {k_pe + 2})")
        print(f"    [laplacian] ERROR: {msg}")
        stats['error'] = 'too_few_cells'
        return np.zeros((n_cells, k_pe)), np.zeros(n_cells, dtype=bool), stats

    # --- Build symmetric normalised Laplacian ---
    degrees = np.array(adj_work.sum(axis=1)).flatten()
    n_isolated = int((degrees == 0).sum())
    if n_isolated > 0:
        print(f"    [laplacian] WARNING: {n_isolated} isolated nodes — adding epsilon")
        degrees[degrees == 0] = 1e-10

    inv_sqrt_deg = np.power(degrees, -0.5)
    D_inv_sqrt = diags(inv_sqrt_deg)
    identity = diags(np.ones(n_work))
    laplacian = identity - D_inv_sqrt @ adj_work @ D_inv_sqrt

    print(f"    [laplacian] Laplacian shape: {laplacian.shape}, "
          f"nnz={laplacian.nnz:,}")

    # --- Eigendecomposition ---
    try:
        eigenvalues, eigenvectors = eigsh(
            laplacian, k=k_pe + 1, sigma=sigma_offset, which='LM'
        )
        print(f"    [laplacian] Shift-invert eigsh succeeded")
    except Exception as e:
        print(f"    [laplacian] Shift-invert failed ({e}), trying standard mode...")
        try:
            eigenvalues, eigenvectors = eigsh(laplacian, k=k_pe + 1, which='SM')
            print(f"    [laplacian] Standard eigsh succeeded")
        except Exception as e2:
            msg = f"Eigendecomposition failed: {e2}"
            print(f"    [laplacian] ERROR: {msg}")
            stats['error'] = 'eigendecomposition_failed'
            return np.zeros((n_cells, k_pe)), np.zeros(n_cells, dtype=bool), stats

    # Sort by eigenvalue ascending
    idx = np.argsort(eigenvalues)
    eigenvalues = eigenvalues[idx]
    eigenvectors = eigenvectors[:, idx]

    # Drop the trivial (smallest ≈ 0) eigenvector, keep next k_pe
    trivial_eigenvalue = eigenvalues[0]
    eigenvalues = eigenvalues[1:]
    eigenvectors = eigenvectors[:, 1:]

    # Fix sign convention: max absolute value in each column is positive
    eigenvectors = _fix_sign_convention(eigenvectors)

    stats.update({
        'n_cells_encoded': n_work,
        'trivial_eigenvalue': float(trivial_eigenvalue),
        'eigenvalue_min': float(eigenvalues.min()),
        'eigenvalue_max': float(eigenvalues.max()),
        'eigenvalue_range': [float(v) for v in eigenvalues],
    })

    print(f"    [laplacian] Trivial ev_0 = {trivial_eigenvalue:.2e} (should be ~0)")
    print(f"    [laplacian] PE eigenvalues: min={eigenvalues.min():.4f}, max={eigenvalues.max():.4f}")
    print(f"    [laplacian] PE shape: ({n_work}, {k_pe})")

    # Orthonormality check
    orth_error = np.max(np.abs(eigenvectors.T @ eigenvectors - np.eye(k_pe)))
    print(f"    [laplacian] Orthonormality deviation: {orth_error:.2e} "
          f"({'OK' if orth_error < 1e-6 else 'WARNING: > 1e-6'})")

    # --- Map back to full cell set ---
    if n_work != n_cells:
        full_pe = np.zeros((n_cells, k_pe))
        full_pe[component_mask] = eigenvectors
        return full_pe, component_mask, stats

    return eigenvectors, component_mask, stats


def _fix_sign_convention(pe_vectors: np.ndarray) -> np.ndarray:
    """Force the entry with max absolute value in each PE dim to be positive.

    This makes eigenvectors deterministic across runs (they are otherwise
    sign-ambiguous).
    """
    pe = pe_vectors.copy()
    for i in range(pe.shape[1]):
        col = pe[:, i]
        max_abs_idx = np.argmax(np.abs(col))
        if col[max_abs_idx] < 0:
            pe[:, i] = -col
    return pe