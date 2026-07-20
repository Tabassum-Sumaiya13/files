"""
src/graph_construction.py — Build k-NN spatial graphs from cell coordinates.

Each sample's cells are treated as nodes; edges connect spatially
nearby cells via a k-NN graph.  The resulting adjacency matrix is
symmetric and has no self-loops.
"""

import numpy as np
from sklearn.neighbors import NearestNeighbors
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
from typing import Tuple, Dict, Optional


def build_spatial_graph(
    coords: np.ndarray,
    k: int = 10,
    symmetrization: str = "union",
) -> Tuple[csr_matrix, Dict]:
    """
    Build a k-NN spatial graph from (N, 2) cell coordinates.

    Parameters
    ----------
    coords : (N, 2) array of X, Y
    k : number of nearest neighbours per cell
    symmetrization : "union" (i->j OR j->i) or "mutual" (AND)

    Returns
    -------
    adj : (N, N) sparse CSR adjacency matrix (binary, symmetric, no self-loops)
    stats : dict with graph statistics
    """
    n_cells = coords.shape[0]
    print(f"    [graph] Building k-NN graph: {n_cells} cells, k={k}, sym={symmetrization}")

    if n_cells < k + 1:
        raise ValueError(f"Too few cells ({n_cells}) for k={k} neighbours")

    # Fit k+1 neighbours (includes self), then drop self
    nbrs = NearestNeighbors(n_neighbors=k + 1, metric='euclidean').fit(coords)
    distances, indices = nbrs.kneighbors(coords)

    # Explicitly drop self-index per row (handles duplicate coordinates)
    self_col = np.argmax(indices == np.arange(n_cells)[:, None], axis=1)
    keep_mask = np.ones_like(indices, dtype=bool)
    keep_mask[np.arange(n_cells), self_col] = False
    indices = indices[keep_mask].reshape(n_cells, k)
    distances = distances[keep_mask].reshape(n_cells, k)

    # Build sparse adjacency
    rows = np.repeat(np.arange(n_cells), k)
    cols = indices.ravel()
    data = np.ones(n_cells * k)
    adj = csr_matrix((data, (rows, cols)), shape=(n_cells, n_cells))

    # Symmetrise
    if symmetrization == "union":
        adj = adj.maximum(adj.T)
    elif symmetrization == "mutual":
        adj = adj.multiply(adj.T)
    else:
        raise ValueError(f"Unknown symmetrization: {symmetrization}")

    adj.setdiag(0)
    adj.eliminate_zeros()

    # Compute stats
    degrees = np.array(adj.sum(axis=1)).flatten()
    n_components, _ = connected_components(adj, directed=False)
    stats = {
        'n_nodes': n_cells,
        'n_edges': adj.nnz // 2,
        'mean_degree': float(degrees.mean()),
        'min_degree': int(degrees.min()),
        'max_degree': int(degrees.max()),
        'n_isolated': int((degrees == 0).sum()),
        'n_components': n_components,
    }

    print(f"    [graph] Edges: {stats['n_edges']:,}, "
          f"degree: {stats['mean_degree']:.1f} (min={stats['min_degree']}, max={stats['max_degree']}), "
          f"components: {stats['n_components']}")

    return adj, stats


def get_largest_component(adj: csr_matrix) -> Tuple[csr_matrix, np.ndarray]:
    """
    Extract the largest connected component.

    Returns
    -------
    adj_lcc : adjacency of the largest component
    mask : boolean array (True = cell is in LCC)
    """
    n_components, labels = connected_components(adj, directed=False)

    if n_components == 1:
        return adj, np.ones(adj.shape[0], dtype=bool)

    component_sizes = np.bincount(labels)
    largest = np.argmax(component_sizes)
    mask = labels == largest
    adj_lcc = adj[mask][:, mask]

    print(f"    [graph] Largest component: {mask.sum()}/{adj.shape[0]} cells "
          f"(dropped {adj.shape[0] - mask.sum()} from {n_components - 1} smaller components)")

    return adj_lcc, mask