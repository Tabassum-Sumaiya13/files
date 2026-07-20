"""
src/neighbor_features.py — Celltype neighbourhood matrices under different graphs.

Readout (fixed, validated against the base paper's data/raw/k10/*.npy):
    M[i][j] = number of type-j cells among the neighbours of type-i cells
    then ROW-NORMALISED -> "of type-i's neighbours, what fraction are type-j"
    then flattened -> num_clusters**2 features per sample.

Row-normalising divides by neighbour count, which cancels the density term.
That is why the graph choice is a second-order effect for this readout.

Graphs compared (the only variable in Experiment A):
    knn10      - the base paper's choice (directed, k=10). Scale-adaptive.
    radius20   - fixed 20 um. Constant physical meaning across samples.
    radius50   - fixed 50 um. Signalling scale.
    delaunay   - parameter-free Voronoi adjacency.

Coordinates are RAW PIXELS from cell_locations_and_labels.csv — never z-scored.
Per-sample z-scoring would destroy real distance and warp each sample differently,
which makes any radius meaningless.
"""

import numpy as np
import pandas as pd
from pathlib import Path
from typing import Dict, Tuple
from sklearn.neighbors import NearestNeighbors
from scipy.spatial import Delaunay

# Derived from the data, cross-checked against cell biology:
#   median cell area 429 px^2 -> diameter 23.4 px -> 8.8 um at this scale
#   median nearest-neighbour distance 20.6 px      -> 7.8 um
UM_PER_PX = 0.377

N_CLUSTERS = 16


def _drop_self(idx: np.ndarray, n: int, k: int) -> np.ndarray:
    """Remove each row's self-index (robust to duplicate coordinates)."""
    self_col = np.argmax(idx == np.arange(n)[:, None], axis=1)
    keep = np.ones_like(idx, dtype=bool)
    keep[np.arange(n), self_col] = False
    return idx[keep].reshape(n, k)


def _accumulate(pairs_src: np.ndarray, pairs_dst: np.ndarray,
                clusters: np.ndarray, n_clusters: int) -> np.ndarray:
    """M[i][j] += 1 for every (centre of type i, neighbour of type j) pair."""
    M = np.zeros((n_clusters, n_clusters))
    np.add.at(M, (clusters[pairs_src], clusters[pairs_dst]), 1)
    return M


def neighbor_mat_knn(coords, clusters, k=10, n_clusters=N_CLUSTERS) -> np.ndarray:
    """Directed kNN — reproduces the base paper exactly (row sum = count_i * k)."""
    n = len(coords)
    nbrs = NearestNeighbors(n_neighbors=k + 1, metric="euclidean").fit(coords)
    _, idx = nbrs.kneighbors(coords)
    idx = _drop_self(idx, n, k)
    src = np.repeat(np.arange(n), k)
    return _accumulate(src, idx.ravel(), clusters, n_clusters)


def neighbor_mat_radius(coords, clusters, radius_px, n_clusters=N_CLUSTERS) -> np.ndarray:
    """Fixed-radius graph. Constant physical meaning in every sample."""
    n = len(coords)
    nbrs = NearestNeighbors(radius=radius_px, metric="euclidean").fit(coords)
    nbr_idx = nbrs.radius_neighbors(coords, return_distance=False)
    src_list, dst_list = [], []
    for i, nb in enumerate(nbr_idx):
        nb = nb[nb != i]  # drop self
        if len(nb):
            src_list.append(np.full(len(nb), i))
            dst_list.append(nb)
    if not src_list:
        return np.zeros((n_clusters, n_clusters))
    return _accumulate(np.concatenate(src_list), np.concatenate(dst_list),
                       clusters, n_clusters)


def load_prebuilt_delaunay(sid: str, graphs_dir: Path) -> np.ndarray:
    """
    Load a pre-built Delaunay edge list from data/step1_delaunay/graphs/.

    Stored as an (2, E) edge-index array listing BOTH directions, so the
    undirected edge count is E/2. Verified identical to scipy's triangulation of
    the same raw pixel coords (exact edge-set match on every sample checked), so
    these are used directly rather than recomputed.
    """
    e = np.load(graphs_dir / f"{sid}_edges.npy")
    if e.shape[0] == 2 and e.shape[1] != 2:
        e = e.T
    return np.unique(np.sort(e, axis=1), axis=0)   # -> (E_undirected, 2)


def neighbor_mat_from_edges(edges: np.ndarray, clusters: np.ndarray,
                            n_clusters=N_CLUSTERS) -> np.ndarray:
    """Neighbour matrix from an undirected edge list (counts both directions)."""
    src = np.concatenate([edges[:, 0], edges[:, 1]])
    dst = np.concatenate([edges[:, 1], edges[:, 0]])
    return _accumulate(src, dst, clusters, n_clusters)


def neighbor_mat_delaunay(coords, clusters, n_clusters=N_CLUSTERS,
                          max_edge_px: float = None) -> np.ndarray:
    """
    Parameter-free Voronoi adjacency via Delaunay triangulation.

    max_edge_px optionally prunes spurious long edges that Delaunay draws across
    tissue gaps / along the convex hull. Left None = fully parameter-free.
    """
    tri = Delaunay(coords)
    # every simplex (triangle) contributes its 3 edges
    s = tri.simplices
    edges = np.vstack([s[:, [0, 1]], s[:, [1, 2]], s[:, [0, 2]]])
    edges = np.unique(np.sort(edges, axis=1), axis=0)

    if max_edge_px is not None:
        d = np.linalg.norm(coords[edges[:, 0]] - coords[edges[:, 1]], axis=1)
        edges = edges[d <= max_edge_px]

    # undirected -> count both directions
    src = np.concatenate([edges[:, 0], edges[:, 1]])
    dst = np.concatenate([edges[:, 1], edges[:, 0]])
    return _accumulate(src, dst, clusters, n_clusters)


# --- the graph menu for Experiment A ---
GRAPHS = {
    "knn10":    lambda c, cl: neighbor_mat_knn(c, cl, k=10),
    "radius20": lambda c, cl: neighbor_mat_radius(c, cl, radius_px=20.0 / UM_PER_PX),
    "radius50": lambda c, cl: neighbor_mat_radius(c, cl, radius_px=50.0 / UM_PER_PX),
    "delaunay": lambda c, cl: neighbor_mat_delaunay(c, cl),
}

# Settled by Experiment A (20 seeds, both cohorts). See config.graph_type.
DEFAULT_GRAPH = "delaunay"


def row_normalize(M: np.ndarray) -> np.ndarray:
    """Row-normalise then flatten. This is the base paper's exact readout."""
    out = M / M.sum(1)[:, np.newaxis]
    return np.nan_to_num(out, nan=0.0).flatten()


def build_all(locations_csv: Path, out_dir: Path, graphs=None) -> Dict[str, pd.DataFrame]:
    """Build row-normalised neighbour features for every sample, for each graph."""
    graphs = graphs or GRAPHS
    out_dir.mkdir(parents=True, exist_ok=True)

    loc = pd.read_csv(locations_csv, usecols=["ACQUISITION_ID", "X", "Y", "CLUSTER_ID"])
    loc["CLUSTER_ID"] = loc["CLUSTER_ID"].astype(int)
    samples = sorted(loc["ACQUISITION_ID"].unique())
    print(f"  {len(samples)} samples, {len(loc):,} cells, {loc.CLUSTER_ID.nunique()} celltypes")
    print(f"  Coordinates: RAW PIXELS (never z-scored). Scale: {UM_PER_PX} um/px")

    results = {}
    for gname, gfun in graphs.items():
        rows, degrees = {}, []
        for sid in samples:
            s = loc[loc.ACQUISITION_ID == sid]
            coords = s[["X", "Y"]].values.astype(float)
            clusters = s["CLUSTER_ID"].values
            if len(coords) < 5:
                continue
            M = gfun(coords, clusters)
            degrees.append(M.sum() / len(coords))
            rows[sid] = row_normalize(M)

        df = pd.DataFrame.from_dict(rows, orient="index")
        df.index.name = "acquisition_id"
        df.columns = [f"nbr_{gname}_{i}" for i in range(df.shape[1])]
        path = out_dir / f"neighbor_{gname}.parquet"
        df.to_parquet(path)
        results[gname] = df
        print(f"  [{gname:9s}] {df.shape[0]} samples x {df.shape[1]} feats | "
              f"mean degree {np.mean(degrees):5.2f} +/- {np.std(degrees):4.2f} -> {path.name}")

    return results


if __name__ == "__main__":
    import sys
    sys.path.insert(0, str(Path(__file__).parent.parent))
    from config import get_config

    cfg = get_config()
    build_all(cfg.locations_path, cfg.processed_dir / "neighbor_features")
