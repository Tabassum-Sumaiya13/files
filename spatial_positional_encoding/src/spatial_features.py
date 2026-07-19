"""
src/spatial_features.py — Portable spatial features from a canonical sample.

Generalises src/enrichment_features.py (which was hardwired to UPMC's 16
integer CLUSTER_IDs and prebuilt Delaunay *.npy edges) so it runs on ANY
ingested cohort, keyed on whichever label column you pass:

    label_col = "lineage"        -> 3-category, PORTABLE, comparable across cohorts
    label_col = "cluster_label"  -> native per-dataset resolution (16 UPMC, 29 CRC)

Two readouts per sample:
  * composition proportions   — the abundance-only BASELINE every spatial
                                feature must beat (p_j = fraction of category j)
  * 5 enrichment scalars      — abundance-corrected spatial organisation
                                (0 == exactly the random-mix null), identical
                                definitions to enrichment_features.py:
        kl_mean, kl_tumor, self_enrich, immune_tumor, stroma_tumor

The graph is Delaunay, built here from X/Y (parameter-free, mean degree ~6),
so no prebuilt edge files are needed.
"""
from __future__ import annotations

from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

ENRICH_FEATURE_NAMES = ["kl_mean", "kl_tumor", "self_enrich", "immune_tumor", "stroma_tumor"]
_FLOOR = 0.5  # half-count floor: keeps log2 finite when a block has zero contacts


# ---------------------------------------------------------------------------
# Graph
# ---------------------------------------------------------------------------
def delaunay_edges(coords: np.ndarray) -> np.ndarray:
    """Undirected, de-duplicated Delaunay edges -> (E, 2) int array.

    Returns an empty (0, 2) array when a triangulation is impossible
    (< 3 points, or all points collinear/coincident).
    """
    from scipy.spatial import Delaunay
    from scipy.spatial.qhull import QhullError  # type: ignore

    n = coords.shape[0]
    if n < 3:
        return np.empty((0, 2), dtype=int)
    try:
        tri = Delaunay(coords)
    except Exception:  # QhullError for degenerate (collinear) inputs
        return np.empty((0, 2), dtype=int)

    s = tri.simplices  # (T, 3)
    e = np.vstack([s[:, [0, 1]], s[:, [1, 2]], s[:, [0, 2]]])
    e = np.unique(np.sort(e, axis=1), axis=0)
    return e.astype(int)


def _count_matrix(edges: np.ndarray, idx: np.ndarray, K: int) -> np.ndarray:
    """Raw K×K neighbour counts from an undirected edge list (both directions).

    Any contact touching an out-of-vocabulary cell (idx < 0 — a cell whose label
    is not one of the K categories, e.g. a new/unmapped cell type) is dropped, so
    such cells never contribute rather than being mis-indexed into M[-1]. This
    keeps the count matrix correct even under label permutation, where an OOV
    label can otherwise be shuffled onto an edge endpoint.
    """
    M = np.zeros((K, K), dtype=float)
    if edges.shape[0] == 0:
        return M
    src = np.concatenate([edges[:, 0], edges[:, 1]])
    dst = np.concatenate([edges[:, 1], edges[:, 0]])
    si, di = idx[src], idx[dst]
    ok = (si >= 0) & (di >= 0)
    np.add.at(M, (si[ok], di[ok]), 1)
    return M


# ---------------------------------------------------------------------------
# Enrichment scalars (abundance-corrected; 0 == chance)
# ---------------------------------------------------------------------------
def _kl_bits(q: np.ndarray, p: np.ndarray) -> float:
    ok = (q > 0) & (p > 0)
    if not ok.any():
        return 0.0
    return float(np.sum(q[ok] * np.log2(q[ok] / p[ok])))


def _log2_enrich(num: float, den: float, p_expected: float) -> float:
    if den <= 0 or p_expected <= 0:
        return 0.0
    frac = max(num, _FLOOR) / den
    return float(np.log2(frac / p_expected))


def _block_enrich(M: np.ndarray, p: np.ndarray, rows: List[int], cols: List[int]) -> float:
    if not rows or not cols:
        return 0.0
    num = float(M[np.ix_(rows, cols)].sum())
    den = float(M[rows, :].sum())
    return _log2_enrich(num, den, float(p[cols].sum()))


def enrichment_scalars(
    M: np.ndarray,
    counts: np.ndarray,
    immune: List[int],
    tumour: List[int],
    stromal: List[int],
) -> np.ndarray:
    """Collapse a K×K neighbour-count matrix into the 5 abundance-corrected scalars.

    `immune`/`tumour`/`stromal` are lists of category indices belonging to each
    lineage (for label_col='lineage' each is a single index; for native labels
    they group the native categories by their lineage).
    """
    total = counts.sum()
    if total <= 0:
        return np.zeros(len(ENRICH_FEATURE_NAMES))
    p = counts.astype(float) / total

    row_sums = M.sum(1, keepdims=True)
    P = np.divide(M, row_sums, out=np.zeros_like(M, dtype=float), where=row_sums > 0)
    valid = row_sums.ravel() > 0

    kl_mean = (
        float(np.mean([_kl_bits(P[i], p) for i in np.flatnonzero(valid)]))
        if valid.any() else 0.0
    )

    tumour_row = M[tumour, :].sum(0) if tumour else np.zeros(M.shape[1])
    ts = tumour_row.sum()
    kl_tumor = _kl_bits(tumour_row / ts, p) if ts > 0 else 0.0

    diag = [
        _log2_enrich(M[i, i], row_sums[i, 0], p[i])
        for i in np.flatnonzero(valid) if p[i] > 0
    ]
    self_enrich = float(np.mean(diag)) if diag else 0.0

    immune_tumor = _block_enrich(M, p, immune, tumour)
    stroma_tumor = _block_enrich(M, p, stromal, tumour)
    return np.array([kl_mean, kl_tumor, self_enrich, immune_tumor, stroma_tumor])


# ---------------------------------------------------------------------------
# Sample-level drivers
# ---------------------------------------------------------------------------
def _lineage_grouping(labels: np.ndarray, lineage: np.ndarray, categories: List[str]):
    """Map each category index -> lineage bucket via the majority lineage of its cells."""
    immune, tumour, stromal = [], [], []
    for j, cat in enumerate(categories):
        lin_vals = lineage[labels == cat]
        if len(lin_vals) == 0:
            continue
        maj = pd.Series(lin_vals).mode()
        lin = str(maj.iloc[0]) if len(maj) else ""
        if lin == "immune":
            immune.append(j)
        elif lin == "tumour":
            tumour.append(j)
        elif lin == "stromal":
            stromal.append(j)
    return immune, tumour, stromal


def composition_proportions(labels: np.ndarray, categories: List[str]) -> np.ndarray:
    counts = pd.Series(labels).value_counts().reindex(categories, fill_value=0).values.astype(float)
    total = counts.sum()
    return counts / total if total > 0 else counts


def sample_features(
    df: pd.DataFrame,
    label_col: str,
    categories: List[str],
    lineage_col: str = "lineage",
    edges: Optional[np.ndarray] = None,
) -> Dict[str, np.ndarray]:
    """Composition proportions + 5 enrichment scalars for one sample.

    Parameters
    ----------
    df : sample dataframe (needs X, Y, `label_col`, `lineage_col`)
    label_col : "lineage" (portable) or "cluster_label" (native)
    categories : fixed category order shared across the cohort (so columns align)
    edges : precomputed Delaunay edges; built from X/Y if None (enables the
            permutation-null test to reuse one graph across label shuffles)
    """
    cat_index = {c: i for i, c in enumerate(categories)}
    labels = df[label_col].astype(str).values
    idx = np.fromiter((cat_index.get(l, -1) for l in labels), dtype=int, count=len(labels))
    K = len(categories)

    counts = np.bincount(idx[idx >= 0], minlength=K).astype(float)
    p = counts / counts.sum() if counts.sum() > 0 else counts

    if edges is None:
        edges = delaunay_edges(df[["X", "Y"]].values.astype(float))
    # edges index cells; drop any edge touching an out-of-vocabulary label
    if edges.shape[0] and (idx < 0).any():
        keep = (idx[edges[:, 0]] >= 0) & (idx[edges[:, 1]] >= 0)
        edges = edges[keep]

    M = _count_matrix(edges, idx, K)
    immune, tumour, stromal = _lineage_grouping(labels, df[lineage_col].astype(str).values, categories)
    enrich = enrichment_scalars(M, counts, immune, tumour, stromal)
    return {"proportions": p, "enrichment": enrich}
