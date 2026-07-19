"""
src/enrichment_features.py — Abundance-corrected spatial features (5 per sample).

WHY THIS EXISTS
---------------
The base paper's readout is P[i][j] = "of type i's neighbours, what fraction are
type j", flattened to 256 columns. Under a random spatial arrangement every row
of P collapses to the global composition:

    P[i][j] -> p_j        (p = global proportion of each celltype)

The celltype-proportion baseline already hands the model p. So most of those 256
columns re-encode information the baseline has, which is why
`Celltype + delaunay` gains only +0.012 over `Celltype` alone — it is largely the
same information twice, spread over 256 noisy columns.

Dividing by p removes the shared part:

    E[i][j] = P[i][j] / p_j       E = 1 -> exactly chance
                                  E > 1 -> genuinely enriched
                                  E < 1 -> genuinely depleted

Every feature below is 0 under the random-mix null, and none is computable from
abundance alone. Five features against ~27 patient-level events is EPV ~5 —
still under the >=10 rule of thumb, but the same order of magnitude rather than
90x off, which is where the 272-wide blocks sit.

GRAPH: Delaunay (Experiment A winner, config default, parameter-free).
"""

import numpy as np
import pandas as pd
from pathlib import Path
from typing import Dict, Optional

N_CLUSTERS = 16

# Verified against the data, NOT the summary doc (which lists four celltypes
# that do not exist: Regulatory/exhausted, Squamous epithelium, Endothelial,
# Smooth muscle). Real CLUSTER_ID -> CLUSTER_LABEL mapping:
#   0 APC              6 Macrophage            12 Tumor (CD21+)
#   1 B cell           7 Naive immune cell     13 Tumor (Ki67+)
#   2 CD4 T cell       8 Stromal / Fibroblast  14 Tumor (Podo+)
#   3 CD8 T cell       9 Tumor                 15 Vessel
#   4 Granulocyte     10 Tumor (CD15+)
#   5 Lymph vessel    11 Tumor (CD20+)
IMMUNE = [0, 1, 2, 3, 4, 6, 7]
TUMOR = [9, 10, 11, 12, 13, 14]
STROMA = [5, 8, 15]

FEATURE_NAMES = [
    "kl_mean",       # how far the whole tissue sits from a random mix
    "kl_tumor",      # how distinctive the tumour niche is
    "self_enrich",   # do celltypes clump with their own kind beyond chance
    "immune_tumor",  # immune infiltration beyond chance
    "stroma_tumor",  # stromal/vessel interface with tumour beyond chance
]

# Half-count floor: a block with zero observed contacts would otherwise give
# log2(0) = -inf. 0.5 means "fewer than one edge", keeping the value finite and
# on the same scale as the data.
_FLOOR = 0.5


def _kl_bits(q: np.ndarray, p: np.ndarray) -> float:
    """KL(q || p) in bits, over the support where both are positive.

    q_j > 0 with p_j == 0 cannot occur: if type j has no cells in this sample it
    cannot be anyone's neighbour either.
    """
    ok = (q > 0) & (p > 0)
    if not ok.any():
        return 0.0
    return float(np.sum(q[ok] * np.log2(q[ok] / p[ok])))


def _log2_enrich(num: float, den: float, p_expected: float) -> float:
    """log2(observed fraction / expected fraction). 0 == exactly chance."""
    if den <= 0 or p_expected <= 0:
        return 0.0
    frac = max(num, _FLOOR) / den
    return float(np.log2(frac / p_expected))


def _block_enrich(M: np.ndarray, p: np.ndarray, rows, cols) -> float:
    """Enrichment of celltype-block `cols` among the neighbours of block `rows`.

        P[A->B] = sum(M[A][B]) / sum(M[A][:])
        p_B     = sum(p[B])
        feature = log2( P[A->B] / p_B )
    """
    num = float(M[np.ix_(rows, cols)].sum())
    den = float(M[rows, :].sum())
    return _log2_enrich(num, den, float(p[cols].sum()))


def enrichment_features(M: np.ndarray, counts: np.ndarray) -> np.ndarray:
    """Collapse a 16x16 neighbour-count matrix into the 5 abundance-corrected scalars.

    Parameters
    ----------
    M : (16, 16) RAW neighbour counts (not row-normalised)
    counts : (16,) number of cells of each type in this sample

    Returns
    -------
    (5,) array ordered as FEATURE_NAMES
    """
    total = counts.sum()
    if total <= 0:
        return np.zeros(len(FEATURE_NAMES))
    p = counts.astype(float) / total

    row_sums = M.sum(1, keepdims=True)
    P = np.divide(M, row_sums, out=np.zeros_like(M, dtype=float), where=row_sums > 0)
    valid = row_sums.ravel() > 0

    # 1. kl_mean — mean over celltypes of KL(neighbourhood || global mix).
    #    0 bits = this type sits in a random mix; higher = distinctive niche.
    if valid.any():
        kl_mean = float(np.mean([_kl_bits(P[i], p) for i in np.flatnonzero(valid)]))
    else:
        kl_mean = 0.0

    # 2. kl_tumor — same, for the pooled tumour row (all 6 tumour variants).
    tumor_row = M[TUMOR, :].sum(0)
    ts = tumor_row.sum()
    kl_tumor = _kl_bits(tumor_row / ts, p) if ts > 0 else 0.0

    # 3. self_enrich — mean log2(P[i][i] / p_i): self-segregation beyond chance.
    diag_vals = []
    for i in np.flatnonzero(valid):
        if p[i] > 0:
            diag_vals.append(_log2_enrich(M[i, i], row_sums[i, 0], p[i]))
    self_enrich = float(np.mean(diag_vals)) if diag_vals else 0.0

    # 4/5. block enrichments
    immune_tumor = _block_enrich(M, p, IMMUNE, TUMOR)
    stroma_tumor = _block_enrich(M, p, STROMA, TUMOR)

    return np.array([kl_mean, kl_tumor, self_enrich, immune_tumor, stroma_tumor])


def _M_from_edges(edges: np.ndarray, clusters: np.ndarray) -> np.ndarray:
    """Raw neighbour counts from an undirected edge list (both directions counted).

    Note this matrix is EXACTLY symmetric — measured max|M - M^T| = 0.0 on every
    sample checked. Delaunay adjacency is mutual by construction, so a
    'bi-directional' readout would return identical numbers here; asymmetry only
    exists under directed kNN.
    """
    src = np.concatenate([edges[:, 0], edges[:, 1]])
    dst = np.concatenate([edges[:, 1], edges[:, 0]])
    M = np.zeros((N_CLUSTERS, N_CLUSTERS))
    np.add.at(M, (clusters[src], clusters[dst]), 1)
    return M


def _load_edges(sid: str, graphs_dir: Path) -> np.ndarray:
    """Load prebuilt Delaunay edges -> (E_undirected, 2)."""
    e = np.load(graphs_dir / f"{sid}_edges.npy")
    if e.shape[0] == 2 and e.shape[1] != 2:
        e = e.T
    return np.unique(np.sort(e, axis=1), axis=0)


def build(
    locations_csv: Path,
    graphs_dir: Path,
    out_path: Path,
    limit: Optional[int] = None,
) -> pd.DataFrame:
    """Build the 5-feature table for every sample with a prebuilt Delaunay graph."""
    loc = pd.read_csv(locations_csv,
                      usecols=["ACQUISITION_ID", "X", "Y", "CLUSTER_ID"])
    loc["CLUSTER_ID"] = loc["CLUSTER_ID"].astype(int)

    available = {p.name.replace("_edges.npy", "")
                 for p in graphs_dir.glob("*_edges.npy")}
    samples = sorted(set(loc.ACQUISITION_ID.unique()) & available)
    if limit:
        samples = samples[:limit]

    print(f"  {len(samples)} samples with prebuilt Delaunay graphs")
    print(f"  Groups: IMMUNE={IMMUNE}  TUMOR={TUMOR}  STROMA={STROMA}")

    rows: Dict[str, np.ndarray] = {}
    skipped = []
    for n, sid in enumerate(samples, 1):
        s = loc[loc.ACQUISITION_ID == sid]
        clusters = s["CLUSTER_ID"].values
        edges = _load_edges(sid, graphs_dir)

        if edges.max() >= len(clusters):
            # Edge indices must address this sample's own cell rows.
            skipped.append((sid, f"edge idx {edges.max()} >= {len(clusters)} cells"))
            continue

        M = _M_from_edges(edges, clusters)
        counts = np.bincount(clusters, minlength=N_CLUSTERS)
        rows[sid] = enrichment_features(M, counts)

        if n % 50 == 0:
            print(f"    {n}/{len(samples)}")

    if skipped:
        print(f"  [WARNING] skipped {len(skipped)}: {skipped[:3]}")

    df = pd.DataFrame.from_dict(rows, orient="index", columns=FEATURE_NAMES)
    df.index.name = "acquisition_id"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(out_path)

    print(f"\n  Built: {df.shape[0]} samples x {df.shape[1]} features -> {out_path.name}")
    print(f"\n  Distribution (0 == exactly chance for every feature):")
    print(df.describe().T[["mean", "std", "min", "max"]].to_string())
    return df


if __name__ == "__main__":
    import sys
    sys.path.insert(0, str(Path(__file__).parent.parent))
    from config import get_config

    cfg = get_config()
    build(
        locations_csv=cfg.locations_path,
        graphs_dir=cfg.processed_dir.parent / "step1_delaunay" / "graphs",
        out_path=cfg.processed_dir / "neighbor_features" / "enrichment.parquet",
    )
