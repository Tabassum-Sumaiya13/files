# Spatial PE Pipeline — Full Analysis, Interpretation & Improvement Ideas

## Part 1: How Each Stage Works (Mechanically)

### Stage 1: Preprocessing

```
Raw CSVs → Clean → Merge → Normalize → Per-sample parquets
```

**What happens step by step:**

| Step | Input | Operation | Output | Why |
|---|---|---|---|---|
| Load | 3 CSVs | Read `cell_locations_and_labels.csv` (2M rows × 10 cols), `labeled_arcsinh_norm_data.csv` (2M rows × 44 cols), `marker_names.csv` (39 names) | 2 DataFrames + marker list | Separate files for spatial coords vs protein expression |
| QC Filter | All samples | Keep only 307 samples from `qc_acq_ids_labeled.csv` | 308 samples (one extra matched) | Paper's QC-passing set ensures fair comparison |
| Clean | Both DFs | Drop duplicates, NaN in X/Y/markers, align cell IDs | Same DFs, fewer rows | Prevent graph construction failures |
| Merge | Locations + Expression | Inner join on `(acquisition_id, cell_id)` | Single DF: 2,061,102 rows × 45 cols | Each cell now has coords + markers + labels |
| Filter | Merged | Drop samples with < 50 cells | All 308 pass | Too-small samples can't produce stable eigenvectors |
| Normalize | X, Y per sample | Per-sample z-score: `X = (X - mean) / std` | X, Y ~ N(0,1) per sample | Makes k-NN distances comparable across tissues of different physical sizes |
| Export | Normalized DF | One `.parquet` per sample + manifest | 308 files in `data/processed/samples/` | Enables parallel per-sample encoding |

**Key data at this point:**
- 308 samples, 2,061,102 total cells
- Each cell: `(acquisition_id, cell_id, X, Y, cluster_id, cluster_label, 39 markers)`
- Coordinates are z-scored per sample

---

### Stage 2: Positional Encoding

```
Per-sample parquet → k-NN graph → Laplacian → Eigenvectors → Feature vector
```

**For EACH of the 308 samples:**

#### Step A: Build k-NN Graph
```
Input:  N cells with (X, Y) coordinates
Output: N × N sparse adjacency matrix

Algorithm:
1. For each cell, find its k=10 nearest neighbors by Euclidean distance
2. Create directed edge: cell_i → cell_j for each neighbor
3. Symmetrize with "union": edge exists if EITHER i→j OR j→i
4. Remove self-loops
```

**What this produces:** A graph where each cell is connected to ~10-11 nearby cells. The graph captures the **spatial neighbourhood structure** — which cells are physically close to each other in the tissue.

**Typical stats from our run:**
- Edges per sample: 8,000–105,000 (depending on sample size)
- Mean degree: 10.8–11.8 (slightly above k=10 due to union symmetrization)
- Components: Usually 1 (fully connected), sometimes 2-5 fragments

#### Step B: Compute Laplacian Eigenvectors

```
Input:  N × N adjacency matrix A
Output: N × 8 positional encoding matrix

Algorithm:
1. If graph has multiple components, keep only the largest (LCC)
2. Compute degree matrix: D = diag(sum of each row of A)
3. Build symmetric normalized Laplacian:
       L_sym = I - D^{-1/2} · A · D^{-1/2}
4. Find the 9 smallest eigenvalues/eigenvectors using shift-invert eigsh
5. Discard the trivial eigenvector (eigenvalue ≈ 0, constant vector)
6. Keep the next 8 eigenvectors → these ARE the positional encoding
7. Fix sign convention: make max-absolute entry positive (reproducibility)
```

**What the eigenvectors represent:**
- **Eigenvector 1** (smallest non-trivial eigenvalue): the **smoothest** spatial partition — like cutting the tissue into two halves. High values on one side, low on the other.
- **Eigenvector 2**: the next smoothest partition — orthogonal to ev1. Like cutting into quadrants.
- **Eigenvectors 3-8**: progressively finer-grained spatial patterns.

Think of it like **spatial frequency decomposition**: ev1 = lowest frequency (global structure), ev8 = higher frequency (more local).

**Typical stats from our run:**
- Trivial eigenvalue: ~10^{-17} (essentially 0, good)
- PE eigenvalue range: 0.0003 – 0.016
- Orthonormality deviation: ~10^{-15} (machine precision, excellent)
- 99.7% of cells got valid PE (only 15 samples had disconnected cells)

#### Step C: Assemble Features

```
Input:  PE vectors (N × 8) + Marker values (N × 39)
Output: Feature matrix (N × 47)

Simply: features = [pe_0, ..., pe_7 | CD31, CD57, ..., CD3e]
```

Each cell now has a **47-dimensional feature vector**: 8 spatial encoding dims + 39 protein marker dims.

---

### Stage 3: Survival Evaluation

```
Per-cell features → Pool to per-sample → RSF → C-index
```

#### Step A: Pool Cell Features → Sample Features
```
Input:  N cells × 47 features per sample
Output: 1 vector per sample

For PE (mean_std pooling):
  - Take mean of each PE dim across all cells → 8 values
  - Take std of each PE dim across all cells → 8 values
  - Concatenate → 16 features per sample
```

**This is the critical step** — we collapse ~6,000 cells into 16 numbers. This is a massive information bottleneck.

#### Step B: Train RSF
```
- RandomSurvivalForest(n_estimators=100, random_state=1029)
- GroupKFold(n_splits=10) by patient_id (81 unique patients)
- Metric: concordance_index_censored (C-index)
```

**GroupKFold by patient** prevents leakage: one patient can have multiple tissue samples — if patient A's samples appear in both train and test, the model could cheat by memorizing patient-level patterns.

---

## Part 2: What We Got — Interpretation

### Results Table

| Rank | Feature Set | C-index | Features |
|---|---|---|---|
| 1 | **Celltype Proportions** | **0.619 ± 0.146** | 16 |
| 2 | PE + Celltype | 0.602 ± 0.184 | 32 |
| 3 | PE (Laplacian) | 0.536 ± 0.162 | 16 |
| 4 | Raw Coordinates | 0.526 ± 0.157 | 4 |
| 5 | PE + Raw Coords | 0.515 ± 0.145 | 20 |
| — | Spatsurv baseline (paper) | **0.704** | ~256+ |
| — | Random guessing | 0.500 | — |

### What This Means

#### 1. PE barely beats raw coordinates (0.536 vs 0.526)

The graph Laplacian eigenvectors capture *slightly* more spatial information than raw X,Y — but the gap is tiny (0.01 C-index). This means:

> **The spatial structure encoded by the Laplacian is not strongly predictive of survival on its own.**

This makes biological sense: whether cell A is at coordinate (100, 200) or (300, 400) in the tissue doesn't directly predict patient death. What matters is **what types of cells are near each other** — which is exactly what celltype proportions capture.

#### 2. Celltype proportions dominate (0.619)

The composition of the tissue (what fraction of cells are T-cells, tumor cells, fibroblasts, etc.) is far more predictive than spatial coordinates. This is consistent with the spatsurv paper's finding.

#### 3. PE + Celltype doesn't improve over celltype alone (0.602 < 0.619)

Adding PE to celltype proportions actually **hurts**. This is a classic sign that:
- The PE features add noise without adding signal
- With 32 features and only ~30 samples per test fold, the model overfits
- The pooling method (mean/std) destroys the spatial information that PE was supposed to capture

#### 4. Our celltype proportions (0.619) vs paper's combined (0.704)

The paper gets 0.704 using **celltype proportions + neighbourhood matrices** (which cell types are physically adjacent). The neighbourhood matrix is a 16×16 matrix (256 features) that captures **spatial relationships between cell types** — exactly the kind of information our PE was supposed to encode.

#### 5. Fold 2 always fails

One of the 10 folds has zero events (all censored patients). This wastes 10% of our evaluation and inflates variance. This is a side-effect of having only 81 patients with uneven event rates.

#### 6. Massive variance (±0.15)

Standard deviations of 0.14–0.18 mean the C-index swings between ~0.4 and ~0.8 across folds. With only 30 samples per test fold and 81 patients total, the evaluation is quite noisy.

---

## Part 3: WHY PE Underperforms — Root Cause Analysis

### Problem 1: Information Destruction via Pooling

```
6,000 cells × 8 PE dims  →  mean/std pooling  →  16 numbers
```

The PE encodes **where each cell is** within the tissue. When we average across all cells, we lose all of this. It's like knowing the GPS coordinate of every house in a city, then reducing it to "average latitude" and "average longitude" — you lose all neighbourhood structure.

**The spatsurv neighbourhood matrix doesn't have this problem** because it directly encodes pairwise celltype relationships (e.g., "fraction of T-cell neighbours that are tumor cells").

### Problem 2: PE is Purely Geometric

The Laplacian eigenvectors encode **spatial position on the graph** — but they are blind to **what type of cell** is at each position. Two samples with identical spatial layouts but different celltype arrangements would get identical PE vectors.

The survival-relevant signal is not "where are the cells" but **"which cell types are near which other cell types"** — a joint spatial-functional question that pure PE doesn't answer.

### Problem 3: Sign/Orientation Ambiguity Across Samples

Even with sign-fixing, the Laplacian eigenvectors of different samples are not aligned. PE dim 0 in sample A might encode "left-right" while PE dim 0 in sample B encodes "top-bottom". When we pool across samples and train an RSF, the model can't learn consistent patterns because the same feature means different things in different samples.

### Problem 4: The Pooling-then-Predict Architecture

```
Current:   cells → PE → pool → RSF
Spatsurv:  cells → spatial statistics (neighbourhood matrix) → RSF
```

Spatsurv computes **domain-specific spatial statistics** that are invariant to rotation, translation, and sample size. Our pipeline computes generic spatial coordinates, then throws away the structure by pooling.

---

## Part 4: Ideas to Build a Better Pipeline

### Tier 1: Quick Fixes (modify existing code)

#### Idea 1: Better Pooling — Celltype-Stratified PE Statistics
Instead of pooling PE across ALL cells, pool separately by cell type:

```python
# For each celltype (16 types) × each PE dim (8) → compute mean
# Result: 16 × 8 = 128 features per sample
for celltype in unique_celltypes:
    mask = cluster_labels == celltype
    celltype_pe_mean = pe_vectors[mask].mean(axis=0)  # 8 values
```

**Why this helps:** It captures WHERE each cell type tends to be. "T-cells have high PE_0" means T-cells are concentrated on one side of the tissue. This is closer to what the neighbourhood matrix captures.

#### Idea 2: PE-Based Neighbourhood Matrix
Use PE vectors to compute a **similarity-weighted neighbourhood matrix**:

```python
# For each pair of celltypes (i, j):
#   Compute average PE distance between type-i and type-j cells
# Result: 16 × 16 = 256 features (same shape as spatsurv's neighbor_mat)
```

This directly competes with spatsurv's neighbourhood matrix but uses the continuous PE space instead of discrete k-NN.

#### Idea 3: Increase k_pe
Currently using 8 eigenvectors. Try 16 or 32 — more dimensions capture finer spatial structure. But beware: more features → more pooled features → more overfitting risk.

#### Idea 4: Use Existing Neighbourhood Matrices as Features
The raw data directory already has `celltype_neighborhoods/k10/*.npy` files! Load those directly and combine with PE:

```python
# data/raw/celltype_neighborhoods/k10/{sample_id}.npy
# These are the exact features spatsurv uses for their 0.704 baseline
```

### Tier 2: Architecture Changes

#### Idea 5: Graph Neural Network (GNN) Instead of Eigenvectors
Replace the Laplacian PE → pooling → RSF pipeline with:

```
cells → k-NN graph → GNN (message passing) → graph-level readout → survival head
```

A GNN can learn to aggregate spatial and marker information jointly, without the information-destroying mean/std pooling. Use a simple architecture like GCN or GraphSAGE.

#### Idea 6: Attention-Based Pooling
Instead of mean/std, use a learnable attention mechanism:

```python
# Learn which cells matter most for survival
attention_weights = softmax(MLP(pe_vectors || markers))  # per-cell importance
sample_embedding = sum(attention_weights * features)      # weighted sum
```

This lets the model focus on clinically relevant cells (e.g., immune cells at the tumor margin) instead of averaging everything.

#### Idea 7: Multi-Scale Graphs
Build multiple graphs with different k values (k=5, 10, 20, 50) and compute PE at each scale. Different scales capture different biological phenomena:
- k=5: immediate microenvironment (cell-cell contact)
- k=50: mesoscale tissue regions

#### Idea 8: PE as Node Features for GNN
Use PE as **input features** to a GNN, alongside markers:

```
Input: node features = [pe_0..7 | marker_0..38]  (47 dims per cell)
Model: 2-layer GCN → global mean pooling → MLP → survival prediction
```

The GNN can learn to use PE to distinguish cells that have similar markers but different spatial contexts.

### Tier 3: Fundamentally Different Approaches

#### Idea 9: Contrastive Pre-training
Pre-train the PE encoder with a contrastive objective:
- **Positive pairs**: cells from the same spatial neighbourhood (same PE cluster)
- **Negative pairs**: cells from distant regions
- Then fine-tune the learned representations for survival

#### Idea 10: Spatial Point Process Features
Instead of graph PE, compute statistical spatial features:
- **Ripley's K function** per cell type (already in spatsurv data!)
- **Mark correlation functions**
- **Cross-type nearest-neighbour distances**
- These are well-established in spatial statistics and directly meaningful

#### Idea 11: Cell Spatial Motif Detection
Identify recurring spatial patterns ("motifs"):
- Immune cell clusters around tumor nests
- Tertiary lymphoid structures
- Stromal barriers
Encode these as binary/count features per sample.

#### Idea 12: Combine With Spatsurv Features
The highest-impact short-term improvement:

```python
features = concat([
    celltype_proportions,      # 16 dims (our 0.619 baseline)
    neighbourhood_matrix,      # 256 dims (from data/raw/k10/*.npy)
    celltype_stratified_PE,    # 128 dims (idea 1 above)
])
# Total: 400 features → RSF with feature selection
```

This combines the proven spatsurv features with the novel PE features and lets the RSF decide what's useful.

---

## Part 5: Recommended Next Steps (Priority Order)

| Priority | Action | Expected Impact | Effort |
|---|---|---|---|
| **1** | Load spatsurv neighbourhood matrices from `data/raw/k10/` and add as baseline | Reproduce 0.704 exactly | Low |
| **2** | Implement celltype-stratified PE pooling (Idea 1) | Moderate — tests if PE captures celltype spatial patterns | Low |
| **3** | Implement PE-based neighbourhood matrix (Idea 2) | High — directly comparable to spatsurv's best feature | Medium |
| **4** | Combine all features (Idea 12) with proper feature selection | High — uses everything | Medium |
| **5** | Build a simple GNN baseline (Idea 5) | Potentially high — avoids pooling bottleneck entirely | High |

> [!IMPORTANT]
> **The core insight:** PE itself is not the problem — the **pooling step** is. The Laplacian eigenvectors contain rich spatial information, but collapsing 6,000 cells × 8 dims into 16 numbers via mean/std destroys it. Every improvement idea above is fundamentally about **preserving more spatial structure** during the cell→sample aggregation step.
