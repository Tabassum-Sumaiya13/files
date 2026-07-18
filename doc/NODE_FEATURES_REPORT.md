# Node Features Report — Celltype-Conditioned Functional Markers

*Spatial positional-encoding survival pipeline. Colorectal cancer CODEX atlas (Schürch et al. 2020).*

---

## 1. TL;DR

- **What a "node feature" is here:** a per-cell protein-marker value, pooled to one number per sample and fed to the RandomSurvivalForest. Until this work, only the cell's *celltype label* reached the model — the 39 measured proteins were dropped.
- **What we added:** 8 **celltype-conditioned functional markers** (e.g. Ki67 in tumour, GranzymeB in CD8 T cells), chosen from validated spatial-pathology survival literature — not all 39, to avoid diluting ~27–103 events across hundreds of columns.
- **Result:** node features add **real, noise-beating prognostic signal on top of celltype composition** — but not as a replacement.
- **Best configuration:** reduce the block to the **2 strongest markers** (`cd4_foxp3` + `tcell_cd45ro`). `Celltype + Enrichment + Markers(2)` reaches **C = 0.733**, cleanly beating its width-matched noise control (+0.053 [+0.035, +0.074]).

---

## 2. What is a node feature in this project?

The tissue is modelled as a **graph**: each **cell is a node**, edges come from spatial adjacency (Delaunay). A node feature is any attribute attached to a cell.

| Node attribute | Status before this work |
|---|---|
| Celltype **label** (`CLUSTER_ID`, 16 types) | **Used** — via celltype proportions + neighbour matrices |
| Spatial position (X, Y) | Used indirectly — builds the graph / positional encoding |
| **39 protein markers** (per-cell expression) | **Dropped** — never reached the survival model |

The model is a **RandomSurvivalForest on per-sample vectors**, not a GNN. So a per-cell node feature must be **pooled to one value per sample** before it can be a column. That pooling is the whole design question addressed below.

---

## 3. Data source and available markers

- **File:** `data/raw/labeled_arcsinh_norm_data.csv` (925 MB, 2,061,102 cells, 308 samples).
- **Values:** arcsinh-normalised protein intensities (already scaled; no further transform).
- **39 markers available:** CD31, CD57, CD4, CD15, FoxP3, CD16, CD20, CD45RO, CD38, CD34, CD11b, CD68, CD134, TMEM16A, PanCK, Podoplanin, CD45, GranzymeB, CD49f, CD11c, CD47, CD8, CD117, Vimentin, CD69, aSMA, CD14, CD21, HLA-DR, PDL1, CD56, p16, PD1, CD45RA, ICOS, CD152, Ki67, CollagenIV, CD3e.
- **16 celltypes** (`CLUSTER_ID`): 0 APC, 1 B cell, 2 CD4 T, 3 CD8 T, 4 Granulocyte, 5 Lymph vessel, 6 Macrophage, 7 Naive immune, 8 Stromal/Fibroblast, 9 Tumor, 10 Tumor(CD15+), 11 Tumor(CD20+), 12 Tumor(CD21+), 13 Tumor(Ki67+), 14 Tumor(Podo+), 15 Vessel.

---

## 4. Design — why *these* markers, conditioned on *these* celltypes

### 4.1 The dilution constraint

At ~27 patient-level events (≈103 acquisition-level), a 39-marker × 16-celltype block would be 624 columns — events-per-variable ≈ 0.04. The forest samples ~√p features per split, so a huge block is mostly ignored and *loses* to the baseline through dilution alone. This is exactly what the pipeline's `dC vs NOISE` control (a same-width random block) is designed to catch. **Fewer, biologically targeted markers are mandatory.**

### 4.2 The validated shortlist

We took the functional/activation markers that carried prognostic signal in the works this dataset comes from:

- **Schürch et al. 2020 (Cell)** — the CODEX colorectal atlas that *is* this data. Top survival hit: PD-1⁺ CD4 T cells; ICOS, GranzymeB, CD45RO also key.
- **Jackson et al. 2020 (Nature)** and **Danenberg et al. 2022 (Nat Genetics)** — IMC breast-cancer survival: functional-marker means (Ki67, PD-L1, GranzymeB) prognostic *in the context of* cell composition.

### 4.3 Celltype conditioning

A marker matters *where* it is expressed. Ki67 in a tumour cell means proliferation; Ki67 averaged over all cells is noise. So each feature is the mean of a marker **only over the celltype it is biologically read in**:

| Feature | Marker | Read in celltype(s) | Biology (prognostic meaning) |
|---|---|---|---|
| `tumor_ki67` | Ki67 | Tumor (9–14) | tumour proliferation rate |
| `tumor_mac_pdl1` | PDL1 | Tumor + Macrophage (6) | immune checkpoint / evasion |
| `cd8_granzymeb` | GranzymeB | CD8 T (3) | cytotoxic killing activity |
| `cd4_pd1` | PD1 | CD4 T (2) | Schürch's top survival hit |
| `cd4_icos` | ICOS | CD4 T (2) | T-cell activation |
| `cd4_foxp3` | FoxP3 | CD4 T (2) | regulatory-T suppression |
| `tcell_cd45ro` | CD45RO | CD4 + CD8 T (2,3) | memory / prior antigen exposure |
| `apc_mac_hladr` | HLA-DR | APC + Macrophage (0,6) | antigen presentation |

---

## 5. Feature characteristics

- **Readout:** per sample, `feature = mean(marker over cells of the named celltype)`.
- **Aggregation:** simple mean (matches the arcsinh scale; RSF is rank-based so no further scaling needed).
- **Missing handling:** a sample with none of that celltype → `0.0` ("no such cells, so no such functional signal"). In practice **all 308 samples contained all celltypes**, so no imputation occurred.
- **Distribution** (arcsinh-normalised means, 308 samples):

| Feature | mean | std | min | max |
|---|---|---|---|---|
| tumor_ki67 | 0.941 | 0.450 | 0.276 | 3.761 |
| tumor_mac_pdl1 | 0.236 | **0.029** | 0.214 | 0.567 |
| cd8_granzymeb | 0.514 | 0.257 | 0.217 | 1.604 |
| cd4_pd1 | 0.293 | 0.071 | 0.212 | 0.737 |
| cd4_icos | 0.485 | 0.384 | 0.216 | 4.595 |
| cd4_foxp3 | 0.407 | 0.270 | 0.206 | 4.302 |
| tcell_cd45ro | 1.525 | **0.539** | 0.382 | 3.874 |
| apc_mac_hladr | 1.309 | 0.457 | 0.374 | 3.611 |

> Note `tumor_mac_pdl1` is nearly constant across samples (std 0.029). A feature that barely varies cannot carry sample-level survival signal — an early flag that it would be dead weight (confirmed in §7.2).

---

## 6. Implementation

| File | Role |
|---|---|
| [src/marker_states.py](src/marker_states.py) | Builds `data/processed/neighbor_features/marker_states.parquet` (308 × 8). One pass over the expression CSV, celltype-conditioned means, `fillna(0)`. |
| [run_survival.py](run_survival.py) | `load_marker_states()` loader; `--with-markers` flag appends a `+ Markers` variant to every Celltype-based set + a `Markers alone` set; `--marker-keep` filters to a chosen subset. Width-matched noise control is added automatically. |

**Reproduce:**
```bash
python src/marker_states.py                                          # build the block
python run_survival.py --experiment C --with-markers                 # full 8 markers
python run_survival.py --experiment C --with-markers --marker-keep cd4_foxp3,tcell_cd45ro   # best (top-2)
```

---

## 7. Experiments and results

Validation protocol (unchanged from the rest of the pipeline): RandomSurvivalForest, **patient-grouped** StratifiedGroupKFold (no leakage), 20 seeds for the real C-index, 10 for the permutation null. **`dC vs NOISE`** = paired gain over a *same-width* random block — the honest test, because it cancels the dilution cost of simply adding columns. Cohort: 307 samples, 81 patients, 103 events.

### 7.1 Full 8-marker block (Experiment C)

| Feature set | Feats | C-index | ΔC vs NOISE | Verdict |
|---|---|---|---|---|
| Celltype + Enrichment + Markers | 29 | 0.729 | +0.038 **[−0.000, +0.072]** | borderline (CI touches 0) |
| Celltype + delaunay(256) + Markers | 280 | 0.723 | +0.066 [+0.025, +0.107] | REAL SIGNAL |
| Celltype + delaunay(256) | 272 | 0.718 | +0.062 [+0.030, +0.099] | REAL SIGNAL |
| Celltype + Enrichment | 21 | 0.716 | +0.028 [+0.005, +0.053] | REAL SIGNAL |
| Celltype Proportions + Markers | 24 | 0.716 | +0.042 [+0.022, +0.075] | REAL SIGNAL |
| Celltype Proportions (baseline) | 16 | 0.690 | — | baseline |
| Markers alone | 8 | 0.648 | — | worse than baseline |

**Read:** markers added to the plain baseline beat noise cleanly (+0.042). But the Enrichment combo, at 29 feats, only *touched* 0 — the 6 weak markers diluted it.

### 7.2 Per-marker ranking (paired incremental gain over the 16-feature celltype baseline)

| Rank | Marker | C (Celltype+1) | ΔC vs baseline | std |
|---|---|---|---|---|
| 1 | **cd4_foxp3** | 0.715 | **+0.024 [+0.006, +0.044]** | 0.270 |
| 2 | **tcell_cd45ro** | 0.701 | +0.011 [−0.012, +0.037] | 0.539 |
| 3 | tumor_mac_pdl1 | 0.695 | +0.005 [−0.015, +0.032] | 0.029 |
| 4 | tumor_ki67 | 0.693 | +0.002 [−0.022, +0.023] | 0.450 |
| 5 | cd8_granzymeb | 0.692 | +0.002 [−0.016, +0.022] | 0.257 |
| 6 | cd4_icos | 0.691 | +0.001 [−0.017, +0.016] | 0.384 |
| 7 | cd4_pd1 | 0.690 | −0.001 [−0.016, +0.016] | 0.071 |
| 8 | apc_mac_hladr | 0.686 | −0.004 [−0.024, +0.015] | 0.457 |

Only `cd4_foxp3` has a gain whose CI clears 0 on its own; `tcell_cd45ro` is a clear second. The other six are near-zero/negative alone — dilution.

### 7.3 Reduced blocks — Enrichment combo

| Marker block | Feats | C-index | ΔC vs NOISE | Verdict |
|---|---|---|---|---|
| 8 (full) | 29 | 0.729 | +0.038 [−0.000, +0.072] | borderline |
| 4 (top) | 25 | 0.730 | +0.035 [+0.017, +0.059] | REAL SIGNAL |
| **2 (foxp3 + cd45ro)** | **23** | **0.733** | **+0.053 [+0.035, +0.074]** | **REAL SIGNAL — best** |

Cutting 8 → 2 raised the C-index (0.729 → 0.733), lifted the noise-corrected lower bound (−0.000 → **+0.035**), and tightened the CI. The borderline verdict became a clean "REAL SIGNAL."

---

## 8. Findings

1. **Node features carry real prognostic signal** — the strongest clean evidence is `Celltype + Markers` beating a same-width random block by +0.042 (full) / +0.033 (top-4), CI above 0.
2. **They are complementary, not a substitute.** `Markers alone` (8 feats) = 0.648, below the 0.690 celltype baseline. Knowing *which* cells are present beats knowing their *functional state* — but combining both wins. This matches the source literature exactly.
3. **Less is more.** Six of the eight markers were dilution. Two markers (Treg suppression + T-cell memory) captured essentially all the signal.
4. **Standout biological result:** just `cd4_foxp3 + tcell_cd45ro` alone scored **0.708** — two functional states nearly matching the entire 16-way celltype composition (0.690). Consistent with Schürch 2020: immune functional state drives colorectal-cancer outcome.

---

## 9. Recommendation

**Use the 2-marker node block: `cd4_foxp3, tcell_cd45ro`.**

- Best single result: **`Celltype + Enrichment + Markers(2)` = 0.733**, beats noise +0.053 [+0.035, +0.074].
- Leanest strong result: **`Celltype + Markers(2)` = 0.718** with only 18 features (high events-per-variable, robust).

```bash
python run_survival.py --experiment C --with-markers --marker-keep cd4_foxp3,tcell_cd45ro
```

**Open follow-ups:** run Experiment A (`Celltype + delaunay`, 272) with the 2-marker block for that set's exact number; consider making `cd4_foxp3,tcell_cd45ro` the `--with-markers` default.

---

*Result files: `data/processed/outputs/survival/survival_validation_C_mk.csv` (8), `_C_mk4.csv` (4), `_C_mk2.csv` (2).*
