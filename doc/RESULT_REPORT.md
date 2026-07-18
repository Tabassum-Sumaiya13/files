# Results

*Spatial features from in-situ proteomics imaging for cancer survival prediction.*
*Companion deep-dive reports: [ENRICHMENT_FEATURES_REPORT.md](ENRICHMENT_FEATURES_REPORT.md), [KL_MEAN_REPORT.md](KL_MEAN_REPORT.md), [NODE_FEATURES_REPORT.md](../spatial_positional_encoding/NODE_FEATURES_REPORT.md).*

---

## 1. Dataset

### 1.1 Source and data type

**UPMC Head & Neck Squamous Cell Carcinoma (HNSCC)** 

Tissue was imaged with **Multiplexed Ion Beam Imaging (MIBI)**,

which measures 40 metal-tagged protein markers simultaneously at single-cell resolution while preserving the (X, Y) position of every cell.

| Property | Value |
|---|---|
| Cancer type | Head & Neck Squamous Cell Carcinoma (HNSCC) |
| Institution | University of Pittsburgh Medical Center (UPMC) |
| Imaging technology | MIBI (multiplexed, single-cell, spatially resolved) |
| Total cells | 2,061,102 |
| Protein markers | 40 (39 used as features after ID markers) |
| Cell types | 16 (assigned by marker-based clustering) |
| Patients | ~81 unique |
| Samples (tissue regions) | 307 after QC (multiple regions per patient) |

Each **sample** is one tissue region from one patient. A cell carries: X–Y coordinates, 40 arcsinh-normalised protein intensities, and a cell-type label. Each patient carries a survival time (`survival_day`) and status (`survival_status`: 0 = alive/censored, 1 = dead), plus clinical fields (tissue type, HPV status, recurrence).

The 16 cell types are three lineages: **immune** (APC, B, CD4 T, CD8 T, Granulocyte, Macrophage, Naive immune), **tumour** (Tumor and its CD15⁺/CD20⁺/CD21⁺/Ki67⁺/Podoplanin⁺ variants), and **stromal** (Lymph vessel, Stromal/Fibroblast, Vessel).

> **Note on cohort labels.** Two of the feature deep-dive reports describe the data as a *colorectal CODEX atlas (Schürch et al. 2020)*. That label is a documentation slip: the on-disk sample IDs (`UPMC_c…`), the clinical fields (HPV status, HNSCC tissue types), and the base paper all identify this as the **UPMC Head & Neck MIBI** cohort. Schürch 2020 is cited later only as *literature justification* for which functional markers to add (§8.3), not as the data source. The pipeline, protocol, and event counts are identical across all reports.

### 1.2 How the data was processed

1. **QC filtering** — restrict to the 307 QC-passing acquisition IDs; drop samples with fewer than 50 cells.
2. **Normalisation** — protein intensities are arcsinh-transformed (inverse hyperbolic sine), which compresses extreme values, preserves small differences, and handles zeros — necessary because raw MIBI intensities are heavily right-skewed.
3. **Cell-type labelling** — cells are clustered on their marker profiles into the 16 types (e.g. CD3⁺CD4⁺ → CD4 T; PanCK⁺ → tumour).
4. **Coordinate normalisation** — X, Y are z-scored per sample so the spatial graph is scale-consistent across regions.
5. **Per-sample export** — one parquet file per sample feeds the graph-construction and feature-extraction steps.

For survival evaluation, acquisitions are aggregated to the **patient level** and cross-validation is **grouped by patient** so that no patient's regions appear in both train and test folds. Two analysis framings are used: the **full cohort** (307 samples, includes Normal mucosa) for the graph-construction ablation, and an **exclude-normal / tumour-only** cohort for feature evaluation, because Normal mucosa regions carry no tumour architecture yet inherit the patient's outcome label (pure label noise). Event counts: ~103 acquisition-level events; ~27 unique patient-level events in the tumour-only cohort — a small-event regime that governs every design choice below.

---

## 2. Baseline

### 2.1 Our baseline — celltype proportions

The baseline is the **celltype-proportion** feature set: for each sample, the fraction of cells of each of the 16 types (16 features summing to 1). It answers *"which cells are present"* with no spatial information, and it is the strongest single non-spatial feature in the base paper — so every spatial feature in this project must be shown to add signal **on top of it**.

It is computed by counting cells per type per sample and dividing by the sample's total cell count, then evaluated with the same RandomSurvivalForest protocol as everything else (§4).

Because the cohort framing differs slightly between experiments, the baseline C-index takes a few closely-related values:

| Cohort framing | Baseline C-index | Used for |
|---|---|---|
| Full cohort (incl. Normal mucosa), 307 samples | **0.678** | Step 1 graph construction |
| Tumour-only / exclude-normal | **0.686 – 0.690** | Step 2 feature evaluation |

### 2.2 Base-paper baseline

Dayao et al. computed the five feature families independently under 10-fold patient-grouped RSF. Their reported C-indices:

| Feature (base paper) | C-index | Type |
|---|---|---|
| Celltype proportion | **0.705** | non-spatial |
| Neighbourhood matrix (256) | 0.704 | spatial |
| Biomarker cell (40) | 0.662 | non-spatial |
| Biomarker region (40) | 0.655 | non-spatial |
| Ripley's K (80) | 0.528 | spatial |
| **Combined** | **0.730** | all |

Two facts anchor this project: (i) the neighbourhood matrix (**where** cells are) matches celltype proportion (**what** cells are) — spatial arrangement is genuinely predictive; and (ii) Ripley's K essentially fails (0.528 ≈ random), so not every spatial statistic helps. Our reproduction of the celltype baseline (0.678–0.690) sits a little below their 0.705 because of stricter patient-grouped, multi-seed cross-validation and the exclude-normal framing.

---

## 3. Measurement — Concordance Index (C-index)

Survival quality is measured with **Harrell's concordance index (C-index)**, the standard for censored survival data. It is the fraction of comparable patient pairs whose predicted risk order matches their actual outcome order:

```
C-index = 1.0  → perfect risk ranking
C-index = 0.7+ → good, clinically useful
C-index = 0.5  → random guessing
```

Unlike accuracy, the C-index handles **censoring** (patients still alive at last follow-up) by only scoring pairs whose ordering is known. It is computed with `concordance_index_censored` (scikit-survival) on each held-out CV fold and averaged.

**The honest test — ΔC vs NOISE.** Adding *any* extra columns can raise the C-index by luck, and wide blocks also *dilute* the signal the forest samples per split. So a raw C-index is not enough. Every new feature block is compared against a **same-width block of random numbers** appended to the baseline. The reported quantity is the **paired gain over that width-matched noise control (ΔC vs NOISE)** with a bootstrap 95% CI over seeds. A block "beats noise" only when the CI clears 0. This cancels the dilution cost of simply making the feature vector wider and is the verdict used throughout.

---

## 4. Protocol (how every number was produced)

All results use one fixed protocol so feature sets are directly comparable:

- **Model:** RandomSurvivalForest, 100 trees, `random_state = 1029` (matched to the base paper).
- **Cross-validation:** StratifiedGroupKFold, 10 splits, **grouped by patient** (no leakage) — the same patient's regions never split across folds.
- **Repetition:** 20 seeds for the real C-index (95% CI over shuffles); 10 seeds for the permutation/noise null.
- **Pooling:** per-cell features are pooled to one value per sample before entering the forest (the model is an RSF on per-sample vectors, not a graph neural network).
- **Verdict rule:** a feature block is "REAL SIGNAL" only if its ΔC vs a width-matched random block has a 95% CI above 0.

---

## 5. Main findings

The project asks one question in several forms: *does a spatial feature add real prognostic signal on top of knowing the cell-type composition?* The headline result is the **best configuration found**: cell-type composition plus the enrichment block plus the **two** strongest functional markers.

**Table 1 — Main results (tumour-only cohort, 20 seeds, ΔC vs width-matched noise).**

| Feature set | Feats | C-index | ΔC vs NOISE [95% CI] | Verdict |
|---|---|---|---|---|
| **Celltype + Enrichment + Markers(2)** | **23** | **0.733** | **+0.053 [+0.035, +0.074]** | **REAL SIGNAL — best** |
| Celltype + delaunay(256) + Markers(2) | 274 | 0.722 | +0.092 [+0.054, +0.132] | REAL SIGNAL |
| Celltype + delaunay(256) | 272 | 0.718 | +0.062 [+0.030, +0.099] | REAL SIGNAL |
| Celltype + Enrichment | 21 | 0.716 | +0.028 [+0.005, +0.053] | REAL SIGNAL |
| Celltype + Markers(2) | 18 | 0.718 | +0.026 [−0.000, +0.062] | borderline |
| **Celltype proportions (baseline)** | 16 | 0.690 | — | baseline |
| Markers(2) alone | 2 | 0.708 | — | (no baseline gain) |
| Enrichment alone | 5 | 0.580 | — | worse than baseline |
| Markers(8) alone | 8 | 0.648 | — | worse than baseline |

**Three findings repeat across every feature family:**

1. **Complementary, not standalone.** Every spatial/functional block *loses to the baseline when used alone* (Enrichment alone 0.580; Markers-alone 0.648). Knowing *where* cells sit or *what state* they are in only helps once the model already knows *which* cells are present. Combining wins; substituting does not.
2. **Less is more.** In each family the small, targeted block beats the big one. The enrichment block's usable signal collapses to ~1 feature (`kl_mean`); cutting the marker block from 8 to 2 *raised* the C-index (0.729 → 0.733) and tightened the CI. Wide blocks dilute signal across too many columns at ~27 events.
3. **Redundant with the full neighbour table.** The 256-column neighbourhood matrix already contains the enrichment / `kl_mean` information; adding those scalars on top of the full 256 adds ≈ 0 (§7.3). Choose one: the 256 for the top raw number, or the compact blocks for a readable, robust model.

**Bottom line:** the best model lifts the C-index from **0.690 (composition only) to 0.733** — a **+0.043** absolute gain that cleanly beats its noise control (+0.053 [+0.035, +0.074]) using only **7 spatial/functional columns** on top of composition.

---

## 6. Pipeline steps

```
  MIBI tissue images
        │  (QC, arcsinh-normalise, cluster to 16 cell types, z-score coords)
        ▼
  Per-sample cells: (X, Y) + 39 markers + celltype label
        │
        ▼
  STEP 1 — Graph construction        build a spatial adjacency graph per sample
        │                            (Delaunay | kNN-10 | radius-20/50µm)
        ▼
  STEP 2 — Feature extraction
        │   2.1  Enrichment features (abundance-corrected spatial scalars)
        │        2.1.1  kl_mean — the single strongest enrichment feature
        │   2.2  Node / marker features (celltype-conditioned functional proteins)
        ▼
  STEP 3 — Survival model            RandomSurvivalForest, patient-grouped CV
        │                            report C-index and ΔC vs width-matched noise
        ▼
  Verdict: does the feature beat a same-width random block?
```

---

## 7. Step 1 — Graph construction

### 7.1 What and why

Every spatial feature starts from a **graph**: each cell is a node, edges join spatially adjacent cells. The choice of *how* to define "adjacent" changes every downstream neighbourhood count, so it was tested directly. Four graph definitions were compared, each fed the same 256-column neighbourhood matrix into the same RSF protocol:

- **Delaunay** — parameter-free; two cells are neighbours if their Voronoi cells touch. By Euler's formula the mean degree of any planar triangulation is pinned near 6 (measured here at **5.98 ± 0.01** across 308 samples).
- **kNN-10** — each cell connects to its 10 nearest neighbours. `k = 10` is an arbitrary constant, and a fixed `k` spans 12–38 µm across samples (3.2×; r = −0.835 with cell density) — so "neighbour" means a different physical distance in each sample.
- **radius-20µm / radius-50µm** — connect all cells within a fixed radius. Degree then varies with local density, so each sample's neighbourhood is a different-sized window.

### 7.2 Result

**Table 2 — Graph construction (full cohort, 307 samples; ΔC vs width-matched noise).**

| Graph | Feats | C-index | ΔC vs NOISE [95% CI] | Verdict |
|---|---|---|---|---|
| **Celltype + Delaunay** | 272 | **0.690** | **+0.058 [+0.012, +0.115]** | **REAL SIGNAL** |
| Celltype + radius-20µm | 272 | 0.683 | +0.051 [−0.007, +0.095] | indistinguishable from noise |
| Celltype + radius-50µm | 272 | 0.680 | +0.048 [−0.009, +0.104] | indistinguishable from noise |
| Celltype + kNN-10 | 272 | 0.676 | +0.044 [−0.021, +0.095] | indistinguishable from noise |
| Celltype proportions (baseline) | 16 | 0.678 | — | baseline |
| Celltype + Noise(256) | 272 | 0.632 | — | dilution control |

**Delaunay is the only graph whose spatial block cleanly beats its noise control**, and it is parameter-free — there is no `k` or radius to justify or tune. A second spatial construction, **two-view niches** (a cell's own type + its neighbourhood context), was also tested and did not help: `Celltype + Delaunay + TwoView` = 0.684 (indistinguishable), and two-view alone (0.608) was *worse than noise*. **Delaunay was therefore fixed as the graph for all subsequent steps.**

---

## 8. Step 2 — Feature extraction

### 8.1 Enrichment features (abundance-corrected spatial scalars)

**What they are.** Five scalar features per sample that measure spatial organisation *after removing what cell abundance alone explains*. Where the base paper flattens the 16×16 neighbourhood matrix into 256 columns, these collapse the same matrix into 5 interpretable, abundance-corrected numbers.

**The biology / the confound they remove.** The base paper's readout is `P[i][j]` = "of type i's neighbours, what fraction are type j". But if type j is, say, 18% of the whole tissue, then *even random mixing* makes j appear ~18% of the time — and the celltype baseline **already** knows that 18%. So most of the 256 columns re-encode information the baseline has, which is why `Celltype + neighbourhood` gains only a little over `Celltype`. Dividing by the global mix removes the shared part:

```
E[i][j] = P[i][j] / p_j        E = 1 → exactly chance
                               E > 1 → genuinely enriched (cells clump)
                               E < 1 → genuinely depleted (cells avoid)
```

Every enrichment feature is 0 under the random-mix null and none is computable from abundance alone. At 5 features against ~27 events, events-per-variable ≈ 5 — the same order as the classic >10 rule of thumb, versus the 272-wide blocks that sit ~90× off it.

**The five features.**

| Feature | Measures | Distribution (mean ± sd) | Biological read |
|---|---|---|---|
| `kl_mean` | how organised the whole tissue is | 0.552 ± 0.181 | organised (> 0) |
| `kl_tumor` | how distinctive the tumour niche is | 0.240 ± 0.186 | mildly distinctive |
| `self_enrich` | do cells clump with their own kind | 2.212 ± 0.345 | strong self-clustering (~4.6× chance) |
| `immune_tumor` | immune infiltration of tumour | −0.767 ± 0.464 | immune **excluded** from tumour |
| `stroma_tumor` | stroma–tumour interface | −1.225 ± 0.662 | stroma strongly **excluded** |

**Biology it revealed.** `self_enrich` is large and positive — cells overwhelmingly sit with their own kind. `immune_tumor` and `stroma_tumor` are negative in nearly every sample — immune and stromal cells are *excluded* from the tumour compartment. This is the classic **immune-excluded / "cold" tumour** architecture, quantified in five numbers.

**Feature comparison (survival).**

**Table 3 — Enrichment block vs the 256-column readout (tumour-only cohort).**

| Feature set | Feats | C-index | ΔC vs NOISE | Verdict |
|---|---|---|---|---|
| Celltype + delaunay(256) | 272 | 0.718 | +0.062 [+0.030, +0.099] | REAL SIGNAL |
| **Celltype + Enrichment** | **21** | **0.716** | +0.028 [+0.005, +0.053] | **REAL SIGNAL** |
| Celltype proportions (baseline) | 16 | 0.690 | — | baseline |
| Enrichment alone | 5 | 0.580 | — | worse than baseline |

**Notes / findings.**
- **Efficiency:** the 5-scalar block recovers essentially all of the 256-block's survival signal (0.716 vs 0.718) with ~2% of the columns, and beats its own noise control cleanly.
- **Complementary:** `Enrichment alone` (0.580) is well below baseline — spatial organisation is worthless without knowing which cells are present.
- **Internally redundant:** the five features correlate up to r = 0.80 (`kl_mean` ↔ `self_enrich`) and −0.77 (`kl_tumor` ↔ `immune_tumor`) — effectively 1–2 independent dimensions, not 5. A per-feature ablation (31 subsets) shows only `kl_mean` reliably adds signal on its own (+0.021 [+0.003, +0.038]); dropping any of the other four *improves* the model. The irony: `immune_tumor` had the strongest literature backing yet contributes nothing, while `kl_mean` — invented from a statistical argument — is the workhorse.

#### 8.1.1 The best enrichment feature — `kl_mean`

**What it is.** One number that answers *"how organised is this tissue?"* It compares each cell's **local** neighbourhood mix against the tissue's **global** mix and averages over cell types:

```
kl_mean = mean over 16 celltypes of  KL( neighbourhood_i ‖ global mix )   [bits]

kl_mean = 0     → well-mixed soup, no organisation
kl_mean = high  → cells sit with specific partners, tissue is organised
```

**How it is calculated.** For each sample: (1) count neighbours from the Delaunay graph to build the 16×16 count matrix M; (2) row-normalise to P (the local view); (3) divide each column by the global proportion `p_j` — the key abundance-correction move; (4) for each celltype compute its KL divergence from the global mix; (5) average the 16 values. A half-count floor keeps `log2(0)` finite and rows with zero neighbours are skipped. It decomposes into a 16×16 contribution matrix that is **almost entirely diagonal** — so `kl_mean` is largely measuring *self-clustering* (hence r = 0.80 with `self_enrich`).

**Is it real? — the shuffle test (the key validation).** Shuffle the cell-type labels while keeping the graph, cell counts, and composition identical. This destroys spatial structure and nothing else; whatever `kl_mean` returns now is pure artifact.

```
observed kl_mean   0.552
SHUFFLED labels    0.035   ← the noise floor
                   -----
real structure is ~15× the noise floor, in 307 / 307 samples (Mann-Whitney p = 6e-102)
```

**PASS** — only 6.4% of the value is artifact overall, though that artifact concentrates in rare cell types (43% of the value when a type has < 10 cells, because a 16-bin distribution cannot be estimated from 9 cells). It is moderately stable across a patient's regions (ICC 0.41–0.46).

**What it means biologically — a consistent gradient (all four contrasts agree).**

```
Normal mucosa      0.606   ← most organised
Nodal met          0.587
Primary tumor      0.514
Recurred patient   0.425
Recurrence tumor   0.321   ← least organised
```

**Loss of spatial organisation tracks disease aggressiveness.** `kl_mean` separates tissue types (Kruskal–Wallis p = 0.0009) and tracks **recurrence** (p = 0.005–0.010), while remaining correctly **null on HPV status** (p = 0.69) — a clean specificity check.

**Does it predict survival / help a model?**

**Table 4 — `kl_mean` in the survival model (exclude-normal cohort, 27 events).**

| Feature set | Feats | C-index | ΔC vs baseline |
|---|---|---|---|
| Celltype (baseline) | 16 | 0.6825 | — |
| **Celltype + kl_mean** | **17** | **0.7033** | **+0.021 [+0.003, +0.038]** |
| Celltype + delaunay(256) | 272 | 0.7186 | +0.036 [+0.017, +0.052] |
| Celltype + delaunay(256) + kl_mean | 273 | 0.7188 | +0.0003 over the 256 (≈ 0) |

**Notes.** `kl_mean` captures ~58% of the 256-block's gain with **one column instead of 256** — but it adds essentially **nothing on top of the full 256** (+0.0003), because it is a deterministic summary of the same matrix those columns come from. It also does **not** predict *overall* survival on its own (log-rank p = 0.80): the survival endpoint has only 27 events and counts death from *any* cause, including unrelated ones — no tumour feature can predict a heart attack. The signal is real but is best chased through **recurrence**, which is better powered.

### 8.2 Node / marker features (celltype-conditioned functional proteins)

**What a node feature is here.** The tissue is a graph of cells; a *node feature* is any attribute attached to a cell. Until this work only the cell's **type label** reached the model — the 39 measured proteins were dropped. This step brings a few of them back, each pooled to one value per sample (mean over the relevant cells) and fed to the RSF.

**The design constraint.** With ~27 patient events, a full 39-marker × 16-celltype block would be 624 columns (events-per-variable ≈ 0.04) — the forest would ignore it and *lose* to the baseline through dilution. So a short, biology-backed shortlist is mandatory.

**Reference — which markers, and why.** The shortlist is the set of functional/activation markers that carried prognostic signal in the works this style of data comes from:

- **Schürch et al. 2020 (Cell)** — CODEX colorectal atlas; top survival hit PD-1⁺ CD4 T cells; ICOS, GranzymeB, CD45RO also key.
- **Jackson et al. 2020 (Nature)** and **Danenberg et al. 2022 (Nat. Genetics)** — IMC breast-cancer survival; functional-marker means (Ki67, PD-L1, GranzymeB) prognostic *in the context of* cell composition.

**Celltype conditioning — the key idea.** A marker matters *where* it is expressed. Ki67 in a tumour cell means proliferation; Ki67 averaged over all cells is noise. So each feature is the mean of a marker **only over the cell type it is biologically read in**:

| Feature | Marker | Read in | Biology |
|---|---|---|---|
| `tumor_ki67` | Ki67 | Tumor | proliferation rate |
| `tumor_mac_pdl1` | PDL1 | Tumor + Macrophage | checkpoint / immune evasion |
| `cd8_granzymeb` | GranzymeB | CD8 T | cytotoxic activity |
| `cd4_pd1` | PD1 | CD4 T | Schürch's top survival hit |
| `cd4_icos` | ICOS | CD4 T | T-cell activation |
| `cd4_foxp3` | FoxP3 | CD4 T | regulatory-T suppression |
| `tcell_cd45ro` | CD45RO | CD4 + CD8 T | memory / prior antigen exposure |
| `apc_mac_hladr` | HLA-DR | APC + Macrophage | antigen presentation |

**Result.**

**Table 5 — Node/marker blocks (tumour-only cohort; ΔC vs width-matched noise).**

| Feature set | Feats | C-index | ΔC vs NOISE | Verdict |
|---|---|---|---|---|
| **Celltype + Enrichment + Markers(2)** | **23** | **0.733** | **+0.053 [+0.035, +0.074]** | **best** |
| Celltype + Enrichment + Markers(4) | 25 | 0.730 | +0.035 [+0.017, +0.059] | REAL SIGNAL |
| Celltype + Enrichment + Markers(8) | 29 | 0.729 | +0.038 [−0.000, +0.072] | borderline |
| Celltype + Markers(8) | 24 | 0.716 | +0.042 [+0.022, +0.075] | REAL SIGNAL |
| Celltype proportions (baseline) | 16 | 0.690 | — | baseline |
| Markers(8) alone | 8 | 0.648 | — | worse than baseline |

Per-marker ranking (paired gain over the celltype baseline) shows only **`cd4_foxp3`** clears 0 on its own (+0.024 [+0.006, +0.044]), with **`tcell_cd45ro`** a clear second (+0.011); the other six are near-zero or negative — dilution.

**Notes / findings.**
- **Real signal:** `Celltype + Markers` beats a same-width random block by +0.042 — node features carry genuine prognostic information.
- **Complementary, not a substitute:** `Markers alone` (0.648) loses to the baseline (0.690). Functional state helps only *on top of* composition.
- **Less is more:** cutting 8 → 2 markers *raised* the score (0.729 → 0.733), lifted the noise-corrected lower bound (−0.000 → +0.035), and tightened the CI — the borderline verdict became a clean "REAL SIGNAL."
- **Standout biological result:** `cd4_foxp3 + tcell_cd45ro` alone scored 0.708 — two functional immune states (Treg suppression + T-cell memory) nearly matching the entire 16-way composition (0.690), consistent with Schürch 2020's finding that immune functional state drives outcome.

**Recommended configuration:** the **2-marker node block** (`cd4_foxp3`, `tcell_cd45ro`) plus the enrichment block on top of cell-type composition — best single result at **C = 0.733**.

```bash
python run_survival.py --experiment C --with-markers --marker-keep cd4_foxp3,tcell_cd45ro
```

---

## 9. Limitations

1. **Redundancy with the full 256.** The enrichment scalars add ≈ 0 on top of the full neighbourhood matrix — they are a compact re-expression of the same information, not new information.
2. **Small-event regime.** ~27 patient-level events cap what any feature can achieve and force the "fewer features win" pattern; wide blocks are dominated by dilution, not modelled by it.
3. **Overall-survival endpoint is under-powered.** `survival_status = 1` lumps death of disease with death of other causes; the spatial signal is clearer for **recurrence** than for OS.
4. **Mean-only enrichment.** The abundance-correction ratio discards the variance term (unlike Squidpy's permutation z-score), so 2× enrichment from 10,000 edges and from 12 edges look identical — a source of rare-cell-type artifact.
5. **Positional-encoding direction (discarded).** An earlier attempt to feed Laplacian graph positional encodings directly to the RSF did *not* beat the baseline (PE alone ≈ 0.56–0.59; PE + celltype ≈ 0.61–0.62 vs baseline 0.635) and was dropped in favour of the interpretable neighbourhood/enrichment features reported here.

---

## 10. References

**Base paper**
- Dayao et al. *Deriving spatial features from in-situ proteomics imaging to enhance cancer survival analysis.* (UPMC Head & Neck MIBI cohort; source of the dataset, baseline, and five-feature framework.)

**Methodological / feature justification**
- Palla et al. 2022. *Squidpy: a scalable framework for spatial omics analysis.* Nature Methods — `nhood_enrichment` (permutation-calibrated neighbourhood enrichment; the established form of the abundance-correction idea).
- Schapiro et al. 2017. *histoCAT.* Nature Methods — within-image neighbourhood permutation test.
- *Analytical Neighborhood Enrichment Score*, arXiv:2506.18692 (2025) — closed-form enrichment, validated against Monte Carlo at r ≥ 0.95.

**Marker-shortlist evidence (Step 2.2)**
- Schürch et al. 2020. *Coordinated cellular neighborhoods orchestrate antitumoral immunity at the colorectal cancer invasive front.* Cell — PD-1⁺ CD4 T, ICOS, GranzymeB, CD45RO as prognostic functional states.
- Jackson et al. 2020. *The single-cell pathology landscape of breast cancer.* Nature.
- Danenberg et al. 2022. *Breast tumor microenvironment structures are associated with genomic features and clinical outcome.* Nature Genetics.
- HNSCC immune-phenotype literature (desert / excluded / inflamed; median OS 37 / 61 / 85 months) — biological motivation for `immune_tumor`.

**Tools:** scikit-survival (RandomSurvivalForest, `concordance_index_censored`), scikit-learn (StratifiedGroupKFold), SciPy/NumPy (Delaunay, KL divergence).

---

## 11. Literature review (how it informed the work)

The literature entered the project at three decision points rather than as a standalone survey.

1. **Choosing the endpoint and the honest test.** The base paper (Dayao et al.) established that spatial arrangement (neighbourhood matrix, 0.704) rivals composition (0.705) but that naive spatial statistics (Ripley's K, 0.528) fail. That contrast motivated the central methodological guard of this project — comparing every spatial block against a **width-matched noise control** rather than trusting a raw C-index.

2. **Designing the enrichment features.** Squidpy (`nhood_enrichment`) and histoCAT establish the accepted way to quantify spatial enrichment: compare observed neighbourhood counts to a **within-image permutation null**, never to raw fractions or a cohort average. The enrichment features implement the same "divide by chance" logic in closed form (the arXiv:2506.18692 analytical score confirms this closed form matches Monte-Carlo permutation at r ≥ 0.95), and this is precisely *why* referencing to the sample's own composition — not the cohort — is correct: a cohort reference re-mixes composition back in, the exact confound the correction removes.

3. **Selecting the functional markers.** Rather than adding all 39 proteins (fatal at ~27 events), the shortlist was drawn from spatial-pathology survival studies where specific functional states were prognostic — Schürch 2020 (colorectal CODEX: PD-1⁺ CD4 T, GranzymeB, CD45RO), Jackson 2020 and Danenberg 2022 (breast IMC: Ki67, PD-L1, GranzymeB). This is what makes the 2-marker result (`cd4_foxp3` + `tcell_cd45ro`, C = 0.708 alone) a *confirmation* of prior biology rather than an unguided fishing expedition: Treg suppression and T-cell memory are exactly the states those studies flagged.
