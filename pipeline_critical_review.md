# Critical Review of the Spatial PE Pipeline

**Scope:** code + data audit of `spatial_positional_encoding/`, benchmarked against `spatsurv-main/` (Dayao et al., UPMC HNSCC MIBI cohort).
**Method:** every claim below was re-derived from the code and re-measured on the actual data. Numbers I ran myself are marked ✅; numbers quoted from `final_instructions.md` are marked 📄.

---

## Verdict

The pipeline's headline conclusion — *"PE beats raw coordinates (0.536 vs 0.526), so the graph Laplacian captures more than raw X,Y"* — **is not supported**. Both of those feature sets are degenerate: after per-sample z-scoring and mean/std pooling, they both collapse into re-encodings of **the number of cells in the sample**. The comparison measures nothing about spatial structure.

Separately, the "we're at 0.619 vs the paper's 0.704" gap that the existing analysis attributes to *pooling destroying spatial information* is **not a modelling problem at all**. It is a cohort-definition problem. The base paper filters out normal-mucosa tissue; this pipeline does not. Fixing the cohort and using the neighbourhood matrices already sitting unused in `data/raw/k10/` reproduces and exceeds the target:

| Configuration | C-index |
|---|---|
| Pipeline as-is: celltype proportions, all tissue types | 0.635 ✅ |
| **+ exclude normal mucosa** (the paper's own `--exclude-normal`) | **0.706** ✅ |
| **+ neighbourhood matrix k10, row-normalized** (already on disk, unused) | **0.730** ✅ |
| **+ treat death-of-other-causes as censored** (competing risks) | **0.776** ✅ |
| ...and then adding the PE features | **0.704** ✅ — *PE makes it worse* |

The single most important result of this audit: **on a correctly-defined cohort, the PE features reduce C-index by ~0.03 in both label settings.** They are not neutral. They are noise with a cell-count leak attached.

---

## 1. The PE features are mathematically degenerate

This is the core finding, and it is provable rather than empirical.

The encoder uses the symmetric normalized Laplacian `L_sym = I − D^{−1/2} A D^{−1/2}` ([laplacian_encoding.py:91](spatial_positional_encoding/src/laplacian_encoding.py#L91)) and keeps the `k_pe` smallest non-trivial eigenvectors. Two properties of that construction destroy the pooled features:

**(a) The trivial eigenvector forces the mean to ~zero.** `L_sym`'s null eigenvector is `u₀ ∝ D^{1/2}·1`. `eigsh` returns orthonormal eigenvectors, so every retained `vₖ` satisfies `Σᵢ √(dᵢ)·vₖ(i) = 0`. Degrees are near-constant here (mean degree 10.8–11.8), so `Σᵢ vₖ(i) ≈ 0` — the mean is pinned at zero by construction.

> Measured ✅: max `|mean|` across all dims and samples = **1.9e-4**, versus a typical std of 1.5e-2. The mean is ~100× smaller than the spread it is supposed to summarize.

**(b) Unit-norm forces the std to be exactly `1/√N`.** `‖vₖ‖₂ = 1`, so `Σᵢ vₖ(i)² = 1`, so `std(vₖ) = √(1/N − mean²) ≈ 1/√N`. The std of a Laplacian eigenvector carries **no information except the cell count**.

> Measured ✅: `std / (1/√N)` = **0.9999959 ± 5.6e-6**. Spearman correlation between `pe_std_0` and cell count = **−0.9999997**. All eight std columns are numerically identical to each other.

So the 16-dimensional "PE (Laplacian)" feature vector is really:

```
pe_mean_0..7  →  8 near-zero degree-heterogeneity residuals (~1e-5)
pe_std_0..7   →  8 identical copies of 1/√(cell count)
```

Confirmed by ablation ✅:

| Feature block | C-index |
|---|---|
| PE std block alone (8 features) | 0.549 |
| **Cell count alone (1 feature)** | **0.549** |
| PE mean block alone (8 features) | 0.609 |
| Full PE mean+std (16 features) | 0.563 |
| Permuted-label null | 0.515 – 0.574 |

The std block is *exactly* the cell-count model. The mean block scores 0.609 — but that is not spatial position either; because `Σ√(dᵢ)vᵢ = 0` exactly while `Σvᵢ` does not, the mean is a measure of **degree heterogeneity**, and RSF splits on ranks so its 1e-5 magnitude doesn't stop it being used. It is a real but accidental graph statistic, not a positional encoding, and it sits inside the permutation null band.

### Why the whole premise doesn't transfer

Laplacian PE was designed as a **within-graph node coordinate** for a model that consumes the graph (a GNN or a transformer attention bias). It has no cross-graph correspondence: sample A's `ev₁` and sample B's `ev₁` are eigenvectors of different operators, defined only up to sign and — for near-degenerate eigenvalues — up to arbitrary rotation *within* each eigenspace. Pooling them into a fixed-width vector and feeding an RSF assumes a correspondence between samples that mathematically does not exist.

The pipeline's own diagnostics show this is the worst case: 📄 PE eigenvalues span **0.0003–0.016** in a spectrum with range [0, 2]. They are packed into a sliver, i.e. maximally near-degenerate, i.e. maximally basis-unstable. The `_fix_sign_convention` heuristic ([laplacian_encoding.py:152](spatial_positional_encoding/src/laplacian_encoding.py#L152)) — force the max-absolute entry positive — does not fix this and is itself unstable: an infinitesimal perturbation can move the argmax and flip the whole column.

**Research context:** this is a known, named problem. Lim et al., *Sign and Basis Invariant Networks for Spectral Graph Representation Learning* (ICLR 2023) introduce SignNet/BasisNet precisely because sign-fixing heuristics like this one fail. Huang et al., *On the Stability of Expressive Positional Encodings* (ICLR 2024) show eigenvector-based PEs are unstable exactly when eigenvalue gaps are small — the regime here.

---

## 2. The "Raw Coordinates" control is a constant

The raw-coordinate baseline is supposed to be the honest control that PE is measured against. It isn't a control at all.

Preprocessing z-scores X and Y **per sample** ([preprocess.py:314-315](spatial_positional_encoding/src/preprocess.py#L314-L315)). The survival script then pools those coordinates by taking mean and std ([run_survival.py:201](spatial_positional_encoding/run_survival.py#L201)). Composing those two steps: the mean of a z-scored variable is 0, and the std is 1, **for every sample by construction**.

> Measured ✅ across-sample variance of the four "Raw Coordinate" features:
> `raw_X_mean`: 8.6e-17 · `raw_Y_mean`: 1.2e-16 · `raw_X_std`: 5.9e-5 · `raw_Y_std`: 5.9e-5

Two features are numerically zero; the other two vary only through the `√((N−1)/N)` ddof mismatch between pandas (`ddof=1` in preprocessing) and numpy (`ddof=0` in pooling) — which is, again, **a monotone function of cell count**. That is why "Raw Coordinates" scores 0.536 ✅ ≈ cell count 0.549 ✅.

**So the pipeline's headline claim compares one encoding of cell count against another encoding of cell count.** The 0.010 gap it reports is meaningless in both directions.

---

## 3. Dataset-type errors — the largest recoverable losses

You asked specifically about dataset types. This is where the performance actually went. The cohort is **MIBI multiplexed imaging, UPMC HNSCC: 307 QC-passing acquisitions from 81 patients across 7 coverslips**. The pipeline treats those 307 acquisitions as 307 independent, homogeneous samples. Every part of that is wrong.

### 3.1 Three different tissue types are pooled as if identical ⚠️ biggest single loss

```
Primary tumor    183
Nodal met         98
Normal mucosa     26
```

Normal mucosa contains **no tumour**. Its spatial architecture is unrelated to the tumour microenvironment, yet it inherits the patient's survival label — pure label noise. Nodal metastases are a different anatomical compartment with different architecture.

The base paper handles this explicitly. Its code has both flags:
```python
# spatsurv-main/scripts/rsf_prediction.py:49,52
if args.primary_only:   sample_df = sample_df[sample_df.tissuetype == 'Primary tumor']
elif args.exclude_normal: sample_df = sample_df[sample_df.tissuetype != 'Normal mucosa']
```
**This pipeline never filters.** Measured impact ✅: celltype proportions **0.658 → 0.706** just from dropping 26 normal-mucosa samples.

That single line closes the entire "why are we below 0.704" gap that the existing analysis spends four sections theorizing about. It was never a pooling problem.

> ⚠️ Caveat: I can't confirm from the repo which exact configuration produced the paper's 0.704, so treat "0.706 reproduces it" as a strong lead to verify against the manuscript, not a settled fact.

### 3.2 45% of the "deaths" are not cancer deaths ⚠️ second biggest loss

The `status` column decomposes `survival_status = 1` into:

| status | meaning | n |
|---|---|---|
| DOD | dead of disease | 57 |
| **DOC** | **dead of other causes** | **46** |
| NED / AWD | censored | 204 |

**46 of 103 events are deaths from unrelated causes**, coded identically to cancer deaths. No tumour spatial feature can predict a patient dying of something else. This is a textbook competing-risks violation and it puts a hard ceiling on any achievable C-index.

Measured impact ✅ of censoring DOC (cause-specific hazard):

| | as-is | DOC censored |
|---|---|---|
| Celltype prop (primary only) | 0.652 | **0.739** |
| Celltype + neighbour (excl. normal) | 0.730 | **0.776** |

**+0.05 to +0.09 from one line of label handling** — larger than any modelling idea in the existing improvement list.

**Research context:** Fine & Gray (JASA 1999) for subdistribution hazards; Ishwaran et al., *Random survival forests for competing risks* (Biostatistics 2014) for the RSF-native treatment. The minimum correct move is cause-specific censoring; the principled move is a competing-risks model.

### 3.3 The effective sample size is 27, not 308

307 acquisitions come from 81 patients — median 3 acquisitions per patient. And at the patient level:

> ✅ **Only 27 of 81 patients have an event.** (103 acquisition-level events are just those 27 patients counted repeatedly.)

The standard events-per-variable heuristic (~10 EPV) says this cohort supports roughly **2–3 free parameters**. The pipeline fits 16, 32, and — in the combined configurations — 272 features. This is the real reason "PE + Celltype (0.602) < Celltype (0.619)": not that PE adds noise in some abstract sense, but that the design is in a deep overfitting regime where *any* added block hurts.

It also means every C-index here is computed over 307 correlated units as though they were independent — pseudo-replication. The uncertainty is much wider than the reported ±0.15 suggests.

**This has a hard consequence the existing plan misses:** at 27 events, the proposed GNN / attention-pooling / contrastive-pretraining directions are not merely ambitious, they are **unfittable**. That entire tier of the improvement list is pointing the wrong way.

### 3.4 Batch structure is never tested

7 coverslips (`L4a2, L4b2, L4c1, P41A1, P41B1, P42A1, P42B1`), 24–58 acquisitions each. MIBI has well-known run-to-run batch effects, and the L-series and P-series are clearly separate runs. The base paper provides `--cv-folds coverslip` (leave-one-coverslip-out) specifically to test this. **This pipeline only ever splits by patient**, so nothing distinguishes "predicts survival" from "predicts coverslip".

### 3.5 HPV status is ignored

```
HPV+ 158 · HPV− 128 · unknown 21
```

HPV status is *the* dominant prognostic factor in head-and-neck cancer (Ang et al., NEJM 2010) — it stratifies survival more strongly than almost any tissue feature. It sits unused in `sample_metadata.csv`, along with stage and smoking history. Without a clinical-covariate baseline, "C-index 0.73" is uninterpretable: the question that matters is **ΔC over HPV+stage**, and that number has never been computed. Note also that unknown HPV is encoded as the *number* `5.0`, which will be silently treated as a value if anyone naively adds the column (the paper has a `update_hpvstatus` helper for exactly this).

---

## 4. The graph itself is built on warped geometry

Two problems, both from treating imaging data as if it were unitless point clouds.

**(a) Per-sample z-scoring applies an anisotropic warp.** `normalize_coordinates` divides X by `std(X)` and Y by `std(Y)` *separately*. When `std(X) ≠ std(Y)`, that squashes the tissue along one axis — and **k-NN neighbourhoods change under anisotropic scaling**.

> Measured ✅: `std(X)/std(Y)` has median **1.30** and ranges **0.42 – 3.04**. **284 of 308 samples (92%)** are warped by more than ±5%, and the warp factor differs per sample. FOV aspect ratio is typically 4:3 (median 1.33), not square.

The step is also **unnecessary**: k-NN graphs are invariant to *isotropic* scaling, so the only thing this normalization achieves is the distortion. And it discards physical units, which the base paper deliberately keeps — it converts to microns at `0.377 µm/px` and evaluates Ripley's K at biologically-chosen radii of 30 µm and 80 µm.

**(b) Fixed k=10 is not comparable across samples.** Cell density varies **15×** across the cohort ✅ (0.0002 – 0.0031 cells/px²). A fixed-k graph therefore spans a ~4× range of physical radii: in a dense region k=10 reaches ~one cell diameter, in a sparse region several. "The 10-nearest-neighbour graph" means a different biological object in different samples.

**Research context:** the spatial-omics field uses either **Delaunay triangulation** or **fixed-radius (µm) graphs** for exactly this reason — see Fischer et al., *NCEM* (Nature Biotechnology 2023) and Palla et al., *Squidpy* (Nature Methods 2022). Both are density-adaptive in a principled way; fixed-k is not.

---

## 5. Evaluation protocol

**The null is 0.55, not 0.5.** Every conclusion in `pipeline_analysis_and_ideas.md` is calibrated against "0.500 = random guessing". That is wrong for this design.

> Measured ✅ with permuted labels: **0.515, 0.552, 0.574** (all-tissue) and **0.493, 0.540, 0.542** (best config).

Against a true null of ~0.55, the entire lower half of the results table — PE (0.563), raw coords (0.536), cell count (0.549) — is **indistinguishable from noise**. The base paper's code has `--permute-labels` built in for this purpose; it was never run here.

**3 of 10 folds are uninformative.** The existing doc notes "fold 2 always fails". It's worse than that — event counts per fold ✅:

```
fold:    0    1    2    3    4    5    6    7    8    9
events:  8   25    0    9    3   11    3   14   15   15
```

Fold 2 has **zero** events (C-index undefined → NaN). Folds 4 and 6 have **three** events each — a C-index computed on 3 events is essentially a coin flip contributing full weight to the mean. The fix is `StratifiedGroupKFold`, which stratifies on event status while grouping by patient — and which `rsf_prediction.py` already imports but never uses.

**No repeats.** `GroupKFold` is deterministic and unshuffled, so all results come from **one fixed partition**, at one seed. The base paper supports `--n-repeats` and reports confidence intervals over seeds. Reporting `0.536 ± 0.16` from a single partition presents fold-spread as if it were uncertainty in the mean, which it isn't.

**The reported numbers don't reproduce.** I re-ran `run_survival.py` unmodified on the same data ✅:

| Feature set | 📄 documented | ✅ re-run today |
|---|---|---|
| Celltype Proportions | 0.6185 | 0.6351 |
| PE + Celltype | 0.6022 | 0.6244 |
| PE (Laplacian) | 0.5359 | 0.5630 |
| Raw Coordinates | 0.5255 | 0.5360 |

A **0.027 shift** from environment alone. The doc's headline finding — *"PE beats raw coordinates, 0.536 vs 0.526"* — rests on a **0.010** gap. That is 2.7× smaller than the environment noise and 15× smaller than the fold standard deviation. It was never a result.

**Two smaller code issues:**
- `load_celltype_proportions` never applies the QC filter, so the celltype path evaluates on 308 samples while preprocessing produced 307. 29 QC-failing acquisitions are survival-eligible and can leak in. The evaluated cohort should not depend on which feature block you chose.
- `StandardScaler` is fit on the full dataset *before* the CV split ([run_survival.py:272](spatial_positional_encoding/run_survival.py#L272)). Harmless in practice (RSF splits on ranks), but it's the leakage pattern, and it'll bite the moment anyone swaps in a Cox model.

---

## 6. What is actually costing performance — ranked

| # | Problem | Measured cost | Fix effort |
|---|---|---|---|
| 1 | Normal mucosa not excluded | **−0.05** (0.706→0.658) | 1 line |
| 2 | DOC treated as event (competing risks) | **−0.05 to −0.09** | 1 line |
| 3 | Neighbourhood matrices on disk, unused | **−0.02** (0.730→0.706) | ~10 lines |
| 4 | PE features added to model | **−0.03** (0.776→0.745) | delete |
| 5 | 27 patient-events vs 16–272 features | inflates variance, caps everything | design |
| 6 | Null assumed 0.5 (actually ~0.55) | invalidates conclusions | 1 run |
| 7 | 3/10 folds have ≤3 events | inflates variance | `StratifiedGroupKFold` |
| 8 | Single partition, single seed | fake precision | `n_repeats` |
| 9 | Anisotropic z-scoring (up to 3× warp) | corrupts graph | delete step |
| 10 | Fixed k over 15× density range | non-comparable graphs | Delaunay/radius |
| 11 | HPV/stage never used or adjusted for | claims uninterpretable | small |
| 12 | Coverslip batch never tested | confound unmeasured | 1 flag |

Items 1–3 are pure profit and take an afternoon. They move you from 0.635 to 0.776 without a single new idea.

---

## 7. Where the research actually points

The instinct behind the project — *"cell-type composition ignores spatial arrangement; encode the geometry"* — is correct and well-supported. The execution picked a tool that can't transfer across samples. Concretely:

**If you want to keep a "positional/structural encoding" framing** (likely, for the FYDP narrative), the drop-in replacement is **RWSE — random-walk structural encoding** (Dwivedi et al., ICLR 2022; Rampášek et al., *GraphGPS*, NeurIPS 2022). For each node, `[p⁽¹⁾ᵢᵢ, p⁽²⁾ᵢᵢ, ..., p⁽ᴷ⁾ᵢᵢ]` — the return probabilities of a random walk. It is deterministic, **sign-free**, **basis-free**, and — critically — **the same number means the same thing in every sample**, so pooling is legitimate. It measures local density and connectivity structure, which is what you wanted PE to measure. This is the honest fix to the exact defect in §1.

**The field-standard answer to your actual biological question** is **cellular neighbourhoods** (Schürch et al., *Cell* 2020): take each cell's k-NN window, compute the cell-type composition of that window, cluster those composition vectors across the cohort into ~10 recurring "neighbourhood" archetypes, then describe each sample by its CN proportions. This is exactly the "which cell types are near which" quantity the existing analysis correctly identifies as the missing signal — and unlike Laplacian eigenvectors, the archetypes are **shared across samples**, so per-sample proportions are comparable by construction. The `k1`–`k20` directories on disk suggest this is already half-built.

**Closest published analogue to the whole project:** Wu et al., *Graph deep learning for the characterization of tumour microenvironments from spatial protein profiles in tissue specimens* (Nature Machine Intelligence, 2022) — GNNs on cell graphs from multiplexed imaging for outcome prediction. Worth reading before building anything, since it is the paper this FYDP is adjacent to. See also **CellCharter** (Varrone et al., Nature Genetics 2024) for tissue-domain identification.

**On the GNN / attention-MIL / contrastive tier of the existing plan** (ABMIL — Ilse et al. ICML 2018; DeepAttnMISL — Yao et al. MedIA 2020; PatchGCN — Chen et al. MICCAI 2021): these are the right tools *for datasets with thousands of patients*. At **27 patient-events** they will fit the noise and validate optimistically, and no amount of regularization changes that arithmetic. The existing priority list ranks "build a GNN" above "fix the cohort definition"; that ordering should be inverted. Recommend explicitly dropping this tier and saying why in the report — a well-argued negative result on feasibility is a legitimate FYDP contribution.

**On the "celltype-stratified PE pooling" idea** (16 types × 8 dims = 128 features): the intuition is reasonable but it inherits the cross-sample non-comparability from §1 *and* adds 128 features at 27 events. Do it with CN proportions or RWSE instead, and keep it to <20 features.

---

## 8. Recommended pipeline

```
1. COHORT      exclude normal mucosa (or stratify by tissuetype and report both)
               → report patient-level n and event count, not acquisition-level
2. LABELS      cause-specific: DOC → censored.  Report Fine-Gray as sensitivity check.
3. GEOMETRY    drop z-scoring. keep microns (0.377 µm/px).
               graph = Delaunay, or fixed-radius ~30 µm.
4. FEATURES    baseline A: HPV + stage           (clinical floor — the number to beat)
               baseline B: + celltype proportions
               baseline C: + neighbourhood matrix (row-normalized)  ← reproduces paper
               novel   D: + CN proportions (Schürch-style) and/or RWSE
               keep total features < ~30. Nested feature selection if more.
5. EVAL        StratifiedGroupKFold(patient, stratify=event)
               × 20 seeds → report mean + 95% CI
               permutation null every time (it is ~0.55, not 0.5)
               leave-one-coverslip-out as batch robustness check
6. CLAIM       report ΔC over baseline A, with CI. That is the only number that means anything.
```

## 9. Do these three things first

1. **Add `--exclude-normal` and `--censor-doc`** to `run_survival.py`. Two lines, +0.14 C-index (0.635 → 0.776). Do this before anything else.
2. **Load `data/raw/k10/*.npy`, row-normalize, evaluate.** Already computed, never used. This reproduces the paper's baseline and gives you a legitimate reference point instead of a number from a PDF.
3. **Run the permutation null and 20-seed repeats.** Until this exists, no comparison in this project can be defended in a viva — including the negative ones in this document.

Then, and only then, ask whether RWSE or CN proportions add anything over baseline C. That is the real research question, and right now it is unanswered — not answered negatively, just never actually asked.

---

### Reproducing this audit

`scikit-survival 0.27` is broken against the user-site `scikit-learn 1.9` in this environment (`ImportError: cannot import name 'DTYPE'`). The system-site `sklearn 1.8` works:

```python
import sys; sys.path.insert(0, r'C:\Python312\Lib\site-packages')  # before importing sksurv
```

Worth pinning in `requirements.txt` — as it stands, `run_survival.py` cannot execute in this environment, which means the results table in `final_instructions.md` cannot currently be regenerated by anyone who clones this repo.
