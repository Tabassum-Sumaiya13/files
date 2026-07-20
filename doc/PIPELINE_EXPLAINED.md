# The Pipeline, End to End — How It Works, What It Finds, What the Results Mean

## 0. What the project is actually trying to prove

The data is spatial proteomics: for each tissue sample, a list of **cells**, each
with an `(X, Y)` position, a **cell-type label**, and ~40 **protein-marker
intensities**. The scientific claim is:

> *How cells are spatially arranged carries biological information — beyond just
> how many of each cell type are present.*

The hard part is proving that claim isn't an illusion. A feature can look
impressive and actually be (a) noise, (b) a disguised copy of "how much tumour is
in this sample," or (c) a fluke of one dataset. The whole pipeline exists to
**build spatial features and then adversarially test whether they're real**.
Survival was the *original* test on the UPMC cohort; it has been generalized so
any dataset — even one with no survival, like CRC — can be tested.

---

## 1. Two phases, one path for every dataset

```
   RAW COHORT (messy, dataset-specific)
        │
        │  PHASE 1 — INGEST        data_preprocessing/run_ingest.py
        │  "make every dataset look identical"
        ▼
   datasets/<NAME>/processed/         ← canonical format, same for UPMC, CRC, anything
        ├─ samples/sample_*.parquet   (X, Y, cluster_label, lineage, markers, [survival])
        ├─ manifest.parquet           (patient_id, survival if any)
        └─ marker_columns.txt
        │
        │  PHASE 2 — VERIFY        spatial_positional_encoding/run_verify.py
        │  "are the spatial features real?"
        ▼
   verification/verify_<taxonomy>.csv   ← the results matrix
                perturbation_<taxonomy>.csv
                separability_<taxonomy>_<label>.csv
```

The point of Phase 1 is that Phase 2 never has to know which dataset it is looking
at. UPMC and CRC arrive at Phase 2 in the exact same shape.

Positional encoding (Laplacian PE) has been **removed** from the pipeline; the
spatial feature block is the abundance-corrected **enrichment** readout.

**Currently ingested** (`processed/manifest.parquet` exists):

| dataset | samples | survival | what it is |
|---|---|---|---|
| UPMC | 308 | yes | the project's baseline head-and-neck CODEX cohort |
| CRC | 140 | no | external validation cohort (Schürch et al. 2020) |
| CRC_doublets | 140 | no | **measurement variant of CRC**, not a third cohort — see §5.3 |

`keren_tnbc` and `hubmap_intestine_codex` have adapter configs but no
`processed/` directory: they were never successfully ingested.

---

## 2. Phase 1 — Ingest: raw → canonical

Entry point: `data_preprocessing/run_ingest.py --dataset <NAME>`.

### 2.1 The adapter is the only dataset-specific thing

`run_ingest.py` dynamically imports `datasets/<NAME>/adapter_config.py`
(`load_config`, run_ingest.py:29-40). Everything downstream is generic. The
adapter declares only:

| declaration | CRC example |
|---|---|
| `LOCATIONS_/EXPRESSION_/METADATA_PATH` | all three point at the one combined CSV |
| `*_COLUMN_MAP` | `"File Name"→acquisition_id`, `"X:X"→X`, `"ClusterName"→cluster_label` |
| `marker_columns(df)` | "any column carrying a `:Cyc_` channel tag, minus nuclear stains" → 56 |
| `CELLTYPE_MAP` | **not hardcoded** — `registry.celltype_map("CRC")` reads the registry CSV |
| `APPLY_ARCSINH`, `ARCSINH_COFACTOR` | `True`, `5.0` (UPMC: `False`, already normalised) |
| `MIN_CELLS_PER_SAMPLE` | `50` |

### 2.2 Validate — read-only, and it can stop the run

`validator.validate_dataset(cfg)` runs ~15 checks and writes
`validation_report.md`. **Nothing under `processed/` is touched.** In order:

| # | check | code | failure mode |
|---|---|---|---|
| 1 | raw files exist | validator.py:50-64 | FAIL — stops immediately |
| 2 | required canonical columns present after rename | validator.py:71-87 | FAIL — stops immediately |
| 3 | X/Y numeric; survival numeric + binary **if present** | validator.py:96-125 | missing survival = **WARN, not FAIL** |
| 4 | marker columns detected | validator.py:127-132 | FAIL |
| 5 | acquisition_id overlap across the 3 tables | validator.py:134-149 | FAIL if empty; orphans WARN |
| 6 | >1 patient (patient-grouped CV would otherwise leak) | validator.py:151-162 | FAIL |
| 7 | sample-size distribution vs `MIN_CELLS_PER_SAMPLE` | validator.py:164-173 | WARN |
| 8 | cell-type mapping coverage | validator.py:175-200 | WARN — unmapped types get dropped |
| 9 | registry structural checks | `registry.validate_registry()` | see below |
| 10 | marker-panel resolution | `markers.resolve_panel()` | WARN |
| 11 | **falsifiability gate** — registry vs marker evidence | `lineage_evidence.evaluate()` | FAIL if material |

Check 3 is what lets a survival-less cohort through: survival absent is a WARN,
so CRC ingests and is fully feature-verifiable.

Check 9 (`registry.py`) tests things a hardcoded dict could never be tested for:
no duplicate `native_label`; every lineage in the canonical vocabulary; **CL
grounding** — a row's declared lineage must equal
`schema.CL_LINEAGE_ANCHOR[cl_term_id]`, so two rows citing the same ontology term
can never disagree; every mapped row carries a citation; every excluded row
records *why*; and the registry covers exactly the labels present in the raw data
(an unregistered label is a FAIL, not a silent drop).

Check 10 (`markers.py`) exists because the old validator compared marker names
with exact string equality and reported "0/10 lineage markers" for a cohort whose
panel plainly contains `'CD45 - hematopoietic cells:Cyc_4_ch_2'`. Resolution now
strips the CODEX channel tag and free-text description, casefolds, and then
consults a curated alias table (`Cytokeratin ≡ PanCK`, `CD3 ≡ CD3e`). Matching is
on the whole normalised token, so `CD4` can never match `CD45`.

### 2.3 The falsifiability gate (check 11)

`lineage_evidence.py` confronts the registry with the cohort's own expression
data, so a mapping cannot enter the pipeline unchallenged:

```
per-cluster marker means
  → z-score each marker ACROSS clusters
  → lineage score = MAX z over that lineage's CORE markers
        (max, not mean: immune is heterogeneous — a T cell is CD45+/CD20−/CD68−,
         so averaging a lineage's subset anchors penalises every subset for not
         being the others)
  → re-standardise the three lineage scores across clusters
        (so a 2-marker panel and a 5-marker panel are commensurable and the
         argmax is not a panel-size artifact)
  → predicted = argmax;  margin = top − runner-up
  → verdict: AGREE / AMBIGUOUS (margin < 0.25) / CONTRADICTED / EXCLUDED
```

It **refuses to run** if any of the three lineages has no resolvable core marker —
an argmax over the remaining lineages would be rigged against the unscoreable one.

`summarise()` grades the outcome by materiality rather than by mere existence: a
contradiction covering ≥5% of mapped cells **FAILs** the ingest, smaller ones
WARN — unless the registry row carries an `evidence_override` with a written
justification, which downgrades it from blocking to reported. An overridden
contradiction is still counted, still listed, and still carried into
`--perturb-map`. Acceptance is a reason, not a result.

The per-cluster table is written to `datasets/<NAME>/lineage_evidence.csv`
(run_ingest.py:62-65) so Phase 2's flip scenarios can read the predictions
instead of parsing markdown.

**Gate:** processing runs only if the report is READY (zero FAILs) or `--force`
is passed (run_ingest.py:70-74).

### 2.4 Process — the exact transformation, logged step by step

`processor.process_dataset()` performs this sequence, writing before/after counts
for each step into `processing_report.md`:

| # | step | code | note |
|---|---|---|---|
| 1 | load 3 tables, rename to canonical | processor.py:44-49 | |
| 2 | dedup on `(acquisition_id, cell_id)` | processor.py:58-64 | |
| 3 | drop null X/Y, drop null marker values | processor.py:67-73 | |
| 4 | `arcsinh(x / cofactor)` on markers | processor.py:79-84 | skipped when `APPLY_ARCSINH=False` (UPMC) |
| 5 | drop marker columns from the *locations* view | processor.py:86-101 | only fires when a cohort ships locations+expression as one file (CRC) — see §11 |
| 6 | inner-merge locations + expression | processor.py:103-108 | unmatched location rows dropped |
| 7 | drop samples below `MIN_CELLS_PER_SAMPLE` | processor.py:110-116 | |
| 8 | `cluster_label → lineage` | processor.py:119-133 | **unmapped cells are dropped**, with a WARN listing the types |
| 9 | per-sample z-score of X and Y | processor.py:136-144 | |
| 10 | merge patient_id + survival | processor.py:146-157 | survival optional → `NaN` |
| 11 | write one parquet/sample + manifest + marker list | processor.py:160-189 | |

Result on the two production cohorts:

| | UPMC | CRC |
|---|---|---|
| samples | 308 | 140 |
| cells | 2,061,102 | 240,554 |
| patients | 81 | 35 |
| native cell types | 16 | 29 (25 kept) |
| survival | yes (103 events) | **none** |

### 2.5 New cell types

New cell types are handled safely at every layer, but not kept automatically:

- **Ingest** — a native type not in the registry is **dropped** (validator FAILs
  on `registry:coverage` if it has no row at all; WARNs and lists it if the row
  exists with a blank lineage). It never crashes. To KEEP a new type, give it a
  lineage in `celltype_registry.csv`.
- **Verify, native taxonomy** — fully dynamic: the category vocabulary is built
  from whatever labels exist, so a new type is auto-discovered.
- **Verify, lineage taxonomy** — new types were already collapsed to a lineage at
  ingest, so the 3-category space stays complete.
- **Counting internals** — out-of-vocabulary cells are dropped, not mis-indexed.

---

## 3. What the features actually measure

The pipeline builds a **spatial graph** (Delaunay triangulation — connect each cell
to its natural tissue neighbours; parameter-free, ~6 neighbours each) and reads out
**5 "enrichment" numbers per sample**. Each is designed so that **0 = a random
arrangement** — they only become non-zero if cells are organized:

| feature | plain-English meaning | expected sign |
|---|---|---|
| `kl_mean` | How different is a typical cell's neighbourhood from the tissue's overall mix? | `>0` = structured |
| `kl_tumor` | How distinctive is the tumour's neighbourhood specifically? | `>0` = tumour sits in its own niche |
| `self_enrich` | Do cells of the same type cluster together more than chance? | `>0` = they clump |
| `immune_tumor` | Do immune cells touch tumour cells more/less than chance? | sign = infiltration vs exclusion |
| `stroma_tumor` | Do stromal cells touch tumour more/less than chance? | sign = interface structure |

These come in two resolutions:

- **lineage** (3 categories) — portable, comparable across datasets.
- **native** (`cluster_label`, 16 for UPMC / 29 for CRC) — finer detail.

**Both resolutions depend on the cell-type → lineage mapping**, and it is worth
being precise about how much. Under `--taxonomy native` the category *vocabulary*
is built from whatever labels exist, but `build_vocab` still groups those
categories into immune / tumour / stromal buckets by the majority `lineage` of
their cells — and `lineage` was assigned from the registry at ingest. So:

| feature | needs the mapping? |
|---|---|
| `kl_mean` | no — only needs *a* partition of the cells |
| `self_enrich` | no |
| `kl_tumor` | **yes**, at both resolutions |
| `immune_tumor` | **yes**, at both resolutions |
| `stroma_tumor` | **yes**, at both resolutions |

The mapping is therefore not a detail of ingest: three of the five headline
features are defined in terms of it. That is why it is a cited, versioned,
ontology-grounded registry that marker evidence can contradict
(`data_preprocessing/celltype_registry.csv`, see doc/CELLTYPE_MAPPING.md), and
why `run_verify.py --perturb-map` re-runs the whole matrix with every contested
assignment dropped and flipped.

The **baseline** every feature is tested against is *composition* — the fraction of
each cell type in the sample. That is the "boring" explanation: "you only need to
know how much tumour there is, not where it is." A feature only matters if it beats
that baseline.

---

## 4. Phase 2 — Verify: the three tests, and what each *finds*

Entry point: `spatial_positional_encoding/run_verify.py --dataset <NAME> --taxonomy <lineage|native>`.

### 4.1 Code flow

1. **Load** — `cohort.load_cohort()` reads `manifest.parquet` + `marker_columns.txt`;
   `iter_samples()` streams the per-sample parquets, reading only the columns needed
   (`X, Y, lineage, cluster_label`).
2. **Build the category vocabulary** — `build_vocab()` (run_verify.py:61). For
   `--taxonomy lineage` the categories are exactly `(immune, tumour, stromal)`.
   For `--taxonomy native` they are every `cluster_label` in the cohort, each
   bucketed into a lineage by the **majority lineage of its cells** — which is how
   native mode still inherits the registry mapping.
3. **Per sample** (`compute_cohort()`, run_verify.py:121):
   - build the Delaunay graph from X/Y (`spatial_features.delaunay_edges`); edges
     touching an out-of-vocabulary cell are dropped
   - `props` = composition proportions ← **the baseline**
   - `real` = the 5 enrichment scalars from the K×K neighbour-count matrix
     (`enrichment_scalars`)
   - **permutation null**: hold the graph fixed, `rng.permutation(idx)` × `n_perm`,
     recompute each time (run_verify.py:150-154)
   - **split-half**: two disjoint random halves, each gets its *own* Delaunay
     graph, recompute both (run_verify.py:157-166)
4. **Assemble the matrix** — `build_matrix()` (run_verify.py:224) collapses the
   per-sample arrays into one row per feature.

### 4.2 The three tests

**① `null_z_median` — "Is it real spatial signal, or chance?"**
Keep the graph fixed, but randomly **shuffle the cell labels** (default 20 times).
That destroys spatial organization while keeping composition identical. Recompute
the feature on each shuffle → that is the "chance" distribution.
`z = (real − chance_mean) / chance_std`.
→ **Finds:** whether the real value stands apart from what random placement gives.
`|z| ≥ 2` means real. `frac_|z|>2` is the fraction of samples clearing that bar.

**② `spatial_specific` — "Does it add anything beyond composition?"**
Try to predict the feature from the composition baseline alone (5-fold regression).
`spatial_specific = 1 − R²` (R² clamped to `[0, 1]`).
→ **Finds:** the fraction of the feature that composition *cannot* explain. `1.0` =
purely spatial; `≥ 0.5` = it genuinely adds. This is the label-free version of
"beats the celltype baseline."

**③ `stability_r` — "Is it reproducible or a fluke?"**
Split each sample's cells into two random halves, rebuild the graph on each half,
recompute, and correlate the two halves across all samples.
→ **Finds:** test–retest reliability. `≥ 0.5` = stable; near `1.0` = you would
measure the same thing on independent halves of the tissue.

**Verdict** = **STRONG** only if all three pass: real (`|z| ≥ 2`) **and** adds
(`spec ≥ 0.5`) **and** stable (`r ≥ 0.5`).

### 4.3 The two optional add-ons

**`--perturb-map`** (run_verify.py:353-417) takes the contested clusters straight
from the registry — exactly the rows carrying an `evidence_override` — and re-runs
the entire matrix under each of:

| scenario | question it answers |
|---|---|
| `baseline` | the registry as it stands |
| `drop:<cluster>` | is the result carried by that cluster? |
| `flip:<cluster>→<L>` | reassign it to what the marker evidence predicted |
| `evidence-all` | flip every contested cluster at once — the worst case |

It then prints plainly whether **any verdict moved**. Nothing here is
dataset-specific: a new cohort gets its own scenarios with no code change, and a
cohort with no contested rows gets only the baseline.

**`--label-col <col>`** (`separability()`, run_verify.py:265) is the optional 4th
test, available only when a dataset has a categorical label in its manifest (e.g.
CRC's CLR/DII `groups`): GroupKFold-by-patient classification AUC comparing
`Composition baseline` vs `Baseline + spatial block` vs `Baseline + noise` vs
`Spatial block alone`. A real block beats both baseline and noise. This never uses
survival as a proof of the features — see §5.4 for what it actually returned.

---

## 5. Reading the results

### 5.1 The verification matrix

**UPMC — lineage** (`datasets/UPMC/processed/verification/verify_lineage.csv`)

| feature | null_z | spatial_specific | stability_r | meaning |
|---|---|---|---|---|
| kl_mean | **1398** | 0.77 | 0.98 | wildly non-random, mostly spatial, rock-solid |
| kl_tumor | **2141** | 0.92 | 0.99 | tumour niche is extremely distinctive |
| self_enrich | 31 | 0.55 | 0.98 | cells clump by type, reproducibly |
| immune_tumor | −52 | **1.00** | 0.99 | immune–tumour contact is *entirely* spatial (composition explains 0%) |
| stroma_tumor | −52 | 0.83 | 0.97 | stromal–tumour interface is real and specific |

**CRC — lineage, NO survival** (`datasets/CRC/processed/verification/verify_lineage.csv`)

| feature | null_z | spatial_specific | stability_r | meaning |
|---|---|---|---|---|
| kl_mean | 81 | 0.94 | 0.96 | real, specific, stable |
| kl_tumor | 93 | 0.99 | 0.96 | tumour niche distinctive |
| self_enrich | 7.6 | 1.00 | 0.62 | clumping real & purely spatial; noisier split-half |
| immune_tumor | −9.4 | 0.87 | 0.88 | immune exclusion, real |
| stroma_tumor | −4.6 | 0.85 | 0.80 | real, weakest of the five |

**UPMC — native, 16 types** (`verify_native.csv`)

Same story, and `self_enrich` jumps to **2.21** — at fine resolution, cells of a
*specific* subtype self-segregate even harder than at the coarse lineage level. A
sanity check that finer labels reveal finer structure, exactly as they should.

**How to read a row in one sentence:** *"This feature is `null_z`
standard-deviations away from random, composition can't explain `spatial_specific`
of it, and it reproduces at `stability_r` on independent halves of the tissue."*

### 5.2 Map sensitivity — nothing depends on the contested rows

`perturbation_lineage.csv`, both cohorts: **every feature is STRONG in every
scenario**, including `evidence-all`.

| cohort | contested clusters (registry `evidence_override`) |
|---|---|
| UPMC | `APC`, `Naive immune cell`, `Tumor (Podo+)` |
| CRC | `plasma cells`, `CD11b+ monocytes`, `CD4+ T cells GATA3+`, `CD163+ macrophages` |

Magnitudes move — flipping UPMC's `Tumor (Podo+)` to stromal lifts `kl_mean`'s
spatial-specific from 0.77 to 0.92, and CRC's `flip:plasma cells→tumour` pushes
`stroma_tumor`'s z from −5.5 to −8.4 — but no verdict changes. **The conclusions
do not rest on the contested cell-type assignments being right.**

### 5.3 The doublet-exclusion control (CRC_doublets)

Production CRC drops 4 native types (6.9% of cells), two of which are doublets
(`tumor cells / immune cells`, `immune cells / vasculature`). Dropping them is
defensible, but the bias is *not neutral*: those cells sit at the immune–tumour
interface, and three of the five features measure exactly that interface.

`CRC_doublets` reads the same raw file with the doublets mapped instead of
dropped. Its baseline (`kl_mean` z 92, `kl_tumor` 101, `self_enrich` 8.5,
`immune_tumor` −8.6, `stroma_tumor` −4.4) tracks production CRC closely and
**every verdict is identical**. The exclusion does not carry the result.

This is a measurement variant, **not independent evidence** — do not report it as
a third cohort.

### 5.4 Separability against survival status — a negative result

`separability_lineage_survival_status.csv` (UPMC):

| feature set | AUC | n_feats |
|---|---|---|
| Composition baseline | 0.603 | 3 |
| Baseline + spatial block | **0.492** | 8 |
| Baseline + noise | 0.628 | 8 |
| Spatial block alone | 0.460 | 5 |

The spatial block does **not** predict UPMC survival status — it does worse than
adding random noise to the same baseline. State this plainly: the features are
verified as real, reproducible spatial structure; they are **not** verified as
prognostic. Those are two different claims and this file separates them.

---

## 6. What the overall result means

- **Every feature is STRONG on both datasets.** The spatial features are (1) not
  noise, (2) not a disguised copy of cell-type counts, (3) reproducible.
- **The two magnitudes differ, and that is meaningful.** UPMC's z-scores are in the
  thousands, CRC's in the tens. Both are decisively real, but UPMC's tissue is more
  strongly organized / has more cells per sample. CRC's `stroma_tumor` (z −4.6,
  stability 0.80) is the feature to watch if signal ever gets marginal.
- **This conclusion required no survival.** CRC has none, yet the features are
  still shown to work with evidence. That is the entire objective.
- **The conclusion is robust to the mapping** (§5.2) and to the doublet exclusion
  (§5.3) — both demonstrated, not asserted.
- **It does not extend to outcome prediction** (§5.4). Real ≠ prognostic.

---

## 7. Where survival fits now

Survival is no longer part of proving the features work — it is a **separate,
optional, per-dataset downstream question**: *"do these validated features also
predict patient outcome?"* (`spatial_positional_encoding/run_survival.py`,
RandomSurvivalForest + patient-grouped C-index). A dataset with survival gets that
as a bonus; a dataset without it (CRC) is still fully verified by Phase 2. On
UPMC, the separability check in §5.4 says the answer to that downstream question
is currently **no**.

---

## 8. How to check it manually, step by step

Every command below runs against the committed outputs. Expected results are
stated so a mismatch is immediately visible.

### Step 0 — what is ingested
```bash
cd spatial_positional_encoding
python run_pipeline.py --list
```
Expect `CRC 140 / CRC_doublets 140 / UPMC 308`, with `survival=no/no/yes`.

### Step 1 — read the two audit reports first
`datasets/<NAME>/validation_report.md` and `processing_report.md`. The processing
report is a literal ledger of §2.4 with before/after row counts. Every number in
that table should appear there. If a step's counts look wrong, the bug is in that
step and nowhere else.

### Step 2 — coordinates really are z-scored per sample
```python
import pandas as pd
s = pd.read_parquet('data_preprocessing/datasets/CRC/processed/samples/sample_reg001_A.parquet')
print(s['X'].mean(), s['X'].std(), s['Y'].mean(), s['Y'].std())
```
Expect ≈ `-0.0, 1.0, -0.0, 1.0`.

### Step 3 — the lineage in the parquet matches the registry
```python
import pandas as pd
s = pd.read_parquet('data_preprocessing/datasets/CRC/processed/samples/sample_reg001_A.parquet',
                    columns=['cluster_label', 'lineage'])
print(s.groupby(['cluster_label', 'lineage']).size())
reg = pd.read_csv('data_preprocessing/celltype_registry.csv')
print(reg[reg.dataset == 'CRC'][['native_label', 'lineage', 'cl_term_id', 'verified']].to_string())
```
Every `(cluster_label, lineage)` pair must appear in the registry, and no cluster
may carry two different lineages.

### Step 4 — confirm exactly which cell types were dropped
```python
import glob, pandas as pd
raw = pd.read_csv('data_preprocessing/datasets/CRC/raw/CRC_clusters_neighborhoods_markers.csv',
                  usecols=['ClusterName'])
proc = set()
for f in glob.glob('data_preprocessing/datasets/CRC/processed/samples/*.parquet'):
    proc |= set(pd.read_parquet(f, columns=['cluster_label'])['cluster_label'])
print(sorted(set(raw.ClusterName) - proc))
```
Expect exactly `dirt`, `undefined`, `immune cells / vasculature`,
`tumor cells / immune cells` — and each must carry a written reason in the
registry's `notes` column.

### Step 5 — the Delaunay graph is sane
```python
import sys, pandas as pd; sys.path.insert(0, 'spatial_positional_encoding')
from src import spatial_features as sf
s = pd.read_parquet('data_preprocessing/datasets/CRC/processed/samples/sample_reg001_A.parquet',
                    columns=['X', 'Y'])
e = sf.delaunay_edges(s.values.astype(float))
print(len(s), len(e), 2 * len(e) / len(s))
```
Expect mean degree ≈ 6. Observed: 1093 cells, 3245 edges, degree **5.94**.

### Step 6 — recompute a feature and its null by hand (the important one)
```python
import sys, numpy as np, pandas as pd; sys.path.insert(0, 'spatial_positional_encoding')
from src import spatial_features as sf
s = pd.read_parquet('data_preprocessing/datasets/CRC/processed/samples/sample_reg001_A.parquet',
                    columns=['X', 'Y', 'lineage'])
idx = np.array([{'immune': 0, 'tumour': 1, 'stromal': 2}[l] for l in s['lineage']])
e = sf.delaunay_edges(s[['X', 'Y']].values.astype(float))
counts = np.bincount(idx, minlength=3).astype(float)
real = sf.enrichment_scalars(sf._count_matrix(e, idx, 3), counts, [0], [1], [2])
rng = np.random.RandomState(0)
null = np.array([sf.enrichment_scalars(sf._count_matrix(e, rng.permutation(idx), 3),
                                       counts, [0], [1], [2]) for _ in range(20)])
print('real', real.round(3))
print('null', null.mean(0).round(3))
print('z   ', ((real - null.mean(0)) / (null.std(0) + 1e-9)).round(1))
```
Observed: `real [0.019 0.027 0.491 -0.142 -0.151]`, `null [0.003 0.008 -0.055 -0.008 -0.003]`.

Two things to confirm here, and they are the crux of the whole argument:

- **the null means are ≈ 0** — the "0 = random arrangement" claim being true
  empirically, not by assertion
- **the shuffle touches only labels, never `e`** — same graph, same composition, so
  anything that moves is spatial arrangement and nothing else

### Step 7 — reproduce the published matrix
```bash
cd spatial_positional_encoding
python run_verify.py --dataset CRC --taxonomy lineage --n-perm 20 --seed 1029
```
Compare the printed table to `verify_lineage.csv`. Same seed → same numbers. Then
re-run with `--seed 7`: `null_z` will shift (it is a 20-draw null) but **every
verdict must stay STRONG**. A verdict that flips with the seed was never real.

### Step 8 — sanity-check the baseline test
Re-run with `--limit 20`. `spatial_specific` should move noticeably — a 3-column
regression behaves differently on 20 samples than on 140. If it is *identical* at
n=20 and n=140, be suspicious of the regression.

### Step 9 — mapping sensitivity
```bash
python run_verify.py --dataset CRC --taxonomy lineage --perturb-map
```
The final block must print *"No verdict changed under any perturbation."* If it
ever prints `VERDICT CHANGES`, that cluster's registry row is load-bearing and
must be resolved before the result is reported.

### Step 10 — the doublet-exclusion control
```bash
python run_verify.py --dataset CRC_doublets --taxonomy lineage --perturb-map
```
Compare against CRC. Identical verdicts ⇒ dropping the interface doublets did not
manufacture the result (§5.3).

### Step 11 — re-run the falsifiability gate
```bash
cd data_preprocessing
python run_ingest.py --dataset CRC --validate-only
```
Nothing under `processed/` is written. Then diff the regenerated
`datasets/CRC/lineage_evidence.csv` against the committed one — the `predicted`
and `margin` columns should be identical (fixed seed 1029, 300k-cell cap).
Confirm every `CONTRADICTED (accepted)` row has a real written justification in
the registry's `evidence_override`, not a placeholder.

---

## 9. Commands

```bash
# Phase 1 — ingest (once per dataset)
cd data_preprocessing
python run_ingest.py --dataset UPMC      # baseline, has survival
python run_ingest.py --dataset CRC       # survival-less — still ingests
python run_ingest.py --dataset CRC --validate-only   # checks only, writes nothing

# Phase 2 — verify (the common, survival-independent matrix)
cd ../spatial_positional_encoding
python run_verify.py --dataset CRC  --taxonomy lineage
python run_verify.py --dataset UPMC --taxonomy native
python run_verify.py --dataset CRC  --taxonomy lineage --perturb-map
python run_verify.py --dataset UPMC --label-col survival_status   # optional 4th test
python run_pipeline.py --list            # what's ingested, and which have survival

# Optional per-dataset survival downstream (only if the cohort has survival)
python run_survival.py
```

## 10. Key files

| file | role |
|---|---|
| `data_preprocessing/run_ingest.py` | Phase 1 entry point |
| `data_preprocessing/datasets/<NAME>/adapter_config.py` | per-dataset raw→canonical mapping |
| `data_preprocessing/celltype_registry.csv` | the cell-type → lineage mapping: cited, versioned, ontology-grounded |
| `data_preprocessing/registry.py` | loads + structurally validates the registry |
| `data_preprocessing/markers.py` | cross-cohort marker-name resolution (aliases, channel tags) |
| `data_preprocessing/lineage_evidence.py` | the falsifiability gate — registry vs marker data |
| `data_preprocessing/processor.py` / `validator.py` / `schema.py` | ingest engine (survival optional) |
| `spatial_positional_encoding/src/cohort.py` | loads a canonical cohort |
| `spatial_positional_encoding/src/spatial_features.py` | Delaunay graph + 5 enrichment scalars |
| `spatial_positional_encoding/run_verify.py` | Phase 2 — the verification matrix |
| `spatial_positional_encoding/run_pipeline.py` | orchestrator (ingest → verify) |
| `spatial_positional_encoding/run_survival.py` | optional per-dataset survival check |

---

## 11. Fixed — the single-file marker-column collision

**The defect (2026-07-20, now fixed).** CRC's adapter points `LOCATIONS_PATH` and
`EXPRESSION_PATH` at the *same* combined CSV, so both frames carried all 56 marker
columns. The expression view was subset to `[acq, cell] + marker_cols` before the
merge, but the locations view was passed in whole — so pandas suffixed the
collision into `CD44 - stroma:Cyc_2_ch_2_x` (**raw, pre-arcsinh**, from locations)
and `..._y` (normalised, from expression), 112 columns in total, while
`marker_columns.txt` recorded the *unsuffixed* names. Result: it named 56 columns
that existed in no parquet.

The one-sided subsetting was the whole bug. The validator's falsifiability gate
(lineage_evidence.py:132-138) already subsets *both* views correctly, and CRC's
adapter already documents the intended behaviour — *"the processor reads it three
times and each COLUMN_MAP pulls out the columns that view needs."* The processor
simply didn't implement the second half of that sentence.

**The fix** (processor.py, before the merge): drop from `locations` only those
columns that collide with `marker_cols`, and log the drop as a numbered step so it
appears in `processing_report.md` like every other transformation. Splitting the
raw CSV into three files was considered and rejected — it would have fixed CRC
while leaving the engine wrong for the next single-file cohort, and it would have
inserted an unlogged reshaping step upstream of the validator, which is exactly
what the two-report audit trail exists to prevent.

**After re-ingesting CRC and CRC_doublets:**

| cohort | markers listed | present in parquets | suffixed cols | total cols |
|---|---|---|---|---|
| UPMC | 39 | 39 ✅ | 0 | 53 (unchanged — never affected) |
| CRC | 56 | 56 ✅ | 0 | 157 → **101** |
| CRC_doublets | 56 | 56 ✅ | 0 | 157 → **101** |

Verified after the fix:

- the surviving marker values are the **normalised** copy, not the raw one —
  spot-checked exactly against `arcsinh(raw / 5)` on three markers
- every non-marker native column is preserved (`groups`, `neighborhood name`,
  `Region`, the `CD4+ICOS+`-style phenotype flags)
- cell/sample/patient counts unchanged: CRC 240,554 / 140 / 35
- **`verify_lineage.csv` is byte-for-byte identical before and after** on both
  cohorts — all 5 features × 6 metrics. Confirmed by diff, not assumed: the
  verification path reads only `X, Y, lineage, cluster_label`, so §5's numbers
  never depended on the marker columns.

**Still open (separate gap):** `groups` (the CLR/DII label) survives into the
sample parquets but not into `manifest.parquet`, which carries only
`acquisition_id, patient_id, n_cells, n_cell_types, survival_*`. Since
`separability()` reads its label from the manifest, the CRC CLR/DII test described
in §4.3 cannot currently run. Fixing it means carrying selected metadata columns
through to the manifest.
