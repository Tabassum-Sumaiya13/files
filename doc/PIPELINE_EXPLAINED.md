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
```

The point of Phase 1 is that Phase 2 never has to know which dataset it is looking
at. UPMC and CRC arrive at Phase 2 in the exact same shape.

Positional encoding (Laplacian PE) has been **removed** from the pipeline; the
spatial feature block is the abundance-corrected **enrichment** readout.

---

## 2. Phase 1 — Ingest: raw → canonical

Driven by a small per-dataset `adapter_config.py` that says "in *my* raw file, the
cell-type column is called *this*, the patient column is called *that*." The
engine (`data_preprocessing/processor.py`) then does the same sequence for every
dataset:

1. **Rename** raw columns → canonical names (`acquisition_id, cell_id, X, Y, cluster_label`).
2. **Dedup** repeated `(sample, cell)` rows.
3. **Drop** cells with missing coordinates / markers.
4. **Normalize markers** — `arcsinh(x / 5)` compression (skipped if the data is
   already normalized, as for UPMC).
5. **Map cell types → lineage** — every native `cluster_label` becomes one of
   `immune / tumour / stromal`. Unmapped types are dropped (e.g. CRC's
   `dirt` / `undefined` / doublet clusters), with a warning listing them.
6. **Z-score coordinates** per sample (position in comparable units across samples).
7. **Attach patient_id and survival** — survival is OPTIONAL and becomes `NaN`
   if the dataset has none.
8. **Write** one parquet per sample + a manifest + the marker list.

Result on the two datasets:

| | UPMC | CRC |
|---|---|---|
| samples | 308 | 140 |
| cells | 2,061,102 | 240,554 |
| patients | 81 | 35 |
| native cell types | 16 | 29 (25 kept) |
| survival | yes (103 events) | **none** |

### New cell types

New cell types are handled safely at every layer, but not kept automatically:

- **Ingest** — a native type not in the adapter's `CELLTYPE_MAP` is **dropped**
  (validator WARNs and lists it; processor logs the drop). It never crashes. To
  KEEP a new type, add it to that dataset's `CELLTYPE_MAP` with its lineage.
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

For each feature, `run_verify.py` computes three independent numbers. This is the
heart of it.

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

**Optional 4th test** (only when a dataset has a categorical label, e.g. CRC's
CLR/DII `groups`): GroupKFold-by-patient classification AUC, comparing
`Composition baseline` vs `Baseline + spatial block` vs `Baseline + noise`. A real
block beats both baseline and noise. This never uses survival.

---

## 5. Reading the results

### UPMC — lineage (`datasets/UPMC/processed/verification/verify_lineage.csv`)

| feature | null_z | spatial_specific | stability_r | meaning |
|---|---|---|---|---|
| kl_mean | **1398** | 0.77 | 0.98 | wildly non-random, mostly spatial, rock-solid |
| kl_tumor | **2141** | 0.93 | 0.99 | tumour niche is extremely distinctive |
| self_enrich | 31 | 0.55 | 0.98 | cells clump by type, reproducibly |
| immune_tumor | −52 | **1.00** | 0.99 | immune–tumour contact is *entirely* spatial (composition explains 0%) |
| stroma_tumor | −52 | 0.83 | 0.97 | stromal–tumour interface is real and specific |

### CRC — lineage, NO survival (`datasets/CRC/processed/verification/verify_lineage.csv`)

| feature | null_z | spatial_specific | stability_r | meaning |
|---|---|---|---|---|
| kl_mean | 81 | 0.94 | 0.96 | real, specific, stable |
| kl_tumor | 93 | 0.99 | 0.96 | tumour niche distinctive |
| self_enrich | 7.6 | 1.00 | 0.62 | clumping real & purely spatial; noisier split-half |
| immune_tumor | −9.4 | 0.87 | 0.88 | immune exclusion, real |
| stroma_tumor | −4.6 | 0.85 | 0.80 | real, weakest of the five |

### UPMC — native, 16 types (`verify_native.csv`)

Same story, and `self_enrich` jumps to **2.21** — at fine resolution, cells of a
*specific* subtype self-segregate even harder than at the coarse lineage level. A
sanity check that finer labels reveal finer structure, exactly as they should.

**How to read a row in one sentence:** *"This feature is `null_z`
standard-deviations away from random, composition can't explain `spatial_specific`
of it, and it reproduces at `stability_r` on independent halves of the tissue."*

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

---

## 7. Where survival fits now

Survival is no longer part of proving the features work — it is a **separate,
optional, per-dataset downstream question**: *"do these validated features also
predict patient outcome?"* (`spatial_positional_encoding/run_survival.py`,
RandomSurvivalForest + patient-grouped C-index). A dataset with survival gets that
as a bonus; a dataset without it (CRC) is still fully verified by Phase 2.

---

## 8. Commands

```bash
# Phase 1 — ingest (once per dataset)
cd data_preprocessing
python run_ingest.py --dataset UPMC      # baseline, has survival
python run_ingest.py --dataset CRC       # survival-less — still ingests

# Phase 2 — verify (the common, survival-independent matrix)
cd ../spatial_positional_encoding
python run_verify.py --dataset CRC  --taxonomy lineage
python run_verify.py --dataset UPMC --taxonomy native
python run_pipeline.py --list            # what's ingested, and which have survival

# Optional per-dataset survival downstream (only if the cohort has survival)
python run_survival.py
```

## 9. Key files

| file | role |
|---|---|
| `data_preprocessing/run_ingest.py` | Phase 1 entry point |
| `data_preprocessing/datasets/<NAME>/adapter_config.py` | per-dataset raw→canonical mapping |
| `data_preprocessing/processor.py` / `validator.py` / `schema.py` | ingest engine (survival optional) |
| `spatial_positional_encoding/src/cohort.py` | loads a canonical cohort |
| `spatial_positional_encoding/src/spatial_features.py` | Delaunay graph + 5 enrichment scalars |
| `spatial_positional_encoding/run_verify.py` | Phase 2 — the verification matrix |
| `spatial_positional_encoding/run_pipeline.py` | orchestrator (ingest → verify) |
| `spatial_positional_encoding/run_survival.py` | optional per-dataset survival check |
