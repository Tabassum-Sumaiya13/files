# Feature verification — the dataset-agnostic pipeline (post-PE)

What the pipeline now
produces and checks is the **spatial enrichment feature block**, and it runs the
**same way on every cohort** — with or without survival.

## The flow

```
raw cohort
   │   Phase 1 — INGEST      (data_preprocessing/run_ingest.py)
   ▼
datasets/<NAME>/processed/          canonical, identical for every dataset
   ├─ samples/sample_<id>.parquet   X, Y, cluster_label, lineage, <markers>, (survival if any)
   ├─ manifest.parquet              acquisition_id, patient_id, survival_day/status (NaN if none)
   └─ marker_columns.txt
   │   Phase 2 — VERIFY      (run_verify.py)  ← the common check, no survival needed
   ▼
datasets/<NAME>/processed/verification/verify_<taxonomy>.csv
```

Survival is **optional**. A cohort with no survival table (e.g. CRC) still
ingests (survival columns become `NaN`) and is fully verifiable. Survival, when
present, is a **separate per-dataset downstream check** (`run_survival.py`,
RandomSurvivalForest + patient-grouped C-index) — deliberately not part of this
common path.

### Commands
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
```

## Taxonomy: two resolutions, one code path

Features are computed on a **label column**, so the same code runs on either:

- `--taxonomy lineage` — the 3-way `immune / tumour / stromal` column present in
  **every** ingested cohort. Portable and directly comparable across datasets.
- `--taxonomy native` — each dataset's own `cluster_label` (UPMC 16, CRC 29).
  Finer, but dimensionality/meaning differ per dataset.

The spatial graph is **Delaunay**, built on the fly from X/Y (parameter-free,
mean degree ≈ 6) — no prebuilt edge files.

## The feature block (5 abundance-corrected scalars)

Same definitions as the original `enrichment_features.py`, generalised off the
hardcoded 16 cluster-ids. Each is **0 under the random-mix null**:

| feature | meaning |
|---|---|
| `kl_mean` | how far each cell type's neighbourhood sits from the global mix |
| `kl_tumor` | how distinctive the tumour niche is |
| `self_enrich` | do cell types clump with their own kind beyond chance |
| `immune_tumor` | immune infiltration of tumour beyond chance |
| `stroma_tumor` | stromal/tumour interface beyond chance |

## The metrics — what, why, how

**Baseline = composition proportions** (fraction of each lineage/cell type per
sample). This is the abundance-only null of the field: *"spatial arrangement adds
nothing over how much of each cell type is present."* Every spatial feature must
be shown to carry something this baseline does not. It needs no labels and exists
for every dataset.

Three **universal, label-free** checks (the common matrix, always runs):

1. **Null z** — feature vs its own within-sample **spatial-shuffle null** (graph
   fixed, cell labels permuted `--n-perm` times). `|z| ≥ 2` ⇒ the value reflects
   real spatial organisation, not chance. *Answers: is it real?*
2. **Spatial-specific = 1 − R²(feature ~ composition)** (5-fold, R² clamped to
   [0,1]). High ⇒ the feature is **not** just re-encoding abundance. This is the
   label-free analog of "beats the celltype baseline." *Answers: does it add
   beyond the baseline?*
3. **Stability r** — split each sample's cells into two random halves, rebuild
   the graph on each, recompute, correlate across samples. High ⇒ reproducible,
   not a small-sample artifact. *Answers: is it trustworthy?*

Verdict per feature: **STRONG** = real (|z|≥2) **and** adds (spec≥0.5) **and**
stable (r≥0.5).

Optional **separability** add-on (only when a categorical label exists, via
`--label-col`): GroupKFold-by-patient classification AUC/accuracy comparing
`Composition baseline` vs `Baseline + spatial block` vs `Baseline + noise`
(width-matched). A real block beats both baseline and noise. This is where a
dataset's own biological label (e.g. CRC `groups` = CLR/DII) plugs in — it is
never survival.

## Results on the two ingested cohorts (taxonomy = lineage)

Both cohorts: **all 5 features STRONG.** Survival was used for neither.

| cohort | samples | survival | null-z range | spatial-specific | stability r |
|---|---|---|---|---|---|
| UPMC | 308 | yes (unused here) | 31 … 2141 | 0.55 … 1.00 | 0.97 … 0.99 |
| CRC  | 140 | none | 4.6 … 93 | 0.85 … 1.00 | 0.62 … 0.96 |

Per-feature CSVs: `datasets/<NAME>/processed/verification/verify_lineage.csv`.

## New / changed files

- `data_preprocessing/schema.py`, `validator.py`, `processor.py` — survival made
  **optional** (`METADATA_REQUIRED` drops survival; missing survival → WARN, not
  FAIL; manifest survival = `NaN`).
- `data_preprocessing/datasets/UPMC/adapter_config.py` — re-ingests the baseline
  through the same path.
- `spatial_positional_encoding/src/cohort.py` — loads a canonical cohort.
- `spatial_positional_encoding/src/spatial_features.py` — Delaunay + the 5
  enrichment scalars, generalised to any label column.
- `spatial_positional_encoding/run_verify.py` — the verification matrix.
- `spatial_positional_encoding/run_pipeline.py` — rewritten: PE removed, now the
  ingest→verify orchestrator.
