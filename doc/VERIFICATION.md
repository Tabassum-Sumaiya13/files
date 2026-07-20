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

# Phase 2 — verify the enrichment block (survival-independent)
cd ../spatial_positional_encoding
python run_verify.py --dataset CRC  --taxonomy lineage
python run_verify.py --dataset UPMC --taxonomy native

# Phase 2b — verify the node-marker block (also survival-independent)
python run_verify_nodes.py --dataset CRC
python run_verify_nodes.py --dataset UPMC

# Phase 3 — optional; exits cleanly on a survival-less cohort
python run_survival.py --dataset UPMC

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

## A second block, verified the same way — node markers

The enrichment block is not the only feature set. `run_verify_nodes.py` applies the
same three tests to the **node-marker block**: 8 celltype-conditioned functional
markers (the mean of one protein over *only* the cells of the type it is read in).

Two things differ, and both are deliberate:

- **`null_z` becomes `cond_z`.** Node features are per-sample means and involve no
  graph, so permuting cell-type labels makes the conditioning set a size-matched
  random draw and the null becomes ~the bulk mean. `cond_z` therefore asks *is this
  marker actually enriched in the cell type the feature is named for* — the premise
  conditioning rests on. It is **not** a spatial claim.
- **`spatial_specific` becomes `composition_specific`.** Identical computation,
  honest name: node features carry no spatial information to be specific about.

Conditioning is resolved per cohort **by Cell Ontology term** from
`celltype_registry.csv`, so `cd4_pd1` becomes `['CD4 T cell']` on UPMC and the four
CD4/Treg clusters on CRC with no hand mapping. Samples lacking the conditioning cell
type are reported as unsupported, never zero-filled.

A negative `cond_z` gets its own verdict, **CONTRADICTED** — the signal is real and
reproducible but says the marker is *depleted* in its own celltype, falsifying the
feature's premise rather than confirming it.

## Results

Current numbers for both blocks, both cohorts, both taxonomies, plus the optional
survival downstream, live in **[RESULT_REPORT.md](RESULT_REPORT.md)** — kept there
rather than duplicated here so there is one place to update.

Headline: enrichment **5/5 STRONG on both cohorts at both taxonomies**; node markers
**6/8 STRONG on both cohorts**; neither block's survival deltas exceed their own
fold-to-fold spread.

Per-cohort CSVs: `datasets/<NAME>/processed/verification/`.

## The active file set

| file | role |
|---|---|
| `data_preprocessing/schema.py`, `validator.py`, `processor.py` | ingest engine; survival is **optional** (missing → WARN, manifest survival = `NaN`) |
| `data_preprocessing/registry.py`, `markers.py`, `lineage_evidence.py` | the cell-type registry, marker-name resolution, and the falsifiability gate |
| `spatial_positional_encoding/src/cohort.py` | loads a canonical cohort |
| `spatial_positional_encoding/src/spatial_features.py` | Delaunay + the 5 enrichment scalars |
| `spatial_positional_encoding/src/node_features.py` | the 8 celltype-conditioned marker features |
| `spatial_positional_encoding/run_verify.py` | enrichment matrix (+ `--perturb-map`, `--label-col`) |
| `spatial_positional_encoding/run_verify_nodes.py` | node-marker matrix |
| `spatial_positional_encoding/run_survival.py` | optional survival downstream |
| `spatial_positional_encoding/run_pipeline.py` | orchestrator / dataset lister |

Anything else that used to live under `spatial_positional_encoding/` has been
retired to `discarded/legacy_pipeline/`, which documents what each module was and
what replaced it.
