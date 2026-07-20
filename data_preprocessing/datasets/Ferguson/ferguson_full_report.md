# Ferguson Dataset — Full Pipeline Report

**Dataset**: Ferguson et al. 2022 (IMC, head & neck squamous cell carcinoma)
**Source**: Bioconductor `SpatialDatasets::spe_Ferguson_2022()`
**Date**: 2026-07-20

---

## 1. Scripts Executed

| Step | Script | Directory | Purpose | Run? |
|---|---|---|---|---|
| 1 | `schema.py` | `data_preprocessing/` | Canonical column/taxonomy definitions | Yes (imported by run_ingest.py) |
| 2 | `report.py` | `data_preprocessing/` | ValidationReport + ChangeLog classes | Yes (imported by run_ingest.py) |
| 3 | `validator.py` | `data_preprocessing/` | 21 schema checks | Yes (imported by run_ingest.py) |
| 4 | `processor.py` | `data_preprocessing/` | Clean, normalise, harmonise, export | Yes (imported by run_ingest.py) |
| 5 | `run_ingest.py` | `data_preprocessing/` | CLI orchestrator for steps 1–4 | Yes |
| 6 | `config.py` | `spatial_positional_encoding/` | Global config singleton (paths, params) | Not directly — imported by src modules during verification |
| 7 | `run_verify.py` | `spatial_positional_encoding/` | Standalone feature verification | Yes |
| 8 | `run_pipeline.py` | `spatial_positional_encoding/` | Orchestrator → delegates to run_verify.py | Yes |
| 9 | `run_survival.py` | `spatial_positional_encoding/` | RSF survival evaluation (patient-grouped C-index) | **Not run** — Ferguson has no survival data (all NaN) |

---

## 2. Dataset Overview

| Metric | Value |
|---|---|
| Modality | Imaging Mass Cytometry (IMC) |
| Tissue | Head & neck squamous cell carcinoma (HNSCC) |
| Total cells | 155,913 |
| Samples (acquisitions) | 44 |
| Patients | 17 |
| Marker channels | 34 (excluding DNA1, DNA2) |
| Cell types (lineage) | 3 — immune (54.3%), stromal (31.6%), tumour (14.2%) |
| Survival data | Not available |

### 2.1 Sample Size Distribution

| Statistic | Cells/Sample |
|---|---|
| Min | 1,042 |
| 25th percentile | 2,455 |
| Median | 3,407 |
| 75th percentile | 4,629 |
| Max | 6,498 |

### 2.2 Patient Structure

| Samples/Patient | Count |
|---|---|
| 1 | 2 |
| 2 | 4 |
| 3 | 10 |
| 4 | 1 |

### 2.3 Marker Panel (34 channels)

```
panCK, CD20, HH3, CD45RA, CD8a, podoplanin, CD16, CADM1, IDO, PDL1,
CD13, CD68, VISTA, CD31, CXCR3, pSTAT3, CCR7, CD14, FX111A, FoxP3,
PD1, CD45RO, OX40, NFKBp65, CD66a, Ki67, LAG3, CD3, granzB, PDL2,
CD4, HLADR, ICOS, TIM3
```

---

## 3. Phase 1 — Preprocessing

### 3.1 Schema Validation

**Result: READY** — 20 pass, 1 warn, 0 fail

| Check | Status | Detail |
|---|---|---|
| file:locations | PASS | `raw/cell_locations.csv` exists |
| file:expression | PASS | `raw/cell_expression.csv` exists |
| file:metadata | PASS | `raw/sample_metadata.csv` exists |
| columns:locations | PASS | All 5 canonical columns present |
| columns:expression | PASS | All 2 canonical columns present |
| columns:metadata | PASS | All 2 canonical columns present |
| dtype:X | PASS | Numeric |
| dtype:Y | PASS | Numeric |
| survival:present | PASS | `survival_day` + `survival_status` present (all NaN) |
| dtype:survival_day | PASS | Numeric |
| dtype:survival_status | PASS | Binary 0/1 |
| markers:present | PASS | 34 marker columns detected |
| ids:overlap | PASS | 44 acquisition_ids common to all 3 tables |
| grouping:patient_id | PASS | 17 unique patients across 44 samples |
| sample_size | PASS | 0/44 samples below MIN_CELLS_PER_SAMPLE=50 |
| celltype:mapping | PASS | All 3 native types mapped to a lineage |
| markers:lineage_validation | **WARN** | 3/10 canonical lineage markers present: CD31, CD20, CD68 |
| markers:node_features | PASS | 6/8 node-marker features reproducible |
| markers:recommended_pair | PASS | FoxP3 + CD45RO both present |
| missing:coords | PASS | 0 null X/Y values |
| missing:markers | PASS | 0 null marker values |

### 3.2 Processing Steps

| Step | Input → Output | Notes |
|---|---|---|
| load | 155,913 rows (loc + expr), 44 metadata rows | 34 marker columns |
| dedup:locations | 155,913 → 155,913 | 0 duplicates |
| dedup:expression | 155,913 → 155,913 | 0 duplicates |
| drop_missing:coords | 155,913 → 155,913 | 0 dropped |
| drop_missing:markers | 155,913 → 155,913 | 0 dropped |
| normalise_markers | — | arcsinh(x / 1.0) on 34 markers |
| merge_locations_expression | 155,913 → 155,913 | 0 lost |
| filter_small_samples | 44 → 44 | 0 dropped |
| celltype_to_lineage | 155,913 mapped | immune=84,603; stromal=49,231; tumour=22,079 |
| normalise_coords | — | Per-sample z-score on X, Y |
| merge_metadata | 44 → 44 | 0 lost; survival absent |
| export | 44 parquet files + manifest | |

### 3.3 Output Files

```
datasets/Ferguson/processed/
├── samples/                (44 parquet files, one per sample)
├── manifest.parquet        (patient_id, n_cells, file_path)
├── marker_columns.txt      (34 marker names)
└── verification/
    └── verify_lineage.csv  (feature verification matrix)
```

Each parquet: `acquisition_id`, `cell_id`, `X` (z-scored), `Y` (z-scored), `cluster_label`, 34 markers (arcsinh-normalised), `lineage`, `patient_id`.

---

## 4. Phase 2 — Feature Verification

Taxonomy: **lineage** (3-way: immune, tumour, stromal)
Graph type: **Delaunay** (parameter-free planar triangulation)
Baseline: **Composition proportions** (fraction of each lineage per sample)

### 4.1 Verification Metrics

| Metric | Definition |
|---|---|
| null_z_median | Spatial signal — feature vs spatial-shuffle null; \|z\| large = real organisation |
| spatial_specific | 1 − R² vs composition baseline; ≥0.5 = adds beyond abundance |
| stability_r | Split-half reproducibility; ≥0.5 = reproducible |
| frac \|z\|>2 | Fraction of samples with significant spatial signal |

### 4.2 Verification Matrix

| Feature | Real Mean | Null z | Spatial Specific | Stability r | Verdict |
|---|---|---|---|---|---|
| kl_mean | +0.188 | +187.2 | 1.000 | 0.815 | **STRONG** |
| kl_tumor | +0.319 | +254.7 | 1.000 | 0.682 | **STRONG** |
| self_enrich | +1.142 | +12.1 | 1.000 | 0.540 | **STRONG** |
| immune_tumor | −0.422 | −6.8 | 1.000 | 0.743 | **STRONG** |
| stroma_tumor | −0.440 | −3.6 | 0.819 | 0.726 | **STRONG** |

### 4.3 Feature Descriptions

| Feature | Meaning |
|---|---|
| kl_mean | Mean KL divergence from null across all enrichment features |
| kl_tumor | KL divergence for tumour neighbourhood enrichment |
| self_enrich | Self-enrichment — cells surrounded by same lineage |
| immune_tumor | Immune–tumour spatial interaction (negative = avoidance) |
| stroma_tumor | Stromal–tumour spatial interaction (negative = avoidance) |

**All 5 spatial enrichment features pass verification** — every feature shows real spatial organisation, adds information beyond composition, and is reproducible across split halves.

---

## 5. Node-Marker Feature Coverage

6 of 8 node-marker features are reproducible on this cohort:

| Feature | Marker | Status |
|---|---|---|
| tumor_ki67 | Ki67 | Available |
| tumor_mac_pdl1 | PDL1 | Available |
| cd4_pd1 | PD1 | Available |
| cd4_icos | ICOS | Available |
| cd4_foxp3 | FoxP3 | Available |
| tcell_cd45ro | CD45RO | Available |
| cd8_granzymeb | GranzymeB | **Missing** (`granzB` naming mismatch) |
| apc_mac_hladr | HLA-DR | **Missing** (`HLADR` naming mismatch) |

The 2 missing features are fixable by adding rename entries to `EXPRESSION_COLUMN_MAP` in `adapter_config.py`.

---

## 6. Survival Evaluation

`run_survival.py` was **not run** on this cohort. The Ferguson dataset contains no time-to-event survival data — both `survival_day` and `survival_status` are NaN across all 44 samples. The survival downstream check (RandomSurvivalForest + patient-grouped C-index) is not applicable.

---

## 7. Usability Verdict

**VERDICT: USABLE for spatial feature validation, NOT for survival prediction.**

The Ferguson dataset is a clean, successfully ingested external cohort that validates the spatial enrichment feature block on an independent dataset. Specifically:

- **Usable for**: Enrichment features (`kl_mean`, `kl_tumor`, `self_enrich`, `immune_tumor`, `stroma_tumor`), node-marker features (6/8), and any cell-type-proportion analysis at the lineage level (3-way).
- **Not usable for**: Survival prediction (`run_survival.py`), patient-outcome modelling, or any analysis requiring the 16-subtype UPMC taxonomy (Ferguson only has 3 lineages).
- **Caveats**: 2 node-marker features blocked by naming mismatches (fixable). Lineage validation limited (3/10 canonical markers). Coarse cell-type resolution reduces granularity of proportion-based features.

---

## 8. Summary

### What worked
- **Clean ingestion**: No duplicates, no missing data, all 44 samples retained
- **All 5 spatial features pass**: Real signal, adds beyond abundance, reproducible
- **Cross-cohort generalisability**: Features validated on external (non-UPMC) cohort
- **FoxP3 + CD45RO present**: Best-known config (C=0.733) can be reproduced

### Known Limitations
1. **No survival data** — survival downstream check not applicable
2. **Coarse cell-type labels** — 3 lineages (not 16 subtypes); less granular proportion features
3. **2 node-marker features unavailable** — `cd8_granzymeb` and `apc_mac_hladr` blocked by marker naming (`granzB` vs `GranzymeB`, `HLADR` vs `HLA-DR`); fixable via adapter column map
4. **Lineage validation limited** — only 3/10 canonical markers present for independent cross-check

---

*Generated: 2026-07-20*
