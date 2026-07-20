# FYDP — Spatial features from in-situ proteomics imaging

Do the **spatial arrangement** of cells and their **celltype-conditioned marker
states** carry information beyond "how many of each cell type are present"? The
pipeline builds two feature blocks and then adversarially tests whether they are
real — on two independent cohorts, without needing survival data.

## Start here

| document | what it gives you |
|---|---|
| **[doc/RESULT_REPORT.md](doc/RESULT_REPORT.md)** | **the results** — what was run, what held up, what did not |
| [doc/PIPELINE_EXPLAINED.md](doc/PIPELINE_EXPLAINED.md) | how it works end to end, plus a step-by-step manual verification guide |
| [doc/VERIFICATION.md](doc/VERIFICATION.md) | the verification design and metric definitions |
| [doc/CELLTYPE_MAPPING.md](doc/CELLTYPE_MAPPING.md) | the cell-type registry and its ontology grounding |
| [doc/LITERATURE_REVIEW.md](doc/LITERATURE_REVIEW.md) | cohort provenance and verified citations |

## Structure

```
data_preprocessing/           PHASE 1 — raw cohort -> canonical format (any dataset)
   run_ingest.py                validate, then process; writes an audit trail
   celltype_registry.csv        cited, ontology-grounded cell-type -> lineage map
   datasets/<NAME>/             per-cohort adapter, raw/, and processed/ output

spatial_positional_encoding/  PHASE 2 — build features, verify they are real
   run_verify.py                enrichment block matrix (+ map perturbation)
   run_verify_nodes.py          node-marker block matrix
   run_survival.py              PHASE 3 — optional, cohorts with survival only

doc/                          reports
discarded/                    retired code and superseded results, kept traceable
```

> The `spatial_positional_encoding/` name is historical — Laplacian positional
> encoding was removed from the pipeline. The retired code is in
> [discarded/legacy_pipeline/](discarded/legacy_pipeline/).

## Quick start

```bash
cd data_preprocessing && python run_ingest.py --dataset CRC
cd ../spatial_positional_encoding
python run_pipeline.py --list
python run_verify.py       --dataset CRC
python run_verify_nodes.py --dataset CRC
```

Full reproduction commands: [doc/RESULT_REPORT.md](doc/RESULT_REPORT.md) §10.

## Cohorts

| | UPMC | CRC |
|---|---|---|
| tissue | head & neck (HNSCC) | colorectal (Schürch et al. 2020) |
| samples / cells | 308 / 2,061,102 | 140 / 240,554 |
| survival | yes — 103 events | **none** — still fully verified |

A cohort with no outcome data is not a second-class cohort here: all three
verification tests are label-free by design, which is what makes CRC usable as
external validation at all.
