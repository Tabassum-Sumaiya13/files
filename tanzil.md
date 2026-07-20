# Tanzil Study Notes

## What this repo is about
This repository studies whether **cell spatial arrangement** and **celltype-conditioned marker states** contain signal beyond simple cell counts.

The pipeline builds two main feature blocks and checks whether they are meaningful on independent cohorts:
- spatial arrangement features
- marker-state / enrichment features

## Big picture
1. **Preprocess raw cohort data** into a canonical format.
2. **Build spatial and marker-based features** from the processed data.
3. **Verify** that the features are real using label-free tests.
4. Optionally run **survival analysis** when outcome data exists.

## Main folders
- `data_preprocessing/` - ingestion and cohort normalization
- `spatial_positional_encoding/` - feature generation and verification
- `doc/` - reports, summaries, and verification notes
- `discarded/` - retired or superseded code kept for traceability

## Files worth reading first
- `README.md` - best high-level overview
- `doc/PIPELINE_EXPLAINED.md` - full end-to-end explanation
- `doc/RESULT_REPORT.md` - what was run and what worked
- `doc/VERIFICATION.md` - test design and metric definitions
- `doc/CELLTYPE_MAPPING.md` - lineage and cell-type registry details

## Core ideas to understand
- Why spatial structure might matter even when cell counts are similar
- How cell types are mapped to a canonical registry
- How preprocessing turns different cohorts into one shared schema
- Why verification is label-free and useful even without survival data
- What makes a feature block biologically meaningful versus just noisy

## Suggested study order
1. Read `README.md`.
2. Read `doc/PIPELINE_EXPLAINED.md`.
3. Read `doc/VERIFICATION.md`.
4. Inspect `data_preprocessing/run_ingest.py`.
5. Inspect `spatial_positional_encoding/run_pipeline.py` and `run_verify.py`.
6. Skim the cohort-specific adapter configs under `data_preprocessing/datasets/`.

## Useful commands
```bash
cd data_preprocessing && python run_ingest.py --dataset CRC
cd ../spatial_positional_encoding
python run_pipeline.py --list
python run_verify.py --dataset CRC
python run_verify_nodes.py --dataset CRC
```

## Quick self-check questions
- What is the difference between preprocessing and verification?
- Why are there two feature blocks?
- What does the registry standardize?
- Why is CRC useful even without survival data?
- Which files define the cohort-specific input adapters?

## Short summary
If I remember only one thing: this repo is not just about making features, it is about proving the features still carry signal after preprocessing and across cohorts.
