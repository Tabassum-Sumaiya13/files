# data_preprocessing/

Ingest an **external** spatial-proteomics cohort (CODEX, IMC, MIBI) — anything
with per-cell XY coordinates, protein markers, cell-type labels, (a linked
survival table if exist ) and turn it into the exact schema
`spatial_positional_encoding/` expects, so its enrichment / node-marker /
survival code can be evaluated on a second cohort for external validation
(the gap flagged in `doc/RESULT_REPORT.md` §9 — every current result comes
from one cohort).

## Two-step process

1. **Validate** (`validator.py`) — checks the raw files against
   [`schema.py`](schema.py): required columns, dtypes, ID consistency across
   the three tables, patient grouping, marker-panel overlap, cell-type
   mapping coverage. **Nothing is changed on disk.** Writes
   `datasets/<name>/validation_report.md`.

2. **Process** (`processor.py`) — only runs if validation is READY (zero FAIL
   checks) or `--force` is passed. Cleans, deduplicates, arcsinh-normalises,
   maps cell types to the immune/tumour/stromal lineages, z-scores
   coordinates per sample, and exports one parquet per sample in the
   canonical format. Writes `datasets/<name>/processing_report.md` — every
   step performed, in order, with before/after row counts, so the exact
   transformation from raw file to processed cohort is auditable later.

## Adding a new dataset

```bash
# 1. Copy the template
cp -r datasets/_template datasets/<your_dataset_name>

# 2. Drop the raw files
#    datasets/<your_dataset_name>/raw/cell_locations.csv
#    datasets/<your_dataset_name>/raw/cell_expression.csv
#    datasets/<your_dataset_name>/raw/sample_metadata.csv

# 3. Fill in datasets/<your_dataset_name>/adapter_config.py
#    - column-name mapping (native -> canonical)
#    - cell-type -> lineage mapping
#    - arcsinh cofactor

# 4. Run it
python run_ingest.py --dataset <your_dataset_name>
```

Read `validation_report.md` first. Every **FAIL** blocks processing —
usually a missing required column or a broken ID join, fixed by editing
`adapter_config.py`. **WARN**s are fine to proceed with, but note what they
cost you:

- missing `FoxP3` / `CD45RO` → the best-known config on this project
  (`Celltype + Enrichment + Markers(2)`, C = 0.733, RESULT_REPORT.md Table 5)
  can't be reproduced exactly, but the enrichment/`kl_mean` features (which
  only need cell-type labels, not specific markers) still can
- unmapped cell types → those cells are dropped, shrinking the cohort
- samples below `MIN_CELLS_PER_SAMPLE` → dropped, same as the baseline
  pipeline's own QC step

## Output layout

```
datasets/<name>/
├── raw/                     # your input files (gitignored)
├── adapter_config.py        # the mapping you filled in
├── validation_report.md     # what schema.py found, before any changes
├── processing_report.md     # every change made, in order, with counts
└── processed/
    ├── samples/sample_<acquisition_id>.parquet
    ├── manifest.parquet     # 1 row/sample: patient_id, n_cells, survival_day/status, file_path
    └── marker_columns.txt
```

`processed/samples/` is in the same shape as
`spatial_positional_encoding/data/processed/samples/` — point
`run_survival.py` at it (or merge it in) to evaluate the existing
enrichment / node-marker features on the new cohort.

## Files

| File | Purpose |
|---|---|
| `schema.py` | Canonical column names, cell-type taxonomy, and the marker lists used to score how much of the existing feature set a new cohort can reproduce |
| `validator.py` | Schema/consistency checks — read-only |
| `processor.py` | Clean, normalise, harmonise, export — only runs on a validated cohort |
| `report.py` | Shared `ValidationReport` / `ChangeLog` classes that render the two markdown reports |
| `run_ingest.py` | CLI: `python run_ingest.py --dataset <name> [--validate-only] [--force]` |
| `datasets/_template/adapter_config.py` | Fill-in-the-blank config for a new cohort |

## Candidate datasets (from this project's own literature review)

See the assistant's prior answer / `doc/RESULT_REPORT.md` §10-11 for the
reasoning — Schürch et al. 2020 (CODEX, colorectal — source of the marker
shortlist), Jackson et al. 2020 and Danenberg et al. 2022 (IMC, breast —
larger event counts, useful for the "underpowered OS endpoint" limitation),
and Keren et al. 2018 (MIBI-TOF, breast — same modality as the current
cohort, lowest adapter friction).
