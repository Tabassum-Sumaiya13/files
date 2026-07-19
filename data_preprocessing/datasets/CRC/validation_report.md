# Validation report — CRC

**Result: READY** — 14 pass, 5 warn, 0 fail

A FAIL means the schema requirement in schema.py is not met and processing will refuse to run without `--force`. A WARN means processing can proceed but some downstream feature block (a specific marker, a cell-type mapping, sample size) will be degraded or unavailable for this cohort.

| Check | Status | Detail |
|---|---|---|
| file:locations | PASS | D:\Desktop\FYDP\FYDP final works\files\data_preprocessing\datasets\CRC\raw\CRC_clusters_neighborhoods_markers.csv |
| file:expression | PASS | D:\Desktop\FYDP\FYDP final works\files\data_preprocessing\datasets\CRC\raw\CRC_clusters_neighborhoods_markers.csv |
| file:metadata | PASS | D:\Desktop\FYDP\FYDP final works\files\data_preprocessing\datasets\CRC\raw\CRC_clusters_neighborhoods_markers.csv |
| columns:locations | PASS | all 5 canonical columns present |
| columns:expression | PASS | all 2 canonical columns present |
| columns:metadata | PASS | all 2 canonical columns present |
| dtype:X | PASS | numeric |
| dtype:Y | PASS | numeric |
| survival:present | WARN | no survival columns ['survival_day', 'survival_status'] — cohort ingests for feature work, but the survival downstream check won't run (fill METADATA_COLUMN_MAP if this cohort has survival data) |
| markers:present | PASS | 56 marker columns: ['CD44 - stroma:Cyc_2_ch_2', 'FOXP3 - regulatory T cells:Cyc_2_ch_3', 'CD8 - cytotoxic T cells:Cyc_3_ch_2', 'p53 - tumor suppressor:Cyc_3_ch_3', 'GATA3 - Th2 helper T cells:Cyc_3_ch_4', 'CD45 - hematopoietic cells:Cyc_4_ch_2', 'T-bet - Th1 cells:Cyc_4_ch_3', 'beta-catenin - Wnt signaling:Cyc_4_ch_4']... |
| ids:overlap | PASS | 140 acquisition_ids common to all 3 tables (locations=140, expression=140, metadata=140) |
| grouping:patient_id | PASS | 35 unique patients across 140 samples (needed for patient-grouped CV — a single patient_id would leak across folds) |
| sample_size | PASS | 0/140 samples below MIN_CELLS_PER_SAMPLE=50 (median cells/sample=1979) — small samples are dropped during processing |
| celltype:mapping | WARN | 25/29 native types mapped (17,831/258,385 cells unmapped) — unmapped types dropped at processing: ['dirt', 'immune cells / vasculature', 'tumor cells / immune cells', 'undefined'] |
| markers:lineage_validation | WARN | 0/10 canonical lineage markers present: [] (used to sanity-check CELLTYPE_MAP against expression, not required to proceed) |
| markers:node_features | WARN | 0/8 node-marker features reproducible on this cohort: (none) |
| markers:recommended_pair | WARN | this project's best config needs ['FoxP3', 'CD45RO'] (C=0.733, RESULT_REPORT.md Table 5); present here: (none) |
| missing:coords | PASS | 0 null X/Y values (rows dropped at processing) |
| missing:markers | PASS | 0 null marker values (rows dropped at processing) |