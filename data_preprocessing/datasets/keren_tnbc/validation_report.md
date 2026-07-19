# Validation report — keren_tnbc

**Result: READY** — 20 pass, 1 warn, 0 fail

A FAIL means the schema requirement in schema.py is not met and processing will refuse to run without `--force`. A WARN means processing can proceed but some downstream feature block (a specific marker, a cell-type mapping, sample size) will be degraded or unavailable for this cohort.

| Check | Status | Detail |
|---|---|---|
| file:locations | PASS | C:\Users\MARZIA ISLAM\Downloads\files\data_preprocessing\datasets\keren_tnbc\raw\cell_locations.csv |
| file:expression | PASS | C:\Users\MARZIA ISLAM\Downloads\files\data_preprocessing\datasets\keren_tnbc\raw\cell_expression.csv |
| file:metadata | PASS | C:\Users\MARZIA ISLAM\Downloads\files\data_preprocessing\datasets\keren_tnbc\raw\sample_metadata.csv |
| columns:locations | PASS | all 5 canonical columns present |
| columns:expression | PASS | all 2 canonical columns present |
| columns:metadata | PASS | all 2 canonical columns present |
| dtype:X | PASS | numeric |
| dtype:Y | PASS | numeric |
| survival:present | PASS | survival_day + survival_status present |
| dtype:survival_day | PASS | numeric |
| dtype:survival_status | PASS | binary 0/1 |
| markers:present | PASS | 49 marker columns: ['C', 'Na', 'Si', 'P', 'Ca', 'Fe', 'dsDNA', 'Vimentin']... |
| ids:overlap | PASS | 40 acquisition_ids common to all 3 tables (locations=40, expression=40, metadata=40) |
| grouping:patient_id | PASS | 40 unique patients across 40 samples (needed for patient-grouped CV — a single patient_id would leak across folds) |
| sample_size | PASS | 0/40 samples below MIN_CELLS_PER_SAMPLE=50 (median cells/sample=4961) — small samples are dropped during processing |
| celltype:mapping | WARN | 16/17 native types mapped (1,725/197,678 cells unmapped) — unmapped types dropped at processing: ['Unidentified'] |
| markers:lineage_validation | PASS | 5/10 canonical lineage markers present: ['CD45', 'Vimentin', 'CD31', 'CD20', 'CD68'] (used to sanity-check CELLTYPE_MAP against expression, not required to proceed) |
| markers:node_features | PASS | 5/8 node-marker features reproducible on this cohort: ['tumor_ki67', 'cd4_pd1', 'cd4_foxp3', 'tcell_cd45ro', 'apc_mac_hladr'] |
| markers:recommended_pair | PASS | this project's best config needs ['FoxP3', 'CD45RO'] (C=0.733, RESULT_REPORT.md Table 5); present here: ['FoxP3', 'CD45RO'] |
| missing:coords | PASS | 0 null X/Y values (rows dropped at processing) |
| missing:markers | PASS | 0 null marker values (rows dropped at processing) |