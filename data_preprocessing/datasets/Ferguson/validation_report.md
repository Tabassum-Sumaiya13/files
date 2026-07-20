# Validation report — Ferguson

**Result: READY** — 20 pass, 1 warn, 0 fail

A FAIL means the schema requirement in schema.py is not met and processing will refuse to run without `--force`. A WARN means processing can proceed but some downstream feature block (a specific marker, a cell-type mapping, sample size) will be degraded or unavailable for this cohort.

| Check | Status | Detail |
|---|---|---|
| file:locations | PASS | C:\Users\User\Downloads\files-main\files-main\data_preprocessing\datasets\Ferguson\raw\cell_locations.csv |
| file:expression | PASS | C:\Users\User\Downloads\files-main\files-main\data_preprocessing\datasets\Ferguson\raw\cell_expression.csv |
| file:metadata | PASS | C:\Users\User\Downloads\files-main\files-main\data_preprocessing\datasets\Ferguson\raw\sample_metadata.csv |
| columns:locations | PASS | all 5 canonical columns present |
| columns:expression | PASS | all 2 canonical columns present |
| columns:metadata | PASS | all 2 canonical columns present |
| dtype:X | PASS | numeric |
| dtype:Y | PASS | numeric |
| survival:present | PASS | survival_day + survival_status present |
| dtype:survival_day | PASS | numeric |
| dtype:survival_status | PASS | binary 0/1 |
| markers:present | PASS | 34 marker columns: ['panCK', 'CD20', 'HH3', 'CD45RA', 'CD8a', 'podoplanin', 'CD16', 'CADM1']... |
| ids:overlap | PASS | 44 acquisition_ids common to all 3 tables (locations=44, expression=44, metadata=44) |
| grouping:patient_id | PASS | 17 unique patients across 44 samples (needed for patient-grouped CV — a single patient_id would leak across folds) |
| sample_size | PASS | 0/44 samples below MIN_CELLS_PER_SAMPLE=50 (median cells/sample=3407) — small samples are dropped during processing |
| celltype:mapping | PASS | all 3 native types mapped to a lineage |
| markers:lineage_validation | WARN | 3/10 canonical lineage markers present: ['CD31', 'CD20', 'CD68'] (used to sanity-check CELLTYPE_MAP against expression, not required to proceed) |
| markers:node_features | PASS | 6/8 node-marker features reproducible on this cohort: ['tumor_ki67', 'tumor_mac_pdl1', 'cd4_pd1', 'cd4_icos', 'cd4_foxp3', 'tcell_cd45ro'] |
| markers:recommended_pair | PASS | this project's best config needs ['FoxP3', 'CD45RO'] (C=0.733, RESULT_REPORT.md Table 5); present here: ['FoxP3', 'CD45RO'] |
| missing:coords | PASS | 0 null X/Y values (rows dropped at processing) |
| missing:markers | PASS | 0 null marker values (rows dropped at processing) |