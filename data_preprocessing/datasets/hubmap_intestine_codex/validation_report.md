# Validation report — hubmap_intestine_codex

**Result: READY** — 19 pass, 1 warn, 0 fail

A FAIL means the schema requirement in schema.py is not met and processing will refuse to run without `--force`. A WARN means processing can proceed but some downstream feature block (a specific marker, a cell-type mapping, sample size) will be degraded or unavailable for this cohort.

| Check | Status | Detail |
|---|---|---|
| file:locations | PASS | D:\files\data_preprocessing\datasets\hubmap_intestine_codex\raw\cell_locations.parquet |
| file:expression | PASS | D:\files\data_preprocessing\datasets\hubmap_intestine_codex\raw\cell_expression.parquet |
| file:metadata | PASS | D:\files\data_preprocessing\datasets\hubmap_intestine_codex\raw\sample_metadata.csv |
| columns:locations | PASS | all 5 canonical columns present |
| columns:expression | PASS | all 2 canonical columns present |
| columns:metadata | PASS | all 4 canonical columns present |
| dtype:X | PASS | numeric |
| dtype:Y | PASS | numeric |
| dtype:survival_day | PASS | numeric |
| dtype:survival_status | PASS | binary 0/1 |
| markers:present | PASS | 47 marker columns: ['MUC2', 'SOX9', 'MUC1', 'CD31', 'Synapto', 'CD49f', 'CD15', 'CHGA']... |
| ids:overlap | PASS | 66 acquisition_ids common to all 3 tables (locations=66, expression=66, metadata=66) |
| grouping:patient_id | PASS | 8 unique patients across 66 samples (needed for patient-grouped CV — a single patient_id would leak across folds) |
| sample_size | PASS | 0/66 samples below MIN_CELLS_PER_SAMPLE=50 (median cells/sample=34388) — small samples are dropped during processing |
| celltype:mapping | PASS | all 25 native types mapped to a lineage |
| markers:lineage_validation | PASS | 9/10 canonical lineage markers present: ['CD45', 'PanCK', 'Vimentin', 'aSMA', 'CD31', 'Podoplanin', 'CD3e', 'CD21', 'CD68'] (used to sanity-check CELLTYPE_MAP against expression, not required to proceed) |
| markers:node_features | PASS | 3/8 node-marker features reproducible on this cohort: ['tumor_ki67', 'tcell_cd45ro', 'apc_mac_hladr'] |
| markers:recommended_pair | WARN | this project's best config needs ['FoxP3', 'CD45RO'] (C=0.733, RESULT_REPORT.md Table 5); present here: ['CD45RO'] |
| missing:coords | PASS | 0 null X/Y values (rows dropped at processing) |
| missing:markers | PASS | 0 null marker values (rows dropped at processing) |