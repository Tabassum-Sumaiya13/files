# Validation report — UPMC

**Result: READY** — 20 pass, 0 warn, 0 fail

A FAIL means the schema requirement in schema.py is not met and processing will refuse to run without `--force`. A WARN means processing can proceed but some downstream feature block (a specific marker, a cell-type mapping, sample size) will be degraded or unavailable for this cohort.

| Check | Status | Detail |
|---|---|---|
| file:locations | PASS | D:\Desktop\FYDP\FYDP final works\files\data\raw\dataset_info\cell_locations_and_labels.csv |
| file:expression | PASS | D:\Desktop\FYDP\FYDP final works\files\data\raw\dataset_info\labeled_arcsinh_norm_data.csv |
| file:metadata | PASS | D:\Desktop\FYDP\FYDP final works\files\data\raw\dataset_info\sample_metadata.csv |
| columns:locations | PASS | all 5 canonical columns present |
| columns:expression | PASS | all 2 canonical columns present |
| columns:metadata | PASS | all 4 canonical columns present |
| dtype:X | PASS | numeric |
| dtype:Y | PASS | numeric |
| dtype:survival_day | PASS | numeric |
| dtype:survival_status | PASS | binary 0/1 |
| markers:present | PASS | 39 marker columns: ['CD117', 'CD11b', 'CD11c', 'CD134', 'CD14', 'CD15', 'CD152', 'CD16']... |
| ids:overlap | PASS | 308 acquisition_ids common to all 3 tables (locations=308, expression=308, metadata=379) |
| grouping:patient_id | PASS | 82 unique patients across 379 samples (needed for patient-grouped CV — a single patient_id would leak across folds) |
| sample_size | PASS | 0/308 samples below MIN_CELLS_PER_SAMPLE=50 (median cells/sample=5912) — small samples are dropped during processing |
| celltype:mapping | PASS | all 16 native types mapped to a lineage |
| markers:lineage_validation | PASS | 10/10 canonical lineage markers present: ['CD45', 'PanCK', 'Vimentin', 'aSMA', 'CD31', 'Podoplanin', 'CD3e', 'CD20', 'CD21', 'CD68'] (used to sanity-check CELLTYPE_MAP against expression, not required to proceed) |
| markers:node_features | PASS | 8/8 node-marker features reproducible on this cohort: ['tumor_ki67', 'tumor_mac_pdl1', 'cd8_granzymeb', 'cd4_pd1', 'cd4_icos', 'cd4_foxp3', 'tcell_cd45ro', 'apc_mac_hladr'] |
| markers:recommended_pair | PASS | this project's best config needs ['FoxP3', 'CD45RO'] (C=0.733, RESULT_REPORT.md Table 5); present here: ['FoxP3', 'CD45RO'] |
| missing:coords | PASS | 0 null X/Y values (rows dropped at processing) |
| missing:markers | PASS | 0 null marker values (rows dropped at processing) |