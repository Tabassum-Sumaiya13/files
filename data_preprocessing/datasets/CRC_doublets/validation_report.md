# Validation report — CRC_doublets

**Result: READY** — 24 pass, 4 warn, 0 fail

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
| celltype:mapping | WARN | 27/29 native types mapped (13,881/258,385 cells unmapped) — unmapped types dropped at processing: ['dirt', 'undefined']; every exclusion has a recorded reason in the registry |
| celltype:rows | PASS | 29 rows (v1.0.0, fingerprint de43398abc3a): 27 mapped, 2 deliberately excluded |
| celltype:unique_labels | PASS | no duplicate native_label rows |
| celltype:lineage_vocabulary | PASS | all lineage values in ('immune', 'tumour', 'stromal') |
| celltype:cl_grounding | PASS | all 27 mapped rows carry a CL term whose anchor matches the declared lineage |
| celltype:citation | PASS | all 27 mapped rows carry a source citation |
| celltype:cl_verified | WARN | 27/27 CL terms not yet confirmed against the ontology (verified != yes). They are traceable but unproven — see doc/CELLTYPE_MAPPING.md for the one-time confirmation procedure |
| celltype:exclusion_reasons | PASS | all 2 excluded types record why: ['dirt', 'undefined'] |
| celltype:coverage | PASS | all 29 native types in the raw data have a registry row |
| markers:lineage_validation | PASS | 27/27 canonical lineage markers resolved (vocabulary v1.0.0): ['CD11b', 'CD11c', 'CD138', 'CD20', 'CD31', 'CD34', 'CD38', 'CD3e', 'CD4', 'CD44', 'CD45', 'CD45RO', 'CD56', 'CD57', 'CD68', 'CD8', 'CDX2', 'CollagenIV', 'EGFR', 'FoxP3', 'Ki67', 'MUC1', 'PanCK', 'Podoplanin', 'Vimentin', 'aSMA', 'p53'] |
| markers:node_features | PASS | 8/8 node-marker features reproducible on this cohort: ['apc_mac_hladr', 'cd4_foxp3', 'cd4_icos', 'cd4_pd1', 'cd8_granzymeb', 'tcell_cd45ro', 'tumor_ki67', 'tumor_mac_pdl1'] |
| markers:recommended_pair | PASS | this project's best config needs ['FoxP3', 'CD45RO'] (C=0.733, RESULT_REPORT.md Table 5); resolved here: ['CD45RO', 'FoxP3'] |
| celltype:marker_evidence | WARN | marker evidence CONTRADICTS the registry for 4/27 mapped clusters covering 3.9% of mapped cells. ACCEPTED with a recorded justification in the registry `evidence_override` (3.9% of cells): plasma cells (declared immune, evidence tumour, margin +1.84, 8,510 cells); CD11b+ monocytes (declared immune, evidence tumour, margin +0.41, 815 cells); CD4+ T cells GATA3+ (declared immune, evidence tumour, margin +0.53, 67 cells); CD163+ macrophages (declared immune, evidence stromal, margin +0.59, 38 cells). every contradiction, accepted or not, is carried into run_verify.py --perturb-map — acceptance is a reason, not a result [scored 258,385 cells over 29 native clusters; core markers used: {'immune': ['CD11b', 'CD11c', 'CD20', 'CD38', 'CD3e', 'CD45', 'CD56', 'CD68'], 'tumour': ['CDX2', 'PanCK'], 'stromal': ['CD31', 'CD34', 'CollagenIV', 'Podoplanin', 'aSMA']}] |
| missing:coords | PASS | 0 null X/Y values (rows dropped at processing) |
| missing:markers | PASS | 0 null marker values (rows dropped at processing) |

## Cell-type mapping provenance

- Cell-type registry: `celltype_registry.csv` v1.0.0, fingerprint `de43398abc3a`
- Marker vocabulary: `markers.py` v1.0.0
- Lineage vocabulary: ('immune', 'tumour', 'stromal')
- Lineage is derived from each row's Cell Ontology term via `schema.CL_LINEAGE_ANCHOR`, not chosen per row.

## Marker evidence per native cluster

Each cluster's mean marker profile, z-scored across clusters, scored against the three lineage core panels (`schema.LINEAGE_MARKER_PANELS`). `predicted` is the argmax; `margin` is top minus runner-up. **AGREE** = evidence supports the registry. **AMBIGUOUS** = margin too small to decide either way. **CONTRADICTED** = evidence favours a different lineage by a clear margin — the registry row must be justified in its `notes` or changed, and carried into `run_verify.py --perturb-map`.

| native_label | n_cells | declared | predicted | margin | z_immune | z_tumour | z_stromal | verdict |
|---|---|---|---|---|---|---|---|---|
| tumor cells | 47602 | tumour | tumour | 4.516 | -0.872 | 3.925 | -0.591 | AGREE |
| smooth muscle | 27817 | stromal | stromal | 3.187 | -0.877 | -0.464 | 2.723 | AGREE |
| granulocytes | 22144 | immune | immune | 1.384 | 1.284 | -0.101 | -0.509 | AGREE |
| stroma | 20139 | stromal | stromal | 0.36 | -0.997 | -0.482 | -0.122 | AGREE |
| B cells | 13043 | immune | immune | 2.578 | 2.028 | -0.579 | -0.551 | AGREE |
| vasculature | 11725 | stromal | stromal | 3.253 | -0.969 | -0.535 | 2.718 | AGREE |
| Tregs | 2791 | immune | immune | 0.43 | -0.025 | -0.485 | -0.456 | AGREE |
| CD4+ T cells | 2303 | immune | immune | 1.054 | 0.798 | -0.255 | -0.601 | AGREE |
| immune cells / vasculature | 2153 | stromal | stromal | 0.693 | -0.431 | 0.233 | 0.926 | AGREE |
| CD68+ macrophages | 2108 | immune | immune | 0.498 | 0.131 | -0.367 | -0.503 | AGREE |
| CD11b+CD68+ macrophages | 1500 | immune | immune | 0.797 | 0.715 | -0.225 | -0.082 | AGREE |
| nerves | 659 | stromal | stromal | 0.336 | 0.118 | -0.471 | 0.453 | AGREE |
| CD11c+ DCs | 400 | immune | immune | 2.311 | 1.708 | -0.607 | -0.603 | AGREE |
| lymphatics | 328 | stromal | stromal | 3.231 | -0.964 | -0.459 | 2.772 | AGREE |
| NK cells | 323 | immune | immune | 0.627 | 1.529 | 0.901 | -0.18 | AGREE |
| CD3+ T cells | 189 | immune | immune | 1.812 | 1.427 | -0.386 | -0.581 | AGREE |
| CD68+ macrophages GzmB+ | 183 | immune | immune | 1.685 | 1.572 | -0.113 | -0.433 | AGREE |
| CD68+CD163+ macrophages | 39596 | immune | immune | 0.106 | -0.108 | -0.214 | -0.416 | AMBIGUOUS |
| CD8+ T cells | 16675 | immune | immune | 0.114 | -0.178 | -0.292 | -0.418 | AMBIGUOUS |
| CD4+ T cells CD45RO+ | 16661 | immune | tumour | 0.009 | -0.218 | -0.209 | -0.477 | AMBIGUOUS |
| immune cells | 3127 | immune | tumour | 0.086 | -1.004 | -0.469 | -0.555 | AMBIGUOUS |
| adipocytes | 1811 | stromal | stromal | 0.206 | -1.01 | -0.555 | -0.349 | AMBIGUOUS |
| tumor cells / immune cells | 1797 | tumour | tumour | 0.084 | -1.023 | -0.527 | -0.611 | AMBIGUOUS |
| plasma cells | 8510 | immune | tumour | 1.84 | 1.172 | 3.013 | 0.28 | CONTRADICTED (accepted) |
| CD11b+ monocytes | 815 | immune | tumour | 0.408 | -1.004 | -0.192 | -0.6 | CONTRADICTED (accepted) |
| CD4+ T cells GATA3+ | 67 | immune | tumour | 0.529 | -0.315 | 0.214 | -0.583 | CONTRADICTED (accepted) |
| CD163+ macrophages | 38 | immune | stromal | 0.593 | -0.745 | -0.25 | 0.343 | CONTRADICTED (accepted) |
| dirt | 7357 | (excluded) | tumour | 0.533 | -0.924 | -0.043 | -0.576 | EXCLUDED |
| undefined | 6524 | (excluded) | tumour | 0.414 | -0.818 | -0.005 | -0.419 | EXCLUDED |

## Accepted contradictions

Contradictions kept deliberately, each with a written reason recorded in `celltype_registry.csv` (`evidence_override`). These are downgraded from blocking to reported — never hidden. They are still counted above and are still carried into the perturbation analysis, which is what decides whether any conclusion actually depends on them.

- **plasma cells** — ACCEPTED 2026-07-20. Evidence favours tumour (z=3.01) despite CD38 raising the immune score to +1.17. Plasma cells are unambiguously B-lineage on cell-identity grounds (CL:0000786), so the registry is not changed. The elevated epithelial signal is most consistent with segmentation spillover from adjacent CDX2+/PanCK+ epithelium, which plasma cells sit against in the lamina propria. Largest single contradiction in either cohort - carried into the perturbation analysis.
- **CD11b+ monocytes** — ACCEPTED 2026-07-20. Small margin (+0.41) with all three lineage scores negative - the cluster is dim across every core panel, so the argmax is weakly determined. Retained as immune per the cited source annotation.
- **CD4+ T cells GATA3+** — ACCEPTED 2026-07-20. 67 cells; per-cluster marker means are unstable at this size. Retained as immune per the cited source annotation.
- **CD163+ macrophages** — ACCEPTED 2026-07-20. 38 cells; per-cluster marker means are unstable at this size. Retained as immune per the cited source annotation.
- **immune cells / vasculature** — INCLUDED for the exclusion-bias test. Doublet spanning immune and vasculature; assigned stromal because that is what the marker evidence predicts (z_stromal +0.93). --perturb-map drops and flips it.
- **tumor cells / immune cells** — INCLUDED for the exclusion-bias test. Doublet spanning tumour and immune; assigned tumour because that is what the marker evidence predicts for it in the CRC report. There is no correct answer - that is the point. --perturb-map drops and flips it.