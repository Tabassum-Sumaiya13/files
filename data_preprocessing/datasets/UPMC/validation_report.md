# Validation report — UPMC

**Result: READY** — 27 pass, 2 warn, 0 fail

A FAIL means the schema requirement in schema.py is not met and processing will refuse to run without `--force`. A WARN means processing can proceed but some downstream feature block (a specific marker, a cell-type mapping, sample size) will be degraded or unavailable for this cohort.

| Check | Status | Detail |
|---|---|---|
| file:locations | PASS | D:\Desktop\FYDP\FYDP final works\files\data\raw\dataset_info\cell_locations_and_labels.csv |
| file:expression | PASS | D:\Desktop\FYDP\FYDP final works\files\data\raw\dataset_info\labeled_arcsinh_norm_data.csv |
| file:metadata | PASS | D:\Desktop\FYDP\FYDP final works\files\data\raw\dataset_info\sample_metadata.csv |
| columns:locations | PASS | all 5 canonical columns present |
| columns:expression | PASS | all 2 canonical columns present |
| columns:metadata | PASS | all 2 canonical columns present |
| dtype:X | PASS | numeric |
| dtype:Y | PASS | numeric |
| survival:present | PASS | survival_day + survival_status present |
| dtype:survival_day | PASS | numeric |
| dtype:survival_status | PASS | binary 0/1 |
| markers:present | PASS | 39 marker columns: ['CD117', 'CD11b', 'CD11c', 'CD134', 'CD14', 'CD15', 'CD152', 'CD16']... |
| ids:overlap | PASS | 308 acquisition_ids common to all 3 tables (locations=308, expression=308, metadata=379) |
| grouping:patient_id | PASS | 82 unique patients across 379 samples (needed for patient-grouped CV — a single patient_id would leak across folds) |
| sample_size | PASS | 0/308 samples below MIN_CELLS_PER_SAMPLE=50 (median cells/sample=5912) — small samples are dropped during processing |
| celltype:mapping | PASS | all 16 native types mapped to a lineage |
| celltype:rows | PASS | 16 rows (v1.0.0, fingerprint de43398abc3a): 16 mapped, 0 deliberately excluded |
| celltype:unique_labels | PASS | no duplicate native_label rows |
| celltype:lineage_vocabulary | PASS | all lineage values in ('immune', 'tumour', 'stromal') |
| celltype:cl_grounding | PASS | all 16 mapped rows carry a CL term whose anchor matches the declared lineage |
| celltype:citation | PASS | all 16 mapped rows carry a source citation |
| celltype:cl_verified | WARN | 16/16 CL terms not yet confirmed against the ontology (verified != yes). They are traceable but unproven — see doc/CELLTYPE_MAPPING.md for the one-time confirmation procedure |
| celltype:coverage | PASS | all 16 native types in the raw data have a registry row |
| markers:lineage_validation | PASS | 21/27 canonical lineage markers resolved (vocabulary v1.0.0): ['CD11b', 'CD11c', 'CD20', 'CD31', 'CD34', 'CD38', 'CD3e', 'CD4', 'CD45', 'CD45RO', 'CD56', 'CD57', 'CD68', 'CD8', 'CollagenIV', 'FoxP3', 'Ki67', 'PanCK', 'Podoplanin', 'Vimentin', 'aSMA']; not found: ['CD138', 'CD44', 'CDX2', 'EGFR', 'MUC1', 'p53'] |
| markers:node_features | PASS | 8/8 node-marker features reproducible on this cohort: ['apc_mac_hladr', 'cd4_foxp3', 'cd4_icos', 'cd4_pd1', 'cd8_granzymeb', 'tcell_cd45ro', 'tumor_ki67', 'tumor_mac_pdl1'] |
| markers:recommended_pair | PASS | this project's best config needs ['FoxP3', 'CD45RO'] (C=0.733, RESULT_REPORT.md Table 5); resolved here: ['CD45RO', 'FoxP3'] |
| celltype:marker_evidence | WARN | marker evidence CONTRADICTS the registry for 3/16 mapped clusters covering 13.5% of mapped cells. ACCEPTED with a recorded justification in the registry `evidence_override` (13.5% of cells): Tumor (Podo+) (declared tumour, evidence stromal, margin +0.89, 18,945 cells); APC (declared immune, evidence tumour, margin +1.20, 11,419 cells); Naive immune cell (declared immune, evidence stromal, margin +1.51, 10,059 cells). every contradiction, accepted or not, is carried into run_verify.py --perturb-map — acceptance is a reason, not a result [scored 300,000 cells over 16 native clusters; core markers used: {'immune': ['CD11b', 'CD11c', 'CD20', 'CD38', 'CD3e', 'CD45', 'CD56', 'CD68'], 'tumour': ['PanCK'], 'stromal': ['CD31', 'CD34', 'CollagenIV', 'Podoplanin', 'aSMA']}; core markers not found in this panel (not scored): {'tumour': ['CDX2']}] |
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
| Tumor (Ki67+) | 39796 | tumour | tumour | 1.159 | -0.809 | 0.786 | -0.372 | AGREE |
| Tumor | 32948 | tumour | tumour | 1.01 | -0.964 | 0.151 | -0.859 | AGREE |
| Macrophage | 30743 | immune | immune | 1.627 | 1.302 | -0.598 | -0.325 | AGREE |
| CD4 T cell | 28570 | immune | immune | 1.559 | 0.953 | -0.937 | -0.606 | AGREE |
| CD8 T cell | 23361 | immune | immune | 0.952 | 0.547 | -0.405 | -0.596 | AGREE |
| Stromal / Fibroblast | 23075 | stromal | stromal | 1.617 | -0.915 | -1.195 | 0.702 | AGREE |
| B cell | 20021 | immune | immune | 2.089 | 1.536 | -1.171 | -0.553 | AGREE |
| Tumor (CD15+) | 18351 | tumour | tumour | 2.332 | -0.925 | 1.499 | -0.833 | AGREE |
| Tumor (CD21+) | 16108 | tumour | tumour | 2.736 | -0.831 | 1.958 | -0.778 | AGREE |
| Vessel | 10402 | stromal | stromal | 2.917 | -0.727 | -1.202 | 2.19 | AGREE |
| Granulocyte | 7988 | immune | immune | 1.342 | 1.244 | -0.098 | -0.618 | AGREE |
| Tumor (CD20+) | 6287 | tumour | tumour | 1.09 | 0.059 | 1.15 | -0.755 | AGREE |
| Lymph vessel | 1927 | stromal | stromal | 0.042 | 1.739 | 0.366 | 1.78 | AMBIGUOUS |
| Tumor (Podo+) | 18945 | tumour | stromal | 0.89 | -0.84 | 0.49 | 1.38 | CONTRADICTED (accepted) |
| APC | 11419 | immune | tumour | 1.202 | -0.8 | 0.502 | -0.7 | CONTRADICTED (accepted) |
| Naive immune cell | 10059 | immune | stromal | 1.512 | -0.57 | -1.294 | 0.943 | CONTRADICTED (accepted) |

## Accepted contradictions

Contradictions kept deliberately, each with a written reason recorded in `celltype_registry.csv` (`evidence_override`). These are downgraded from blocking to reported — never hidden. They are still counted above and are still carried into the perturbation analysis, which is what decides whether any conclusion actually depends on them.

- **APC** — ACCEPTED 2026-07-20. Evidence favours tumour (PanCK) over immune (CD45). 'APC' is a FUNCTIONAL grouping (professional antigen-presenting cell), not a marker-derived identity, and the source cohort's annotation is not reproducible from this panel alone. Retained as immune per the cited source annotation. Materiality is settled by the perturbation analysis (run_verify.py --perturb-map), not by this note: if a conclusion moves when APC is reassigned, this row must be revisited.
- **Naive immune cell** — ACCEPTED 2026-07-20. Evidence favours stromal, driven by high CD31 (endothelial) - consistent with endothelial contamination of a generically-named cluster ('naive immune cell' names no subset). Retained as immune per the cited source annotation. Carried into the perturbation analysis.
- **Tumor (Podo+)** — ACCEPTED 2026-07-20. Genuinely ambiguous by construction: this cluster is PanCK-positive AND Podoplanin-high, and Podoplanin is simultaneously the lymphatic-endothelial anchor in the stromal panel. The max-based stromal score therefore wins on a marker the cluster is NAMED for. Retained as tumour per the cited source annotation. Carried into the perturbation analysis.