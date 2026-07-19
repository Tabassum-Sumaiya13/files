# Processing report — hubmap_intestine_codex

What changed, in order, from the raw input files to the processed cohort now sitting in `processed/`. Re-run `python run_ingest.py --dataset <name>` any time to regenerate this after editing `adapter_config.py`.

| Step | Detail |
|---|---|
| load | locations=2,603,217 rows, expression=2,603,217 rows, metadata=66 rows, 47 marker columns detected |
| dedup:locations | 2,603,217 -> 2,603,217 rows (0 duplicate (sample, cell) pairs dropped) |
| dedup:expression | 2,603,217 -> 2,603,217 rows (0 duplicate (sample, cell) pairs dropped) |
| drop_missing:coords | 2,603,217 -> 2,603,217 rows (0 dropped for missing X/Y) |
| drop_missing:markers | 2,603,217 -> 2,603,217 rows (0 dropped for missing marker values) |
| normalise_markers | skipped — adapter_config.APPLY_ARCSINH=False, data already normalised |
| merge_locations_expression | 2,603,217 location rows -> 2,603,217 merged rows (0 location rows had no matching expression row and were dropped) |
| filter_small_samples | 66 -> 66 samples (dropped 0 sample(s) below MIN_CELLS_PER_SAMPLE=50: []) |
| celltype_to_lineage | all 2,603,217 cells mapped to a lineage: {'stromal': 972763, 'tumour': 928849, 'immune': 701605} |
| normalise_coords | per-sample z-score applied to X, Y (matches spatial_positional_encoding/src/preprocess.py Step 5) |
| merge_metadata | 66 -> 66 samples kept a survival label (0 sample(s) had no metadata row and were dropped) |
| export | 66 sample parquet files -> D:\files\data_preprocessing\datasets\hubmap_intestine_codex\processed\samples, manifest -> D:\files\data_preprocessing\datasets\hubmap_intestine_codex\processed\manifest.parquet |
| final_cohort | 2,603,217 cells, 66 samples, 8 patients, 0 events (0.0% event rate) |