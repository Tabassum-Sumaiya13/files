# Processing report — CRC_doublets

What changed, in order, from the raw input files to the processed cohort now sitting in `processed/`. Re-run `python run_ingest.py --dataset <name>` any time to regenerate this after editing `adapter_config.py`.

| Step | Detail |
|---|---|
| load | locations=258,385 rows, expression=258,385 rows, metadata=258,385 rows, 56 marker columns detected |
| dedup:locations | 258,385 -> 258,385 rows (0 duplicate (sample, cell) pairs dropped) |
| dedup:expression | 258,385 -> 258,385 rows (0 duplicate (sample, cell) pairs dropped) |
| drop_missing:coords | 258,385 -> 258,385 rows (0 dropped for missing X/Y) |
| drop_missing:markers | 258,385 -> 258,385 rows (0 dropped for missing marker values) |
| normalise_markers | arcsinh(x / 5.0) applied to 56 marker columns |
| merge_locations_expression | 258,385 location rows -> 258,385 merged rows (0 location rows had no matching expression row and were dropped) |
| filter_small_samples | 140 -> 140 samples (dropped 0 sample(s) below MIN_CELLS_PER_SAMPLE=50: []) |
| celltype_to_lineage | **[WARN]** 13,881 cells (5.4%) had a cell type with no lineage in the cell-type registry and were dropped: ['dirt', 'undefined'] — to keep them, give them a lineage in celltype_registry.csv (see registry.py); the reason each is currently excluded is recorded in that file's `notes` column |
| normalise_coords | per-sample z-score applied to X, Y (matches spatial_positional_encoding/src/preprocess.py Step 5) |
| merge_metadata | 140 -> 140 samples kept a metadata row (0 sample(s) had no metadata row and were dropped); survival ABSENT (survival-less cohort) |
| export | 140 sample parquet files -> D:\Desktop\FYDP\FYDP final works\files\data_preprocessing\datasets\CRC_doublets\processed\samples, manifest -> D:\Desktop\FYDP\FYDP final works\files\data_preprocessing\datasets\CRC_doublets\processed\manifest.parquet |
| final_cohort | 244,504 cells, 140 samples, 35 patients, no survival (feature-verification cohort only) |