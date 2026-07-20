# Processing report — Ferguson

What changed, in order, from the raw input files to the processed cohort now sitting in `processed/`. Re-run `python run_ingest.py --dataset <name>` any time to regenerate this after editing `adapter_config.py`.

| Step | Detail |
|---|---|
| load | locations=155,913 rows, expression=155,913 rows, metadata=44 rows, 34 marker columns detected |
| dedup:locations | 155,913 -> 155,913 rows (0 duplicate (sample, cell) pairs dropped) |
| dedup:expression | 155,913 -> 155,913 rows (0 duplicate (sample, cell) pairs dropped) |
| drop_missing:coords | 155,913 -> 155,913 rows (0 dropped for missing X/Y) |
| drop_missing:markers | 155,913 -> 155,913 rows (0 dropped for missing marker values) |
| normalise_markers | arcsinh(x / 1.0) applied to 34 marker columns |
| merge_locations_expression | 155,913 location rows -> 155,913 merged rows (0 location rows had no matching expression row and were dropped) |
| filter_small_samples | 44 -> 44 samples (dropped 0 sample(s) below MIN_CELLS_PER_SAMPLE=50: []) |
| celltype_to_lineage | all 155,913 cells mapped to a lineage: {'immune': 84603, 'stromal': 49231, 'tumour': 22079} |
| normalise_coords | per-sample z-score applied to X, Y (matches spatial_positional_encoding/src/preprocess.py Step 5) |
| merge_metadata | 44 -> 44 samples kept a metadata row (0 sample(s) had no metadata row and were dropped); survival present |
| export | 44 sample parquet files -> C:\Users\User\Downloads\files-main\files-main\data_preprocessing\datasets\Ferguson\processed\samples, manifest -> C:\Users\User\Downloads\files-main\files-main\data_preprocessing\datasets\Ferguson\processed\manifest.parquet |
| final_cohort | 155,913 cells, 44 samples, 17 patients, 0 events (nan% event rate) |