# Processing report — UPMC

What changed, in order, from the raw input files to the processed cohort now sitting in `processed/`. Re-run `python run_ingest.py --dataset <name>` any time to regenerate this after editing `adapter_config.py`.

| Step | Detail |
|---|---|
| load | locations=2,061,102 rows, expression=2,061,102 rows, metadata=379 rows, 39 marker columns detected |
| dedup:locations | 2,061,102 -> 2,061,102 rows (0 duplicate (sample, cell) pairs dropped) |
| dedup:expression | 2,061,102 -> 2,061,102 rows (0 duplicate (sample, cell) pairs dropped) |
| drop_missing:coords | 2,061,102 -> 2,061,102 rows (0 dropped for missing X/Y) |
| drop_missing:markers | 2,061,102 -> 2,061,102 rows (0 dropped for missing marker values) |
| normalise_markers | skipped — adapter_config.APPLY_ARCSINH=False, data already normalised |
| merge_locations_expression | 2,061,102 location rows -> 2,061,102 merged rows (0 location rows had no matching expression row and were dropped) |
| filter_small_samples | 308 -> 308 samples (dropped 0 sample(s) below MIN_CELLS_PER_SAMPLE=50: []) |
| celltype_to_lineage | all 2,061,102 cells mapped to a lineage: {'tumour': 910312, 'immune': 908224, 'stromal': 242566} |
| normalise_coords | per-sample z-score applied to X, Y (matches spatial_positional_encoding/src/preprocess.py Step 5) |
| merge_metadata | 308 -> 308 samples kept a metadata row (0 sample(s) had no metadata row and were dropped); survival present |
| export | 308 sample parquet files -> D:\Desktop\FYDP\FYDP final works\files\data_preprocessing\datasets\UPMC\processed\samples, manifest -> D:\Desktop\FYDP\FYDP final works\files\data_preprocessing\datasets\UPMC\processed\manifest.parquet |
| final_cohort | 2,061,102 cells, 308 samples, 81 patients, 103 events (33.4% event rate) |