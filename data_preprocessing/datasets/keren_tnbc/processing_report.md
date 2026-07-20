# Processing report — keren_tnbc

What changed, in order, from the raw input files to the processed cohort now sitting in `processed/`. Re-run `python run_ingest.py --dataset <name>` any time to regenerate this after editing `adapter_config.py`.

| Step | Detail |
|---|---|
| load | locations=197,678 rows, expression=197,678 rows, metadata=40 rows, 49 marker columns detected |
| dedup:locations | 197,678 -> 197,678 rows (0 duplicate (sample, cell) pairs dropped) |
| dedup:expression | 197,678 -> 197,678 rows (0 duplicate (sample, cell) pairs dropped) |
| drop_missing:coords | 197,678 -> 197,678 rows (0 dropped for missing X/Y) |
| drop_missing:markers | 197,678 -> 197,678 rows (0 dropped for missing marker values) |
| normalise_markers | skipped — adapter_config.APPLY_ARCSINH=False, data already normalised |
| merge_locations_expression | 197,678 location rows -> 197,678 merged rows (0 location rows had no matching expression row and were dropped) |
| filter_small_samples | 40 -> 40 samples (dropped 0 sample(s) below MIN_CELLS_PER_SAMPLE=50: []) |
| celltype_to_lineage | **[WARN]** 1,725 cells (0.9%) had a cell type not in CELLTYPE_MAP and were dropped: ['Unidentified'] — add them to adapter_config.py to keep them |
| normalise_coords | per-sample z-score applied to X, Y (matches spatial_positional_encoding/src/preprocess.py Step 5) |
| merge_metadata | 40 -> 40 samples kept a metadata row (0 sample(s) had no metadata row and were dropped); survival present |
| export | 40 sample parquet files -> C:\Users\MARZIA ISLAM\Downloads\files\data_preprocessing\datasets\keren_tnbc\processed\samples, manifest -> C:\Users\MARZIA ISLAM\Downloads\files\data_preprocessing\datasets\keren_tnbc\processed\manifest.parquet |
| final_cohort | 195,953 cells, 40 samples, 40 patients, 13 events (39.4% event rate) |