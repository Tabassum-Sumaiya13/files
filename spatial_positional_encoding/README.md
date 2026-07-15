# Spatial Positional Encoding Pipeline

Converts per-cell **(X, Y)** coordinates + **39 protein markers** from CODEX spatial
proteomics data into graph-based **Laplacian positional encodings (PE)**, then
evaluates whether PE features improve patient **survival prediction** over the
spatsurv baseline (celltype proportions → RSF, C-index ≈ 0.704).

## Project Structure

```
spatial_positional_encoding/
├── config.py                  # All paths, column names, hyperparameters
├── run_pipeline.py            # CLI: preprocess → encode (Steps 1-2)
├── run_survival.py            # CLI: survival RSF evaluation (Step 3)
├── requirements.txt           # Python dependencies
│
├── src/
│   ├── __init__.py
│   ├── preprocess.py          # Load raw CSVs, clean, merge, normalise, export
│   ├── graph_construction.py  # k-NN spatial graph from cell coordinates
│   ├── laplacian_encoding.py  # Laplacian eigenvector PE computation
│   ├── feature_assembly.py    # Combine PE + markers into feature vectors
│   └── positional_encoder.py  # Per-sample encoding orchestrator
│
├── data/
│   ├── raw/                   # INPUT — do not modify
│   │   ├── cell_locations_and_labels.csv
│   │   ├── labeled_arcsinh_norm_data.csv
│   │   ├── marker_names.csv
│   │   ├── sample_metadata.csv
│   │   └── qc_acq_ids_labeled.csv
│   ├── processed/             # Step 1 output
│   │   ├── samples/           # One parquet per sample
│   │   ├── manifest.parquet
│   │   └── qc_reports/
│   └── encodings/             # Step 2 output
│       └── encoding_*.parquet
│
└── notebooks/                 # Validation notebooks
```

## Setup

```bash
pip install -r requirements.txt
```

Verify your raw data files exist in `data/raw/`.

## Running — 3 Steps

### Step 1: Preprocess (clean, merge, normalise, export per-sample files)

```bash
python run_pipeline.py --phase preprocess
```

This reads the raw CSVs, filters to 307 QC-passing samples, merges locations
with expression data, normalises coordinates per-sample, and writes one parquet
per sample to `data/processed/samples/`.

### Step 2: Encode (build graphs → compute Laplacian PE)

```bash
python run_pipeline.py --phase encode
```

For each sample: builds a k-NN spatial graph → computes graph Laplacian →
extracts the 8 smallest non-trivial eigenvectors as positional encodings →
concatenates with protein markers → saves to `data/encodings/`.

### Run Steps 1+2 together:

```bash
python run_pipeline.py --phase all
```

### Step 3: Evaluate survival prediction

```bash
python run_survival.py
```

Compares PE features against the spatsurv celltype-proportion baseline using
RandomSurvivalForest with GroupKFold cross-validation (by patient).

## Quick test (debug mode)

Process only 5 samples to verify everything works:

```bash
python run_pipeline.py --phase all --debug 5
```

## Key Parameters (edit in config.py)

| Parameter | Default | Description |
|---|---|---|
| `k_neighbors` | 10 | k-NN graph: edges per cell |
| `k_pe` | 8 | Laplacian eigenvector dimensions |
| `min_cells_per_sample` | 50 | Drop samples smaller than this |
| `normalize_coords` | True | Per-sample z-score of X, Y |
| `use_qc_sample_list` | True | Use only 307 QC-passing samples |
| `n_workers` | 4 | Parallel workers for encoding |
| `rsf_n_estimators` | 100 | RSF trees (matches spatsurv) |
| `rsf_n_splits` | 10 | CV folds (GroupKFold by patient) |

## Feature Vector Per Cell

```
[ pe_0, pe_1, ..., pe_7 | CD31, CD57, CD4, ..., CD3e ]
  ←   8 Laplacian PE  →   ←     39 protein markers    →
  Total: 47 dimensions
```
