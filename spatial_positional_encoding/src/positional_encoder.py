"""
src/positional_encoder.py — Orchestrates the per-sample encoding pipeline:
    load sample parquet -> build graph -> Laplacian PE -> assemble features -> save

Also provides the batch runner with parallel processing.
"""

import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, List
import warnings
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor, as_completed
import traceback

from src.graph_construction import build_spatial_graph
from src.laplacian_encoding import compute_laplacian_encoding
from src.feature_assembly import create_output_dataframe


def process_single_sample(
    sample_file: str,
    output_dir: str,
    marker_cols: List[str],
    k_neighbors: int = 10,
    k_pe: int = 8,
    symmetrization: str = "union",
    keep_largest_component: bool = True,
    sigma_offset: float = 1e-8,
) -> Dict:
    """
    Process one sample: graph -> Laplacian -> PE -> save parquet.

    Returns a result dict with status, stats, and output path.
    """
    try:
        df = pd.read_parquet(sample_file)
        sample_id = df['acquisition_id'].iloc[0]

        print(f"\n  --- Processing: {sample_id} ({len(df)} cells) ---")

        coords = df[['X', 'Y']].values
        markers = df[marker_cols].values

        print(f"    Coords shape: {coords.shape}, dtype: {coords.dtype}")
        print(f"    Markers shape: {markers.shape}, dtype: {markers.dtype}")
        print(f"    Coords range: X=[{coords[:, 0].min():.2f}, {coords[:, 0].max():.2f}], "
              f"Y=[{coords[:, 1].min():.2f}, {coords[:, 1].max():.2f}]")

        # Step A: Build graph
        adj, graph_stats = build_spatial_graph(
            coords, k=k_neighbors, symmetrization=symmetrization
        )

        # Step B: Compute Laplacian PE
        pe_vectors, valid_mask, lap_stats = compute_laplacian_encoding(
            adj, k_pe=k_pe, sigma_offset=sigma_offset,
            keep_largest_component=keep_largest_component,
        )

        # Check for encoding failure
        if lap_stats.get('error') is not None:
            print(f"    [FAILED] {sample_id}: {lap_stats['error']}")
            return {
                'sample_id': sample_id,
                'status': 'failed',
                'error': lap_stats['error'],
                'n_cells': len(df),
            }

        if len(pe_vectors) != len(df):
            raise ValueError(f"PE length {len(pe_vectors)} != data length {len(df)}")

        # Step C: Assemble features + output DataFrame
        output_df, metadata = create_output_dataframe(
            sample_id=sample_id,
            cell_ids=df['cell_id'].values,
            cluster_labels=df['cluster_label'].values,
            pe_vectors=pe_vectors,
            markers=markers,
            marker_names=marker_cols,
            valid_mask=valid_mask,
        )

        n_valid = int(output_df['is_valid_pe'].sum())
        n_invalid = len(df) - n_valid

        # Step D: Save
        output_file = Path(output_dir) / f"encoding_{sample_id}.parquet"
        output_df.to_parquet(output_file, compression='snappy', index=False)

        print(f"    [OK] {sample_id}: {n_valid}/{len(df)} valid PE cells, "
              f"saved to {output_file.name}")

        return {
            'sample_id': sample_id,
            'status': 'success',
            'n_cells': len(df),
            'n_valid_pe': n_valid,
            'n_invalid_pe': n_invalid,
            'file_path': str(output_file),
            'graph_stats': graph_stats,
            'laplacian_stats': lap_stats,
            'feature_dims': metadata['total_dims'],
        }

    except Exception as e:
        sample_id = Path(sample_file).stem.replace('sample_', '', 1)
        tb = traceback.format_exc()
        print(f"    [ERROR] {sample_id}: {e}")
        return {
            'sample_id': sample_id,
            'status': 'failed',
            'error': str(e),
            'traceback': tb,
        }


def run_encoding_pipeline(cfg, marker_cols: List[str]) -> pd.DataFrame:
    """
    Run positional encoding on all preprocessed samples.

    Reads samples from cfg.samples_dir, writes encodings to cfg.encoding_dir.
    Returns a results DataFrame (one row per sample).
    """
    print("\n" + "=" * 70)
    print("  POSITIONAL ENCODING PIPELINE")
    print("=" * 70)

    cfg.ensure_dirs()

    sample_files = sorted(cfg.samples_dir.glob("sample_*.parquet"))
    if not sample_files:
        raise FileNotFoundError(f"No sample files in {cfg.samples_dir}")

    print(f"  Found {len(sample_files)} sample files")
    print(f"  Output directory: {cfg.encoding_dir}")
    print(f"  k_neighbors={cfg.k_neighbors}, k_pe={cfg.k_pe}")
    print(f"  Workers: {cfg.n_workers}")
    print(f"  Markers: {len(marker_cols)}")

    results = []

    if cfg.n_workers > 1:
        print(f"\n  Running with {cfg.n_workers} parallel workers...")
        with ProcessPoolExecutor(max_workers=cfg.n_workers) as executor:
            futures = {
                executor.submit(
                    process_single_sample,
                    str(f), str(cfg.encoding_dir), marker_cols,
                    cfg.k_neighbors, cfg.k_pe, cfg.edge_symmetrization,
                    cfg.keep_largest_component, cfg.sigma_offset,
                ): f
                for f in sample_files
            }
            with tqdm(total=len(futures), desc="Encoding samples") as pbar:
                for future in as_completed(futures):
                    result = future.result()
                    results.append(result)
                    pbar.update(1)
                    if result['status'] == 'failed':
                        pbar.set_description(f"Failed: {result.get('sample_id', '?')}")
    else:
        print(f"\n  Running sequentially...")
        for f in sample_files:
            result = process_single_sample(
                str(f), str(cfg.encoding_dir), marker_cols,
                cfg.k_neighbors, cfg.k_pe, cfg.edge_symmetrization,
                cfg.keep_largest_component, cfg.sigma_offset,
            )
            results.append(result)

    # Summary
    results_df = pd.DataFrame(results)
    n_success = (results_df['status'] == 'success').sum()
    n_failed = (results_df['status'] == 'failed').sum()

    print("\n" + "=" * 70)
    print("  ENCODING COMPLETE")
    print("=" * 70)
    print(f"  Success: {n_success}")
    print(f"  Failed:  {n_failed}")

    if n_failed > 0:
        print(f"\n  Failed samples:")
        for _, row in results_df[results_df['status'] == 'failed'].iterrows():
            print(f"    {row['sample_id']}: {row.get('error', 'unknown')}")

    # Warn about samples that succeeded but dropped cells
    if 'n_invalid_pe' in results_df.columns:
        partial = results_df[
            (results_df['status'] == 'success') &
            (results_df['n_invalid_pe'] > 0)
        ]
        if len(partial) > 0:
            print(f"\n  [WARNING] {len(partial)} samples dropped some cells (disconnected):")
            for _, row in partial.iterrows():
                print(f"    {row['sample_id']}: {row['n_invalid_pe']}/{row['n_cells']} invalid")

    # Save report
    report_path = cfg.encoding_dir / "encoding_report.parquet"
    results_df.to_parquet(report_path, compression='snappy', index=False)
    print(f"\n  Report saved: {report_path}")

    return results_df


def load_all_encodings(encoding_dir: Path) -> pd.DataFrame:
    """Load all per-sample encoding parquets into one DataFrame."""
    files = [
        f for f in sorted(encoding_dir.glob("encoding_*.parquet"))
        if f.name != "encoding_report.parquet"
    ]
    if not files:
        raise FileNotFoundError(f"No encoding files in {encoding_dir}")

    print(f"\n  Loading {len(files)} encoding files...")
    dfs = [pd.read_parquet(f) for f in files]
    combined = pd.concat(dfs, ignore_index=True)
    print(f"  Combined: {combined.shape}, samples: {combined['acquisition_id'].nunique()}")
    return combined