#!/usr/bin/env python3
"""
run_pipeline.py — Single entry point for the full pipeline.

Usage:
    python run_pipeline.py --phase preprocess     # Step 1: clean + export per-sample files
    python run_pipeline.py --phase encode         # Step 2: graph -> Laplacian PE
    python run_pipeline.py --phase all            # Step 1 + 2 together

    python run_pipeline.py --phase preprocess --debug 5   # Fast iteration: 5 samples only
"""

import sys
import argparse
from pathlib import Path

# Ensure project root is on sys.path
sys.path.insert(0, str(Path(__file__).parent))

from config import get_config
from src.preprocess import run_preprocessing, load_marker_names
from src.positional_encoder import run_encoding_pipeline


def main():
    parser = argparse.ArgumentParser(description="Spatial PE Pipeline")
    parser.add_argument(
        '--phase', required=True,
        choices=['preprocess', 'encode', 'all'],
        help="Which phase to run"
    )
    parser.add_argument(
        '--debug', type=int, default=None,
        help="Limit to N samples for fast iteration"
    )
    parser.add_argument(
        '--workers', type=int, default=None,
        help="Number of parallel workers for encoding"
    )
    parser.add_argument(
        '--k-neighbors', type=int, default=None,
        help="k for k-NN graph (default: 10)"
    )
    parser.add_argument(
        '--k-pe', type=int, default=None,
        help="Number of PE dimensions (default: 8)"
    )
    args = parser.parse_args()

    # Build config with any overrides
    overrides = {}
    if args.debug is not None:
        overrides['debug_n_samples'] = args.debug
    if args.workers is not None:
        overrides['n_workers'] = args.workers
    if args.k_neighbors is not None:
        overrides['k_neighbors'] = args.k_neighbors
    if args.k_pe is not None:
        overrides['k_pe'] = args.k_pe

    cfg = get_config(**overrides)
    cfg.print_summary()

    # ----------------------------------------------------------------
    if args.phase in ('preprocess', 'all'):
        normalized, manifest, marker_cols = run_preprocessing(cfg)

    # ----------------------------------------------------------------
    if args.phase in ('encode', 'all'):
        # If we didn't just preprocess, load marker names from file
        if args.phase == 'encode':
            marker_cols = load_marker_names(cfg)

        results_df = run_encoding_pipeline(cfg, marker_cols)

        # Final summary
        n_ok = (results_df['status'] == 'success').sum()
        print("\n" + "=" * 70)
        print("  PIPELINE COMPLETE")
        print("=" * 70)
        print(f"  Encoded samples: {n_ok}")
        if 'n_valid_pe' in results_df.columns:
            total_valid = results_df.loc[results_df['status'] == 'success', 'n_valid_pe'].sum()
            total_cells = results_df.loc[results_df['status'] == 'success', 'n_cells'].sum()
            print(f"  Total cells: {total_cells:,}")
            print(f"  Cells with valid PE: {total_valid:,} "
                  f"({total_valid / max(total_cells, 1):.1%})")
        if 'feature_dims' in results_df.columns:
            dims = results_df.loc[results_df['status'] == 'success', 'feature_dims'].iloc[0]
            print(f"  Feature dims per cell: {dims} "
                  f"({cfg.k_pe} PE + {dims - cfg.k_pe} markers)")
        print(f"\n  Output: {cfg.encoding_dir}")
        print(f"  Next step: python run_survival.py")


if __name__ == "__main__":
    main()
