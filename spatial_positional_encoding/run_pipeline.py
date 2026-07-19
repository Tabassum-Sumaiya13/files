#!/usr/bin/env python3
"""
run_pipeline.py — Orchestrator for the current (post-PE) pipeline.

The pipeline is now two phases, and both are dataset-agnostic — UPMC, CRC and
any future cohort flow through the identical path:

  Phase 1 — INGEST  (data_preprocessing/)
      Bring a raw cohort into the canonical schema (per-sample parquet with
      X/Y, cluster_label, lineage, markers + a manifest). Survival is OPTIONAL:
      a cohort with none still ingests and is fully feature-verifiable.
          cd ../data_preprocessing
          python run_ingest.py --dataset <NAME>

  Phase 2 — VERIFY  (run_verify.py)
      Dataset-agnostic feature verification — does the spatial feature block
      carry real, reproducible signal beyond the composition baseline? Needs no
      survival. See run_verify.py for the metric definitions.
          python run_verify.py --dataset <NAME> --taxonomy lineage

  Survival, when a cohort has it, is a SEPARATE per-dataset downstream check
  (RandomSurvivalForest + patient-grouped C-index) — not part of this common
  path. See run_survival.py.

Positional encoding (Laplacian PE) has been removed from the pipeline; the
spatial feature block is the abundance-corrected enrichment readout
(src/spatial_features.py).

This script is a convenience dispatcher: it checks a dataset is ingested and
runs Phase 2 on it.
"""
import argparse
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from src.cohort import list_datasets, load_cohort


def main():
    parser = argparse.ArgumentParser(
        description="Post-PE pipeline orchestrator (ingest -> verify)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--dataset", help="Ingested dataset folder name (UPMC, CRC, …)")
    parser.add_argument("--taxonomy", default="lineage", choices=["lineage", "native"],
                        help="lineage = portable 3-way (default); native = per-dataset resolution")
    parser.add_argument("--list", action="store_true", help="List ingested datasets and exit")
    parser.add_argument("--label-col", default=None,
                        help="optional manifest column for the separability check")
    args, extra = parser.parse_known_args()

    available = list_datasets()
    if args.list or not args.dataset:
        print("Ingested datasets (data_preprocessing/datasets/<name>/processed/):")
        for d in available:
            co = load_cohort(d)
            print(f"  - {d:12s} {len(co.sample_ids):>4d} samples  "
                  f"survival={'yes' if co.has_survival() else 'no'}")
        if not available:
            print("  (none yet — run:  cd ../data_preprocessing && python run_ingest.py --dataset <NAME>)")
        if not args.dataset:
            print("\nThen:  python run_pipeline.py --dataset <NAME>")
            return

    if args.dataset not in available:
        print(f"[ERROR] '{args.dataset}' is not ingested. Available: {available or '(none)'}")
        print(f"        Ingest it first: cd ../data_preprocessing && "
              f"python run_ingest.py --dataset {args.dataset}")
        sys.exit(1)

    # Phase 2 — delegate to run_verify.py (single source of truth for the matrix)
    cmd = [sys.executable, str(Path(__file__).parent / "run_verify.py"),
           "--dataset", args.dataset, "--taxonomy", args.taxonomy]
    if args.label_col:
        cmd += ["--label-col", args.label_col]
    cmd += extra
    print(f"\n[run_pipeline] Phase 2 verify -> {' '.join(cmd[1:])}\n")
    raise SystemExit(subprocess.call(cmd))


if __name__ == "__main__":
    main()
