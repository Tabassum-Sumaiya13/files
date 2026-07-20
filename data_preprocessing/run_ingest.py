"""
run_ingest.py — Entry point for bringing a new external cohort into this
project's canonical schema.

    python run_ingest.py --dataset schurch_codex
    python run_ingest.py --dataset schurch_codex --validate-only
    python run_ingest.py --dataset schurch_codex --force

Looks for datasets/<name>/adapter_config.py, runs validation against
schema.py, writes datasets/<name>/validation_report.md, and — only if the
dataset is READY (zero FAIL checks) or --force is given — processes it into
datasets/<name>/processed/ and writes datasets/<name>/processing_report.md.

To add a new dataset see datasets/_template/adapter_config.py.
"""
import argparse
import importlib.util
import sys
from pathlib import Path

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))

from validator import validate_dataset   # noqa: E402
from processor import process_dataset    # noqa: E402
from report import ChangeLog             # noqa: E402


def load_config(dataset_name: str):
    path = HERE / "datasets" / dataset_name / "adapter_config.py"
    if not path.exists():
        raise FileNotFoundError(
            f"No adapter config at {path}.\n"
            f"Copy datasets/_template/ to datasets/{dataset_name}/ and fill in "
            f"adapter_config.py first."
        )
    spec = importlib.util.spec_from_file_location(f"adapter_{dataset_name}", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dataset", required=True, help="Folder name under datasets/")
    ap.add_argument("--validate-only", action="store_true", help="Only run validation, do not process")
    ap.add_argument("--force", action="store_true", help="Process even if validation has FAIL checks (not recommended)")
    args = ap.parse_args()

    cfg = load_config(args.dataset)
    dataset_dir = HERE / "datasets" / args.dataset

    report = validate_dataset(cfg)
    report_path = dataset_dir / "validation_report.md"
    report.save(report_path)
    print(f"\n{report.to_markdown()}\n")
    print(f"Validation report written to {report_path}")

    # Machine-readable form of the marker-evidence table, so the downstream
    # sensitivity analysis (run_verify.py --perturb-map) can read which lineage
    # the evidence predicted for each contested cluster instead of parsing markdown.
    if getattr(report, "evidence_table", None) is not None:
        ev_path = dataset_dir / "lineage_evidence.csv"
        report.evidence_table.to_csv(ev_path, index=False)
        print(f"Marker-evidence table written to {ev_path}")

    if args.validate_only:
        return

    if not report.ready and not args.force:
        print(f"\n[STOPPED] {report.n_fail} FAIL check(s) above — fix the raw data or "
              f"adapter_config.py, or re-run with --force to process anyway (not recommended, "
              f"results on a cohort that failed schema checks should be treated with suspicion).")
        return

    print(f"\n{'=' * 70}\n  PROCESSING: {cfg.DATASET_NAME}\n{'=' * 70}")
    log = ChangeLog(cfg.DATASET_NAME)
    out_dir = dataset_dir / "processed"
    manifest = process_dataset(cfg, out_dir, log)
    log_path = dataset_dir / "processing_report.md"
    log.save(log_path)

    print(f"\nProcessing report written to {log_path}")
    print(f"Processed cohort written to {out_dir}")
    print(f"  {len(manifest)} samples, {int(manifest['n_cells'].sum()):,} cells total")


if __name__ == "__main__":
    main()
