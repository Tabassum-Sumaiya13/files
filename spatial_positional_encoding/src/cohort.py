"""
src/cohort.py — Load a cohort that data_preprocessing has already ingested into
the canonical format, so every dataset (UPMC, CRC, …) is read the exact same way.

A canonical cohort lives at:
    data_preprocessing/datasets/<name>/processed/
        samples/sample_<acquisition_id>.parquet   # one per tissue sample
        manifest.parquet                           # per-sample survival + patient
        marker_columns.txt                         # marker column names, one per line

Each sample parquet carries the columns the feature blocks need:
    acquisition_id, cell_id, X, Y, cluster_label, lineage, <markers…>
    (+ patient_id, survival_day, survival_status merged in per row)

`lineage` (immune|tumour|stromal) is the PORTABLE label present in every cohort;
`cluster_label` is the native per-dataset label. Feature builders take the label
column as a parameter so the same code runs on either.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, List, Optional, Tuple

import numpy as np
import pandas as pd

# spatial_positional_encoding/ and data_preprocessing/ are siblings under the repo root.
_REPO_ROOT = Path(__file__).resolve().parents[2]
_DATASETS_DIR = _REPO_ROOT / "data_preprocessing" / "datasets"

LINEAGES = ("immune", "tumour", "stromal")


@dataclass
class Cohort:
    """A single ingested cohort, read from its canonical processed/ folder."""

    name: str
    processed_dir: Path
    manifest: pd.DataFrame
    marker_columns: List[str]

    # --- discovery -------------------------------------------------------
    @property
    def samples_dir(self) -> Path:
        return self.processed_dir / "samples"

    @property
    def sample_ids(self) -> List[str]:
        return self.manifest["acquisition_id"].astype(str).tolist()

    def sample_path(self, sid: str) -> Path:
        return self.samples_dir / f"sample_{sid}.parquet"

    # --- loading ---------------------------------------------------------
    def load_sample(self, sid: str, columns: Optional[List[str]] = None) -> pd.DataFrame:
        return pd.read_parquet(self.sample_path(sid), columns=columns)

    def iter_samples(
        self, columns: Optional[List[str]] = None, limit: Optional[int] = None
    ) -> Iterator[Tuple[str, pd.DataFrame]]:
        ids = self.sample_ids[:limit] if limit else self.sample_ids
        for sid in ids:
            p = self.sample_path(sid)
            if p.exists():
                yield sid, pd.read_parquet(p, columns=columns)

    # --- labels ----------------------------------------------------------
    def has_survival(self) -> bool:
        """True only if the manifest carries usable (non-null, event-bearing) survival."""
        m = self.manifest
        if not {"survival_day", "survival_status"} <= set(m.columns):
            return False
        t = pd.to_numeric(m["survival_day"], errors="coerce")
        e = pd.to_numeric(m["survival_status"], errors="coerce")
        return bool(t.notna().any() and e.isin([0, 1]).any() and e.fillna(0).sum() > 0)

    def native_labels(self) -> List[str]:
        """Sorted set of native cluster_label values across the cohort (for reporting)."""
        seen = set()
        for _, df in self.iter_samples(columns=["cluster_label"]):
            seen.update(df["cluster_label"].astype(str).unique())
        return sorted(seen)


def dataset_dir(name: str) -> Path:
    return _DATASETS_DIR / name


def list_datasets() -> List[str]:
    """Names of every dataset that has been ingested (has a processed/manifest.parquet)."""
    if not _DATASETS_DIR.exists():
        return []
    out = []
    for d in sorted(_DATASETS_DIR.iterdir()):
        if (d / "processed" / "manifest.parquet").exists():
            out.append(d.name)
    return out


def load_cohort(name: str) -> Cohort:
    """Load an ingested cohort by dataset folder name (e.g. 'UPMC', 'CRC')."""
    processed = dataset_dir(name) / "processed"
    manifest_path = processed / "manifest.parquet"
    if not manifest_path.exists():
        raise FileNotFoundError(
            f"No ingested cohort for '{name}' at {manifest_path}.\n"
            f"Run:  cd data_preprocessing && python run_ingest.py --dataset {name}"
        )
    manifest = pd.read_parquet(manifest_path)
    manifest["acquisition_id"] = manifest["acquisition_id"].astype(str)

    markers_txt = processed / "marker_columns.txt"
    marker_columns = (
        [ln.strip() for ln in markers_txt.read_text(encoding="utf-8").splitlines() if ln.strip()]
        if markers_txt.exists()
        else []
    )
    return Cohort(name=name, processed_dir=processed, manifest=manifest, marker_columns=marker_columns)
