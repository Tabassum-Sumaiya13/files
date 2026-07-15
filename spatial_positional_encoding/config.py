"""
config.py — Single source of truth for all paths, column names, and hyperparameters.

Every script in this project imports from here. Edit THIS file to change
paths or tune parameters — never hardcode them elsewhere.
"""

from pathlib import Path
from dataclasses import dataclass, field
from typing import Optional, List

# ---------------------------------------------------------------------------
# Project root: the directory that contains this config.py file
# ---------------------------------------------------------------------------
PROJECT_ROOT = Path(__file__).parent.resolve()


@dataclass
class Config:
    """All configuration for preprocessing, encoding, and survival evaluation."""

    # ------------------------------------------------------------------
    # Paths  (relative to PROJECT_ROOT, resolved in __post_init__)
    # ------------------------------------------------------------------
    raw_data_dir: str = "data/raw"
    processed_data_dir: str = "data/processed"
    encoding_output_dir: str = "data/encodings"

    # ------------------------------------------------------------------
    # Raw-file names  (inside raw_data_dir)
    # ------------------------------------------------------------------
    locations_file: str = "cell_locations_and_labels.csv"
    expression_file: str = "labeled_arcsinh_norm_data.csv"
    marker_names_file: str = "marker_names.csv"
    sample_metadata_file: str = "sample_metadata.csv"
    qc_sample_ids_file: str = "qc_acq_ids_labeled.csv"  # 307 QC-passing IDs

    # ------------------------------------------------------------------
    # Column names in the RAW CSVs  (exactly as they appear on disk)
    # ------------------------------------------------------------------
    # cell_locations_and_labels.csv
    LOC_ACQ_COL: str = "ACQUISITION_ID"
    LOC_CELL_COL: str = "CELL_ID"
    LOC_X_COL: str = "X"
    LOC_Y_COL: str = "Y"
    LOC_CLUSTER_ID_COL: str = "CLUSTER_ID"
    LOC_CLUSTER_LABEL_COL: str = "CLUSTER_LABEL"

    # labeled_arcsinh_norm_data.csv
    EXPR_SAMPLE_COL: str = "sample_id"
    EXPR_CELL_COL: str = "cell_id"
    EXPR_CLUSTER_COL: str = "cluster"
    EXPR_CLUSTER_LABEL_COL: str = "cluster_label"

    # sample_metadata.csv
    META_ACQ_COL: str = "acquisition_id"
    META_PATIENT_COL: str = "patient_id"
    META_COVERSLIP_COL: str = "coverslip_label"
    META_SURVIVAL_TIME_COL: str = "survival_day"
    META_SURVIVAL_STATUS_COL: str = "survival_status"

    # ------------------------------------------------------------------
    # Internal standardised column names  (used after preprocessing)
    # ------------------------------------------------------------------
    ACQ_COL: str = "acquisition_id"
    CELL_COL: str = "cell_id"
    X_COL: str = "X"
    Y_COL: str = "Y"
    CLUSTER_ID_COL: str = "cluster_id"
    CLUSTER_LABEL_COL: str = "cluster_label"

    # ------------------------------------------------------------------
    # Preprocessing parameters
    # ------------------------------------------------------------------
    min_cells_per_sample: int = 50        # Drop samples smaller than this
    normalize_coords: bool = True         # Per-sample z-score X, Y
    use_qc_sample_list: bool = True       # Only keep 307 QC-passing samples
    debug_n_samples: Optional[int] = None # Limit to N samples for fast iteration

    # ------------------------------------------------------------------
    # Graph construction
    # ------------------------------------------------------------------
    k_neighbors: int = 10                 # k in k-NN spatial graph
    edge_symmetrization: str = "union"    # "union" or "mutual"

    # ------------------------------------------------------------------
    # Positional encoding
    # ------------------------------------------------------------------
    k_pe: int = 8                         # Number of Laplacian eigenvector dims
    keep_largest_component: bool = True    # Discard disconnected subgraphs
    sigma_offset: float = 1e-8            # Shift-invert stabiliser for eigsh

    # ------------------------------------------------------------------
    # Parallel processing
    # ------------------------------------------------------------------
    n_workers: int = 4

    # ------------------------------------------------------------------
    # Survival evaluation
    # ------------------------------------------------------------------
    rsf_n_estimators: int = 100
    rsf_n_splits: int = 10
    rsf_random_state: int = 1029          # Same as spatsurv baseline
    pooling_method: str = "mean_std"      # How to collapse cell->sample features

    # ------------------------------------------------------------------
    # Resolved Path objects  (set in __post_init__)
    # ------------------------------------------------------------------
    _raw_dir: Path = field(init=False, repr=False)
    _processed_dir: Path = field(init=False, repr=False)
    _encoding_dir: Path = field(init=False, repr=False)
    _samples_dir: Path = field(init=False, repr=False)
    _qc_dir: Path = field(init=False, repr=False)

    def __post_init__(self):
        self._raw_dir = PROJECT_ROOT / self.raw_data_dir
        self._processed_dir = PROJECT_ROOT / self.processed_data_dir
        self._encoding_dir = PROJECT_ROOT / self.encoding_output_dir
        self._samples_dir = self._processed_dir / "samples"
        self._qc_dir = self._processed_dir / "qc_reports"

    # --- convenience properties ---

    @property
    def raw_dir(self) -> Path:
        return self._raw_dir

    @property
    def processed_dir(self) -> Path:
        return self._processed_dir

    @property
    def encoding_dir(self) -> Path:
        return self._encoding_dir

    @property
    def samples_dir(self) -> Path:
        return self._samples_dir

    @property
    def qc_dir(self) -> Path:
        return self._qc_dir

    @property
    def manifest_path(self) -> Path:
        return self._processed_dir / "manifest.parquet"

    # Raw file paths
    @property
    def locations_path(self) -> Path:
        return self._raw_dir / self.locations_file

    @property
    def expression_path(self) -> Path:
        return self._raw_dir / self.expression_file

    @property
    def marker_names_path(self) -> Path:
        return self._raw_dir / self.marker_names_file

    @property
    def metadata_path(self) -> Path:
        return self._raw_dir / self.sample_metadata_file

    @property
    def qc_ids_path(self) -> Path:
        return self._raw_dir / self.qc_sample_ids_file

    def ensure_dirs(self):
        """Create all output directories."""
        for d in [self._raw_dir, self._processed_dir, self._encoding_dir,
                  self._samples_dir, self._qc_dir]:
            d.mkdir(parents=True, exist_ok=True)

    def print_summary(self):
        """Print the active configuration for debugging."""
        print("\n" + "=" * 70)
        print("  CONFIGURATION SUMMARY")
        print("=" * 70)
        print(f"  Project root        : {PROJECT_ROOT}")
        print(f"  Raw data dir        : {self._raw_dir}")
        print(f"  Processed data dir  : {self._processed_dir}")
        print(f"  Encoding output dir : {self._encoding_dir}")
        print(f"  Samples dir         : {self._samples_dir}")
        print()
        print(f"  use_qc_sample_list  : {self.use_qc_sample_list}")
        print(f"  debug_n_samples     : {self.debug_n_samples}")
        print(f"  min_cells_per_sample: {self.min_cells_per_sample}")
        print(f"  normalize_coords    : {self.normalize_coords}")
        print()
        print(f"  k_neighbors (graph) : {self.k_neighbors}")
        print(f"  k_pe (eigenvectors) : {self.k_pe}")
        print(f"  keep_largest_comp.  : {self.keep_largest_component}")
        print(f"  n_workers           : {self.n_workers}")
        print()
        print(f"  RSF estimators      : {self.rsf_n_estimators}")
        print(f"  RSF CV splits       : {self.rsf_n_splits}")
        print(f"  RSF random_state    : {self.rsf_random_state}")
        print(f"  pooling_method      : {self.pooling_method}")
        print("=" * 70 + "\n")


# ---------------------------------------------------------------------------
# Global singleton — import this everywhere
# ---------------------------------------------------------------------------
_cfg: Optional[Config] = None


def get_config(**overrides) -> Config:
    """Return the global Config, creating it on first call.

    Pass keyword arguments to override defaults:
        cfg = get_config(debug_n_samples=5, k_pe=4)
    """
    global _cfg
    if _cfg is None or overrides:
        _cfg = Config(**overrides)
    return _cfg