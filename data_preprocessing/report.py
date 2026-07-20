"""
report.py — Shared change-log + markdown rendering used by validator.py and
processor.py, so every run leaves behind a human-readable audit trail of
what was checked and what was changed.
"""
from dataclasses import dataclass
from pathlib import Path
from typing import List


@dataclass
class Check:
    name: str
    status: str  # "PASS" | "WARN" | "FAIL"
    detail: str


class ValidationReport:
    """Accumulates the results of validator.py's schema checks."""

    def __init__(self, dataset_name: str):
        self.dataset_name = dataset_name
        self.checks: List[Check] = []
        # Per-cluster marker-evidence table from lineage_evidence.evaluate(),
        # attached by validator.py when the falsifiability gate could run.
        # Rendered in full below the check table so a contradiction can be read
        # off directly rather than inferred from a one-line summary.
        self.evidence_table = None
        # {native_label: justification} for contradictions accepted in the registry.
        self.overrides: dict = {}
        # Provenance of the cell-type mapping used for this run.
        self.provenance: List[str] = []

    def add(self, name: str, status: str, detail: str):
        assert status in ("PASS", "WARN", "FAIL")
        self.checks.append(Check(name, status, detail))
        print(f"  [{status}] {name}: {detail}")

    @property
    def n_fail(self) -> int:
        return sum(c.status == "FAIL" for c in self.checks)

    @property
    def n_warn(self) -> int:
        return sum(c.status == "WARN" for c in self.checks)

    @property
    def ready(self) -> bool:
        """READY = no FAILs. WARNs are fine to proceed with (they usually
        just mean some downstream feature block can't be reproduced)."""
        return self.n_fail == 0 and len(self.checks) > 0

    def to_markdown(self) -> str:
        n_pass = len(self.checks) - self.n_fail - self.n_warn
        lines = [
            f"# Validation report — {self.dataset_name}",
            "",
            f"**Result: {'READY' if self.ready else 'NOT READY'}** "
            f"— {n_pass} pass, {self.n_warn} warn, {self.n_fail} fail",
            "",
            "A FAIL means the schema requirement in schema.py is not met and "
            "processing will refuse to run without `--force`. A WARN means "
            "processing can proceed but some downstream feature block "
            "(a specific marker, a cell-type mapping, sample size) will be "
            "degraded or unavailable for this cohort.",
            "",
            "| Check | Status | Detail |",
            "|---|---|---|",
        ]
        for c in self.checks:
            detail = c.detail.replace("|", "\\|")
            lines.append(f"| {c.name} | {c.status} | {detail} |")

        if self.provenance:
            lines += ["", "## Cell-type mapping provenance", ""] + \
                     [f"- {p}" for p in self.provenance]

        if self.evidence_table is not None and len(self.evidence_table):
            lines += [
                "",
                "## Marker evidence per native cluster",
                "",
                "Each cluster's mean marker profile, z-scored across clusters, scored "
                "against the three lineage core panels (`schema.LINEAGE_MARKER_PANELS`). "
                "`predicted` is the argmax; `margin` is top minus runner-up. "
                "**AGREE** = evidence supports the registry. "
                "**AMBIGUOUS** = margin too small to decide either way. "
                "**CONTRADICTED** = evidence favours a different lineage by a clear margin — "
                "the registry row must be justified in its `notes` or changed, and carried "
                "into `run_verify.py --perturb-map`.",
                "",
                "| " + " | ".join(self.evidence_table.columns) + " |",
                "|" + "---|" * len(self.evidence_table.columns),
            ]
            for row in self.evidence_table.itertuples(index=False):
                lines.append("| " + " | ".join(str(v) for v in row) + " |")

        if self.overrides:
            lines += [
                "",
                "## Accepted contradictions",
                "",
                "Contradictions kept deliberately, each with a written reason recorded in "
                "`celltype_registry.csv` (`evidence_override`). These are downgraded from "
                "blocking to reported — never hidden. They are still counted above and are "
                "still carried into the perturbation analysis, which is what decides whether "
                "any conclusion actually depends on them.",
                "",
            ]
            for label, reason in self.overrides.items():
                lines.append(f"- **{label}** — {reason}")

        return "\n".join(lines)

    def save(self, path: Path):
        Path(path).parent.mkdir(parents=True, exist_ok=True)
        Path(path).write_text(self.to_markdown(), encoding="utf-8")


@dataclass
class LogStep:
    name: str
    detail: str
    level: str = "info"  # "info" | "warn"


class ChangeLog:
    """Accumulates what processor.py did, in order, so the transformation
    from raw input to processed cohort is fully auditable afterwards."""

    def __init__(self, dataset_name: str):
        self.dataset_name = dataset_name
        self.steps: List[LogStep] = []

    def step(self, name: str, detail: str, level: str = "info"):
        self.steps.append(LogStep(name, detail, level))
        tag = "[WARN]" if level == "warn" else "[OK]"
        print(f"  {tag} {name}: {detail}")

    def to_markdown(self) -> str:
        lines = [
            f"# Processing report — {self.dataset_name}",
            "",
            "What changed, in order, from the raw input files to the processed "
            "cohort now sitting in `processed/`. Re-run "
            "`python run_ingest.py --dataset <name>` any time to regenerate "
            "this after editing `adapter_config.py`.",
            "",
            "| Step | Detail |",
            "|---|---|",
        ]
        for s in self.steps:
            detail = s.detail.replace("|", "\\|")
            if s.level == "warn":
                detail = f"**[WARN]** {detail}"
            lines.append(f"| {s.name} | {detail} |")
        return "\n".join(lines)

    def save(self, path: Path):
        Path(path).parent.mkdir(parents=True, exist_ok=True)
        Path(path).write_text(self.to_markdown(), encoding="utf-8")
