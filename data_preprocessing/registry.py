"""
registry.py — The cell-type -> lineage mapping, externalised.

WHY THIS EXISTS
---------------
The mapping used to be a hardcoded Python dict inside each dataset's
adapter_config.py, duplicated in five places across the repo in two
incompatible spellings ('tumour' vs 'TUMOR'), with no citation, no version, and
nothing anywhere asserting the copies agreed. The validator checked only that
every native label was a KEY in the dict — so {"tumor cells": "immune"} passed
every check it had.

The mapping cannot be made automatic: no algorithm resolves a native label like
"tumor cells / immune cells" into a lineage. What CAN be fixed is the status of
that human judgement. This module makes it:

  externally grounded  lineage is a FUNCTION of a Cell Ontology term, looked up
                       in schema.CL_LINEAGE_ANCHOR — not a free choice per row.
                       Two rows with the same cl_term_id can never disagree.
  cited                every row records the publication it came from.
  versioned            REGISTRY_VERSION + a content hash go into every report,
                       so a result is traceable to the exact mapping used.
  falsifiable          lineage_evidence.py scores each cluster against marker
                       data and can CONTRADICT a row. See validator check #11.

WHAT THIS MODULE DOES NOT DO
----------------------------
It does not traverse the Cell Ontology. A real traversal needs the OBO file and
would derive lineage from `is_a` ancestry. Instead schema.CL_LINEAGE_ANCHOR is a
small curated table of exactly the CL terms this project uses. That is weaker: it
catches inconsistency and gives every row a traceable term id, but it does not
prove the term itself is the right one. Every row therefore carries
`verified=no` until a human confirms the term against the ontology (see
docs/CELLTYPE_MAPPING.md for the one-time procedure).

THE FILE
--------
celltype_registry.csv, one row per (dataset, native_label):

    dataset,native_label,lineage,cl_term_id,cl_label,source,verified,notes

An EMPTY `lineage` means the type is intentionally excluded (artifact, unassigned
cluster, or ambiguous doublet); `notes` must then say why. Exclusion used to be a
code comment — now it is data, and the reason ships with the result.
"""
from __future__ import annotations

import hashlib
from pathlib import Path
from typing import Dict, List, Optional

import pandas as pd

REGISTRY_VERSION = "1.0.0"
REGISTRY_PATH = Path(__file__).parent / "celltype_registry.csv"

REQUIRED_COLUMNS = [
    "dataset", "native_label", "lineage", "cl_term_id", "cl_label",
    "source", "verified", "notes", "evidence_override",
]


def registry_fingerprint(path: Path = REGISTRY_PATH) -> str:
    """Short content hash — pins a report to the exact mapping that produced it."""
    return hashlib.sha256(path.read_bytes()).hexdigest()[:12]


def load_registry(path: Path = REGISTRY_PATH) -> pd.DataFrame:
    """Read the registry, normalising blanks to empty strings."""
    if not path.exists():
        raise FileNotFoundError(f"Cell-type registry not found at {path}")
    df = pd.read_csv(path, dtype=str, keep_default_na=False).fillna("")
    missing = [c for c in REQUIRED_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(f"{path.name} is missing required column(s): {missing}")
    for c in REQUIRED_COLUMNS:
        df[c] = df[c].astype(str).str.strip()
    return df


def dataset_rows(dataset: str, path: Path = REGISTRY_PATH) -> pd.DataFrame:
    return load_registry(path)[lambda d: d["dataset"] == dataset].reset_index(drop=True)


def celltype_map(dataset: str, path: Path = REGISTRY_PATH) -> Dict[str, str]:
    """The {native_label: lineage} dict the processor consumes.

    Excluded types (blank lineage) are omitted, which is exactly how the
    processor already treats an unmapped label — so this is a drop-in
    replacement for the old hardcoded CELLTYPE_MAP with the reason for each
    exclusion now recorded in the registry instead of a code comment.
    """
    rows = dataset_rows(dataset, path)
    return {r.native_label: r.lineage for r in rows.itertuples() if r.lineage}


def excluded_types(dataset: str, path: Path = REGISTRY_PATH) -> Dict[str, str]:
    """{native_label: reason} for every type deliberately left unmapped."""
    rows = dataset_rows(dataset, path)
    return {r.native_label: (r.notes or "(no reason recorded)")
            for r in rows.itertuples() if not r.lineage}


def evidence_overrides(dataset: str, path: Path = REGISTRY_PATH) -> Dict[str, str]:
    """{native_label: justification} for contradictions that have been ACCEPTED.

    The escape hatch for the falsifiability gate. `--force` is the wrong one: it
    is global, silent, and suppresses every check at once, so using it once
    teaches you to use it always. An override here is the opposite in each
    respect — it is per cluster, it must carry a written reason, that reason is
    printed in the validation report, and it changes nothing about any other row.

    An overridden contradiction is still COUNTED, still LISTED, and still shown
    in the per-cluster table. It is downgraded from blocking to reported, which
    is the difference between "we did not notice" and "we noticed, and here is
    why we kept it".
    """
    rows = dataset_rows(dataset, path)
    return {r.native_label: r.evidence_override
            for r in rows.itertuples() if r.evidence_override}


def validate_registry(dataset: str, native_labels: Optional[List[str]] = None,
                      path: Path = REGISTRY_PATH) -> List[tuple]:
    """Structural checks on the registry itself.

    Returns a list of (check_name, status, detail) for validator.py to record.
    Checks the things a hardcoded dict could never be checked for:
      - the dataset has rows at all
      - no duplicate native_label
      - every lineage value is in the canonical vocabulary
      - lineage agrees with the CL term's anchor (the grounding check)
      - mapped rows carry a cl_term_id and a source
      - CL terms are marked verified
      - registry covers exactly the labels present in the raw data
    """
    import schema  # local import: schema imports nothing from here

    out: List[tuple] = []
    rows = dataset_rows(dataset, path)

    if rows.empty:
        out.append(("registry:rows", "FAIL",
                    f"no rows for dataset '{dataset}' in {path.name} — add them "
                    f"(one row per native cell-type label)"))
        return out

    mapped = rows[rows["lineage"] != ""]
    excluded = rows[rows["lineage"] == ""]
    out.append(("registry:rows", "PASS",
                f"{len(rows)} rows (v{REGISTRY_VERSION}, fingerprint "
                f"{registry_fingerprint(path)}): {len(mapped)} mapped, "
                f"{len(excluded)} deliberately excluded"))

    dupes = rows["native_label"][rows["native_label"].duplicated()].tolist()
    out.append(("registry:unique_labels",
                "FAIL" if dupes else "PASS",
                f"duplicate native_label rows: {dupes}" if dupes
                else "no duplicate native_label rows"))

    bad_lineage = sorted(set(mapped["lineage"]) - set(schema.LINEAGES))
    out.append(("registry:lineage_vocabulary",
                "FAIL" if bad_lineage else "PASS",
                f"lineage value(s) outside {schema.LINEAGES}: {bad_lineage}"
                if bad_lineage else f"all lineage values in {schema.LINEAGES}"))

    # --- the grounding check: lineage must be a function of the CL term -----
    conflicts, unanchored = [], []
    for r in mapped.itertuples():
        if not r.cl_term_id:
            unanchored.append(r.native_label)
            continue
        anchor = schema.CL_LINEAGE_ANCHOR.get(r.cl_term_id)
        if anchor is None:
            unanchored.append(f"{r.native_label} ({r.cl_term_id})")
        elif anchor != r.lineage:
            conflicts.append(f"{r.native_label}: declared '{r.lineage}' but "
                             f"{r.cl_term_id} anchors to '{anchor}'")
    if conflicts:
        out.append(("registry:cl_grounding", "FAIL",
                    "declared lineage contradicts the CL anchor — "
                    + "; ".join(conflicts)))
    elif unanchored:
        out.append(("registry:cl_grounding", "WARN",
                    f"{len(unanchored)} mapped row(s) have no CL term in "
                    f"schema.CL_LINEAGE_ANCHOR, so their lineage is ungrounded: "
                    f"{unanchored[:6]}"))
    else:
        out.append(("registry:cl_grounding", "PASS",
                    f"all {len(mapped)} mapped rows carry a CL term whose anchor "
                    f"matches the declared lineage"))

    no_source = mapped.loc[mapped["source"] == "", "native_label"].tolist()
    out.append(("registry:citation",
                "WARN" if no_source else "PASS",
                f"{len(no_source)} mapped row(s) have no source citation: {no_source[:6]}"
                if no_source else f"all {len(mapped)} mapped rows carry a source citation"))

    unverified = mapped.loc[mapped["verified"].str.lower() != "yes", "native_label"].tolist()
    out.append(("registry:cl_verified",
                "WARN" if unverified else "PASS",
                f"{len(unverified)}/{len(mapped)} CL terms not yet confirmed against the "
                f"ontology (verified != yes). They are traceable but unproven — see "
                f"doc/CELLTYPE_MAPPING.md for the one-time confirmation procedure"
                if unverified else "all CL terms confirmed against the ontology"))

    no_reason = excluded.loc[excluded["notes"] == "", "native_label"].tolist()
    if len(excluded):
        out.append(("registry:exclusion_reasons",
                    "WARN" if no_reason else "PASS",
                    f"excluded type(s) with no recorded reason: {no_reason}" if no_reason
                    else f"all {len(excluded)} excluded types record why: "
                         f"{sorted(excluded['native_label'])}"))

    # --- coverage against what is actually in the raw data ------------------
    if native_labels is not None:
        present = set(map(str, native_labels))
        known = set(rows["native_label"])
        unregistered = sorted(present - known)
        stale = sorted(known - present)
        if unregistered:
            out.append(("registry:coverage", "FAIL",
                        f"{len(unregistered)} native type(s) in the raw data have NO "
                        f"registry row and would be dropped silently: {unregistered[:10]} "
                        f"— add them to {path.name} (with a lineage, or blank + a reason)"))
        else:
            out.append(("registry:coverage", "PASS",
                        f"all {len(present)} native types in the raw data have a registry row"))
        if stale:
            out.append(("registry:stale_rows", "WARN",
                        f"{len(stale)} registry row(s) name types not present in the raw "
                        f"data (harmless, but likely stale): {stale[:10]}"))

    return out
