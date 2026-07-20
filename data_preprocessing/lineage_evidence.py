"""
lineage_evidence.py — Confront the cell-type -> lineage registry with the marker data.

WHY THIS EXISTS
---------------
Before this, the only expression-aware check on the lineage map lived in
spatial_positional_encoding/src/validate_groups.py, which was:
  - imported by nothing (never ran as part of ingest),
  - hardcoded to absolute paths in a DIFFERENT project directory, so it could
    not execute in this repo at all,
  - carrying a SIXTH private copy of the map, in a different spelling,
  - and asking only "is the lineage's marker above the cross-type mean?" — a
    one-sided threshold that cannot say which lineage the evidence actually
    favours.

This module replaces it and runs inside validator.py, so no cohort can be
ingested without its registry being confronted with its own expression data.

THE METHOD
----------
1. Per native cluster, take the MEAN of each resolved lineage marker.
2. z-score each marker ACROSS clusters — "is this cluster high for this marker
   relative to the other clusters in this cohort?". This is what makes the score
   comparable between cohorts with different dynamic ranges and normalisations.
3. lineage_score(cluster, L) = MAX z over L's resolved CORE markers
   (schema.LINEAGE_MARKER_PANELS). Supporting markers are reported, never scored.

   MAX, not mean, because immune and stromal are heterogeneous lineages: a T cell
   is CD45+/CD20-/CD68-, an endothelial cell is CD31+/aSMA-/CollagenIV-. Averaging
   across a lineage's subset anchors penalises every subset for not being the
   others, which is not what "does this cluster belong to lineage L" asks.

4. RE-STANDARDISE each lineage score across clusters (z again), so the three are
   commensurable before they are compared. Without this the argmax is unsound:
   a score built from 1 marker and a score built from 5 have different spreads,
   and whichever lineage happens to own more markers wins by construction. This
   step is what makes step 5 a fair comparison rather than a panel-size artifact.

5. predicted = argmax over the three lineages; margin = top - runner-up.
6. Verdict vs the registry's declared lineage:
       AGREE        predicted == declared
       AMBIGUOUS    margin < MARGIN_THRESHOLD — evidence does not separate the
                    lineages for this cluster, whichever way it points
       CONTRADICTED predicted != declared AND margin >= MARGIN_THRESHOLD

argmax, not a threshold: "CD45 is above average" cannot distinguish an immune
cluster from a cluster that is simply bright everywhere. "CD45 evidence beats
PanCK and aSMA evidence, by this margin" can.

DOES THIS GENERALISE?
---------------------
It needs at least one resolvable CORE marker for EACH of the three lineages. If
any lineage has none, an argmax over the remaining lineages would be rigged
against the unscoreable one, so the check REFUSES to run and says which lineage
it lacked — rather than reporting a confident wrong answer. Panels are canonical
names resolved through markers.py, so a cohort naming markers differently is
handled automatically; a cohort with a genuinely different panel needs its
markers added to schema.LINEAGE_MARKER_PANELS by hand.

WHAT A CONTRADICTION MEANS
--------------------------
Not automatically that the registry is wrong. Marker means over a mixed cluster,
segmentation spillover between touching cells, and panels missing a lineage's
best marker all produce genuine disagreement. It means the row is not supported
by this cohort's own data and must be either justified in the registry `notes`
or changed — and it should be carried into the perturbation analysis
(run_verify.py --perturb-map) to show whether the final result depends on it.
"""
from __future__ import annotations

from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

import markers as marker_utils
import schema

# Below this separation between the best and second-best lineage score, the
# evidence is treated as not deciding. 0.25 SD of the cross-cluster marker
# distribution: small enough that clear cases still resolve, large enough that
# a coin-flip is never reported as a contradiction.
MARGIN_THRESHOLD = 0.25

# Contradictions covering at least this fraction of mapped cells FAIL the ingest;
# smaller ones WARN. See summarise() for why the gate is graded rather than binary.
MATERIAL_CONTRADICTION_FRAC = 0.05

# Cap on cells used for the per-cluster means. The means are stable well below
# this; the cap keeps validation fast on million-cell cohorts. Fixed seed.
MAX_CELLS = 300_000
SEED = 1029


def _scoreable_panels(marker_cols: List[str]) -> Tuple[Dict[str, Dict[str, str]], Dict[str, List[str]]]:
    """Resolve each lineage's core panel against this cohort's actual columns."""
    resolved, missing = {}, {}
    for lineage, panel in schema.LINEAGE_MARKER_PANELS.items():
        found, miss = marker_utils.resolve_panel(panel["core"], marker_cols)
        resolved[lineage] = found
        missing[lineage] = miss
    return resolved, missing


def evaluate(
    locations: pd.DataFrame,
    expression: pd.DataFrame,
    marker_cols: List[str],
    declared_map: Dict[str, str],
    max_cells: int = MAX_CELLS,
) -> Tuple[Optional[pd.DataFrame], str]:
    """Score every native cluster against the lineage marker panels.

    Both frames must already carry the canonical `acquisition_id` / `cell_id`
    columns; `locations` must carry `cluster_label`.

    Returns (table, note). `table` is None when the check cannot run, and `note`
    always explains what happened.
    """
    resolved, missing = _scoreable_panels(marker_cols)
    unscoreable = [L for L, found in resolved.items() if not found]
    if unscoreable:
        return None, (
            f"cannot run: no core marker resolved for lineage(s) {unscoreable}. "
            f"An argmax over the remaining lineages would be rigged against them. "
            f"Core panels sought: "
            + "; ".join(f"{L}={schema.LINEAGE_MARKER_PANELS[L]['core']}" for L in unscoreable)
            + f". Add this cohort's equivalent marker to schema.LINEAGE_MARKER_PANELS "
              f"or its alias to markers.SYNONYMS."
        )

    keys = [schema.ACQ_COL, schema.CELL_COL]
    loc = locations[keys + [schema.CLUSTER_LABEL_COL]].copy()
    if len(loc) > max_cells:
        loc = loc.sample(n=max_cells, random_state=SEED)

    used_cols = sorted({col for found in resolved.values() for col in found.values()})
    expr = expression[keys + used_cols]
    df = loc.merge(expr, on=keys, how="inner")
    if df.empty:
        return None, "cannot run: no cells survived the locations/expression join"

    df[schema.CLUSTER_LABEL_COL] = df[schema.CLUSTER_LABEL_COL].astype(str)
    grouped = df.groupby(schema.CLUSTER_LABEL_COL)
    means = grouped[used_cols].mean()
    n_cells = grouped.size()

    if len(means) < 3:
        return None, (f"cannot run: only {len(means)} native cluster(s) present; "
                      f"z-scoring markers across clusters needs at least 3")

    # z-score each marker ACROSS clusters
    sd = means.std(ddof=0).replace(0, np.nan)
    z = (means - means.mean()) / sd
    z = z.fillna(0.0)

    # Step 3 — raw lineage score = MAX z over the lineage's resolved core markers
    # ("bright for at least one anchor of this lineage").
    raw = pd.DataFrame(
        {L: z[list(found.values())].max(axis=1) for L, found in resolved.items()},
        index=means.index,
    )
    # Step 4 — re-standardise each lineage score across clusters so a 1-marker
    # score and a 5-marker score are on the same scale before they are compared.
    raw_sd = raw.std(ddof=0).replace(0, np.nan)
    lin_z = ((raw - raw.mean()) / raw_sd).fillna(0.0)

    rows = []
    for cluster in means.index:
        scores = {L: float(lin_z.loc[cluster, L]) for L in resolved}
        order = sorted(scores, key=scores.get, reverse=True)
        predicted, runner_up = order[0], order[1]
        margin = scores[predicted] - scores[runner_up]
        declared = declared_map.get(cluster, "")

        if not declared:
            verdict = "EXCLUDED"          # registry deliberately drops this type
        elif margin < MARGIN_THRESHOLD:
            verdict = "AMBIGUOUS"
        elif predicted == declared:
            verdict = "AGREE"
        else:
            verdict = "CONTRADICTED"

        rows.append({
            "native_label": cluster,
            "n_cells": int(n_cells[cluster]),
            "declared": declared or "(excluded)",
            "predicted": predicted,
            "margin": round(margin, 3),
            "z_immune": round(scores["immune"], 3),
            "z_tumour": round(scores["tumour"], 3),
            "z_stromal": round(scores["stromal"], 3),
            "verdict": verdict,
        })

    table = pd.DataFrame(rows).sort_values(
        ["verdict", "n_cells"], ascending=[True, False]).reset_index(drop=True)

    used = {L: sorted(found) for L, found in resolved.items()}
    note = (f"scored {len(df):,} cells over {len(means)} native clusters; "
            f"core markers used: {used}")
    if any(missing.values()):
        note += f"; core markers not found in this panel (not scored): " \
                f"{ {L: m for L, m in missing.items() if m} }"
    return table, note


def summarise(table: pd.DataFrame,
              overrides: Optional[Dict[str, str]] = None) -> Tuple[str, str]:
    """Collapse the per-cluster table into a (status, detail) validator check.

    Materiality: a contradiction is graded by how much of the cohort it affects,
    not merely by existing. Blocking an ingest on a 38-cell cluster would be
    disproportionate and would train users to reach for --force, which is how a
    gate stops being a gate. Contradictions are ALWAYS listed in full whatever
    the status; only whether they block differs.
    """
    mapped = table[table["verdict"] != "EXCLUDED"]
    bad = mapped[mapped["verdict"] == "CONTRADICTED"]
    amb = mapped[mapped["verdict"] == "AMBIGUOUS"]
    ok = mapped[mapped["verdict"] == "AGREE"]

    if len(bad):
        overrides = overrides or {}
        n_mapped = max(int(mapped["n_cells"].sum()), 1)
        accepted = bad[bad["native_label"].isin(overrides)]
        open_ = bad[~bad["native_label"].isin(overrides)]
        # Only UNJUSTIFIED contradictions can block: an accepted one has been
        # seen and reasoned about in the registry, which is what the gate wants.
        frac_open = int(open_["n_cells"].sum()) / n_mapped
        frac_all = int(bad["n_cells"].sum()) / n_mapped

        def _list(sub):
            return "; ".join(f"{r.native_label} (declared {r.declared}, evidence "
                             f"{r.predicted}, margin {r.margin:+.2f}, {r.n_cells:,} cells)"
                             for r in sub.itertuples())

        status = "FAIL" if frac_open >= MATERIAL_CONTRADICTION_FRAC else "WARN"
        parts = [f"marker evidence CONTRADICTS the registry for {len(bad)}/{len(mapped)} "
                 f"mapped clusters covering {frac_all:.1%} of mapped cells"]
        if len(open_):
            parts.append(f"UNJUSTIFIED ({frac_open:.1%} of cells, "
                         f"{'at or above' if status == 'FAIL' else 'below'} the "
                         f"{MATERIAL_CONTRADICTION_FRAC:.0%} materiality threshold): {_list(open_)}")
        if len(accepted):
            parts.append(f"ACCEPTED with a recorded justification in the registry "
                         f"`evidence_override` ({int(accepted['n_cells'].sum()) / n_mapped:.1%} "
                         f"of cells): {_list(accepted)}")
        parts.append("every contradiction, accepted or not, is carried into "
                     "run_verify.py --perturb-map — acceptance is a reason, not a result")
        detail = ". ".join(parts)
    elif len(amb):
        detail = (f"{len(ok)}/{len(mapped)} clusters confirmed by marker evidence; "
                  f"{len(amb)} not separable (margin < {MARGIN_THRESHOLD}): "
                  + ", ".join(f"{r.native_label} (declared {r.declared}, evidence "
                              f"{r.predicted}, margin {r.margin:+.2f})"
                              for r in amb.itertuples()))
        status = "WARN"
    else:
        detail = f"all {len(mapped)} mapped clusters confirmed by marker evidence"
        status = "PASS"
    return status, detail
