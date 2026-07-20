# The Cell-Type → Lineage Mapping: De-hardcoding, Grounding, and Falsification

**Date:** 2026-07-20 · **Registry version:** 1.0.0 · **Marker vocabulary:** 1.0.0

This report documents a change to how the project decides that a native cell-type
label (`"CD8+ T cells"`, `"Tumor (Podo+)"`) belongs to a lineage
(`immune` / `tumour` / `stromal`), why the previous arrangement was unsafe, what
replaced it, and what the replacement found when run on UPMC and CRC.

The mapping is not a detail of ingest. **Three of the five headline spatial
features are defined in terms of it** — `kl_tumor`, `immune_tumor` and
`stroma_tumor` — at *both* taxonomy resolutions. Only `kl_mean` and
`self_enrich` are independent of it.

---

## 1. What was wrong

### 1.1 Five hardcoded copies, two incompatible spellings

| location | form |
|---|---|
| `schema.py` `UPMC_CANONICAL_CELLTYPES` | `str → "immune"/"tumour"/"stromal"` |
| `datasets/UPMC/adapter_config.py` `CELLTYPE_MAP` | same |
| `datasets/CRC/adapter_config.py` `CELLTYPE_MAP` | same |
| `src/validate_groups.py` `MY_GROUP` | `str → "IMMUNE"/"TUMOR"/"STROMA"` |
| `src/enrichment_features.py`, `src/marker_states.py` | `int CLUSTER_ID` lists |

Nothing asserted the copies agreed, and the validator's spelling
(`tumour`, lowercase) differed from the validator script's (`TUMOR`, uppercase),
so they could not even be diffed automatically.

### 1.2 The validator checked completeness, never correctness

The `celltype:mapping` check tested that the map was non-empty, that its values
were in the vocabulary, and that every native label was *a key in the dict*.
`{"tumor cells": "immune"}` passed all three. The validator never looked at
expression data when scoring the map.

### 1.3 The one expression-aware check was broken by string matching

`markers:lineage_validation` used exact equality (`"CD45" in marker_cols`). On CRC:

```
reported : markers:lineage_validation | WARN | 0/10 canonical lineage markers present: []
actual   : 'CD45 - hematopoietic cells:Cyc_4_ch_2'   ← in the panel
           'FOXP3 - regulatory T cells:Cyc_2_ch_3'   ← in the panel
```

For the only external cohort ingested, the sole safety net reported "cannot
verify" for a spurious reason and the map went in unchallenged. A WARN that fires
wrongly is worse than no check: it trains you to ignore it.

### 1.4 The lineage validator was orphaned and unrunnable

`spatial_positional_encoding/src/validate_groups.py` was imported by nothing,
never ran as part of ingest, carried a sixth private copy of the map, and pointed
at absolute paths in a *different project directory*
(`d:/Desktop/FYDP/FYDP 3/files/...`) — it could not execute in this repo at all.

Re-applying its own rule to its own saved output
(`outputs/spatial_maps/lineage_z.csv`) gave **1 FAIL / 16**: `APC` had CD45 below
the cross-cluster mean and PanCK above it.

### 1.5 Excluded types were a code comment

CRC dropped 17,831 cells (6.9%) across four types. The reason existed only as a
comment in `adapter_config.py` — it did not travel with the result.

---

## 2. What replaced it

```
celltype_registry.csv          ← the mapping, externalised. One cited row per
   (dataset, native_label)        (dataset, native_label). No dict anywhere.
        │
        ├─ registry.py            loader + structural validation
        │     lineage must agree with cl_term_id's anchor  → FAIL on conflict
        │     every raw label needs a row                  → FAIL on gap
        │     excluded rows must record a reason           → WARN if not
        │     evidence_override  = per-row, reasoned acceptance of a contradiction
        │
        ├─ schema.CL_LINEAGE_ANCHOR    Cell Ontology term → lineage
        ├─ markers.py                  normalisation + curated alias vocabulary
        └─ lineage_evidence.py         the falsifiability gate (runs at ingest)
                │
                ▼
   validator.py  →  validation_report.md + lineage_evidence.csv
                │
                ▼
   run_verify.py --perturb-map   →  does any CONCLUSION depend on a contested row?
```

### 2.1 Externalised and cited

`data_preprocessing/celltype_registry.csv`:

```
dataset,native_label,lineage,cl_term_id,cl_label,source,verified,notes,evidence_override
CRC,CD8+ T cells,immune,CL:0000625,"CD8-positive, alpha-beta T cell","Schurch et al. 2020, Cell 182(5):1341-1359",no,,
```

Adapters now read it — `CELLTYPE_MAP = registry.celltype_map(DATASET_NAME)`. The
registry reproduced both existing maps exactly (UPMC 16 mapped / CRC 25 mapped +
4 excluded), so the refactor is behaviour-preserving: both cohorts re-ingest to
identical cell counts.

### 2.2 Externally grounded

Lineage is a **function of the row's Cell Ontology term** via
`schema.CL_LINEAGE_ANCHOR`, not a free per-row choice. Two rows citing the same
CL term cannot be given different lineages; `registry:cl_grounding` FAILs if a
declared lineage disagrees with its term's anchor.

> **Stated limitation.** This is a curated lookup of the ~21 CL terms this
> project uses, **not a traversal of the Cell Ontology**. A real traversal would
> derive lineage from `is_a` ancestry and needs the OBO file. The table catches
> inconsistency and makes every row traceable to a resolvable term; it does not
> prove the term is the right one for that label. Every row therefore carries
> `verified=no` and the validator WARNs about it. See §6 for the one-time
> confirmation procedure.

### 2.3 Versioned

Each validation report carries `REGISTRY_VERSION` and a SHA-256 fingerprint of
the registry file, so a result is traceable to the exact mapping that produced it.

### 2.4 Falsifiable

`lineage_evidence.py` runs inside the validator on every ingest:

1. per native cluster, mean of each resolved lineage marker;
2. z-score each marker **across clusters**;
3. `lineage_score(cluster, L)` = **max** z over `L`'s resolved core markers;
4. **re-standardise each lineage score across clusters**;
5. `predicted` = argmax; `margin` = top − runner-up;
6. verdict vs the declared lineage: `AGREE` / `AMBIGUOUS` (margin < 0.25) /
   `CONTRADICTED`.

---

## 3. Two methodological corrections made during development

Both are recorded because the first version of the gate produced results that
looked like findings and were artifacts.

**(a) `mean` was the wrong aggregator.** Immune and stromal are heterogeneous
lineages — a T cell is CD45+/CD20−/CD68−, an endothelial cell is
CD31+/aSMA−/CollagenIV−. Averaging a lineage's subset anchors penalises every
subset for not being the others. Changed to **max**: "bright for at least one
anchor of this lineage". This alone moved CRC from 10 to 17 `AGREE`, correctly
resolving granulocytes, all macrophage subsets, and DCs.

**(b) Scores from panels of different sizes are not comparable.** A score built
from 1 marker and one built from 5 have different spreads, so the argmax was
partly a panel-size artifact. Fixed by re-standardising each lineage score across
clusters before comparing (step 4).

**A panel-specificity bug found this way:** MUC1 was in the tumour core panel and
is a genuine *plasma cell* marker — it made CRC's plasma-cell cluster read as
tumour. Moved to `supporting`.

**Guarding against tuning-to-agree.** After fixing the aggregator, four CRC
contradictions remained, two of which were *a priori gaps in my panel*: the
immune core covered T/B/myeloid/DC anchors but had **no NK anchor and no plasma
anchor**, while both cohorts carry CD56 and CD38. Adding them was decided from
which leukocyte subsets the label vocabulary contains — not from which fix made a
disagreement disappear. The test that this was not tuning: **NK cells flipped to
AGREE, plasma cells did not.** Adding CD38 raised the plasma-cell immune score
from −0.90 to +1.17 and the verdict still stood, because its tumour score is
+3.01. The gate remained able to disagree with me.

---

## 4. Results

### 4.1 Marker resolution — the broken check, fixed

| check | CRC before | CRC after | UPMC after |
|---|---|---|---|
| `markers:lineage_validation` | **0/10** | **27/27** | 21/27 |
| `markers:node_features` | **0/8** | **8/8** | 8/8 |
| `markers:recommended_pair` | **none** | FoxP3 + CD45RO | FoxP3 + CD45RO |

UPMC's 21/27 is correct, not a failure: its 39-plex panel genuinely lacks CDX2,
MUC1, EGFR, p53, CD44 and CD138. The check now names what is missing instead of
reporting a bare count.

Collision guard verified: `CD4` resolves to CD4 (not CD45), `CD3` to CD3 (not CD31).

### 4.2 The falsifiability gate

**UPMC — 3 of 16 clusters contradicted, 13.5% of cells (materially blocking):**

| cluster | declared | evidence | margin | cells |
|---|---|---|---|---|
| Tumor (Podo+) | tumour | stromal | +0.89 | 18,945 |
| APC | immune | tumour | +1.20 | 11,419 |
| Naive immune cell | immune | stromal | +1.51 | 10,059 |

These are **the same three clusters** flagged by hand from the old z-matrix
before the gate existed (APC's CD45<PanCK; Naive immune's CD31 z=+1.90;
Tumor (Podo+)'s Podoplanin=lymphatic). Independent corroboration by a different
scoring rule.

**CRC — 4 of 25 clusters contradicted, 3.9% of cells (below the 5% threshold):**

| cluster | declared | evidence | margin | cells |
|---|---|---|---|---|
| plasma cells | immune | tumour | +1.84 | 8,510 |
| CD11b+ monocytes | immune | tumour | +0.41 | 815 |
| CD4+ T cells GATA3+ | immune | tumour | +0.53 | 67 |
| CD163+ macrophages | immune | stromal | +0.59 | 38 |

The gate **correctly blocked** the UPMC ingest. Each contradiction was then
accepted with a written per-row justification in `evidence_override`, which
downgrades it from blocking to reported — it is still counted, still listed,
still shown in the per-cluster table, and still carried into the perturbation
analysis. `--force` was deliberately not used: it is global and silent, and using
it once teaches you to use it always.

### 4.3 Sensitivity — does any conclusion depend on a contested row?

`run_verify.py --dataset <NAME> --taxonomy lineage --perturb-map` re-runs the
full matrix with each contested cluster **dropped**, then **flipped** to the
lineage the evidence predicted, then all flipped at once.

| cohort | scenarios | verdict changes |
|---|---|---|
| UPMC | 8 | **0** |
| CRC | 9 | **0** |
| CRC_doublets | 11 | **0** |

**28 scenarios across three cohorts, zero verdict changes.** Every feature stays
`STRONG (real + adds + stable)` under every perturbation. No conclusion in either
matrix rests on a contested cell-type assignment.

That is the claim the mapping work exists to support, and it is stronger than
"the mapping is right": it does not require the mapping to be right. What the
perturbations *do* move is effect sizes — see §4.4 and §4.5.

### 4.4 The doublet-exclusion bias — quantified

CRC drops two ambiguous doublet clusters (3,950 cells): `"tumor cells / immune
cells"` and `"immune cells / vasculature"`. These sit at the immune–tumour
interface, which is exactly what `immune_tumor` measures — so the exclusion is
**not neutral with respect to the quantity being measured.**

A variant cohort `CRC_doublets` (same raw file, doublets mapped instead of
dropped: 244,504 cells vs 240,554 = +3,950 exactly) makes the comparison direct:

| feature | CRC (doublets excluded) | CRC_doublets (included) | \|ratio\| |
|---|---|---|---|
| kl_mean | +118.30 | +137.31 | 0.86 |
| kl_tumor | +133.28 | +160.80 | 0.83 |
| self_enrich | +10.65 | +11.94 | 0.89 |
| **immune_tumor** | **−11.20** | **−8.38** | **1.34** |
| stroma_tumor | −5.49 | −6.80 | 0.81 |

**Verdicts are identical — all five STRONG in both.** But the effect sizes are
not, and the direction is the predicted one: excluding the interface doublets
**inflates \|immune_tumor\| by 34%**, while *deflating* all four other features
(0.81–0.89×). `immune_tumor` is the only feature that moves the wrong way, and it
is the one that measures the compartment the dropped cells occupy.

**Honest reading:** the conclusion "immune–tumour contact is spatially structured
beyond chance" is robust to the exclusion. The *magnitude* of that structure is
not — it is overstated by roughly a third in the production CRC numbers. Any
claim about effect size, or any cross-cohort comparison of `immune_tumor`
magnitudes, must state which convention was used.

### 4.5 UPMC sensitivity

Eight scenarios (baseline, 3 drops, 3 flips, evidence-all). **No verdict changed:**
all five features remain `STRONG (real + adds + stable)` throughout.

`null_z_median`:

| scenario | kl_mean | kl_tumor | self_enrich | immune_tumor | stroma_tumor |
|---|---|---|---|---|---|
| baseline | +1484.2 | +2624.0 | +34.4 | −58.6 | −58.1 |
| drop:APC | +1683.3 | +2494.3 | +36.6 | −62.9 | −63.0 |
| drop:Naive immune cell | +1683.8 | +2540.2 | +37.1 | −60.7 | −62.9 |
| drop:Tumor (Podo+) | +1979.4 | +3014.8 | +38.4 | −61.5 | −59.3 |
| flip:APC→tumour | +1606.9 | +3018.9 | +34.6 | −61.9 | −62.8 |
| flip:Naive immune cell→stromal | +2062.6 | +2676.4 | +42.3 | −54.9 | −72.2 |
| flip:Tumor (Podo+)→stromal | +1509.2 | +2000.7 | +43.5 | −60.9 | **−39.4** |
| evidence-all | +1828.4 | +2490.4 | +46.0 | −53.8 | −51.5 |

Two observations worth stating rather than burying:

- **`stroma_tumor` moves most**, −58.1 → **−39.4** (a 32% reduction) when
  `Tumor (Podo+)` is reassigned to stromal. That is expected — moving a
  19k-cell cluster from tumour to stromal directly changes what a
  stromal↔tumour contact statistic is computed over. The verdict holds
  (\|z\|≫2), but the magnitude is conditional on that one contested row.
- **The closest call in either cohort** is UPMC's `self_enrich`
  `spatial_specific`: **0.551** at baseline against a **≥0.5** threshold, falling
  to **0.523** under `flip:APC→tumour`. It never crosses, but the margin is
  0.023 — this is the single result most exposed to the mapping, and it should
  be reported with that caveat rather than as a flat pass.

Stability is essentially unaffected everywhere (r = 0.971–0.991 across all
scenarios and features).

---

## 5. Does this generalize to a new cohort?

| layer | generalizes automatically? | if not, what is manual |
|---|---|---|
| Marker decoration/case (`CD45 - x:Cyc_4_ch_2`, `FOXP3`) | **yes** — mechanical | — |
| Marker aliases (`Cytokeratin`≡`PanCK`) | **no** | add to `markers.SYNONYMS`, bump version |
| Registry structure, CL grounding, citation, coverage | **yes** | — |
| A native label → lineage decision | **no, and never can be** | one registry row, cited |
| CL term for a novel cell type | **no** | add to `schema.CL_LINEAGE_ANCHOR` |
| Evidence gate scoring | **yes**, if each lineage has ≥1 resolvable core marker | else it **refuses to run** and names the missing lineage rather than reporting a confident wrong answer |
| Lineage core panels for an exotic panel | **no** | add to `schema.LINEAGE_MARKER_PANELS` |
| Perturbation scenarios | **yes** — derived from the registry's own contested rows | — |

The irreducible manual work for a new cohort is: **one registry row per native
label** (lineage + CL term + citation), plus any novel marker aliases. Everything
else is enforced or derived. Critically, an *incomplete* registry now **FAILs**
(`registry:coverage`) instead of silently shrinking the cohort.

---

## 6. Outstanding: confirming the CL terms

All 45 mapped rows carry `verified=no`. The terms were assigned best-effort and
are traceable but unproven. One-time procedure per distinct term (~21 terms, not
45 rows):

1. Open `https://www.ebi.ac.uk/ols4/ontologies/cl/classes?obo_id=CL:0000625`.
2. Confirm the label matches `cl_label` and that the native label is a reasonable
   instance of it.
3. Confirm the term's ancestry is consistent with `schema.CL_LINEAGE_ANCHOR`.
4. Set `verified=yes`; bump `REGISTRY_VERSION`.

Terms most worth checking first, because they are the least certain:

- `CL:0000125` (glial cell) for CRC `nerves` — enteric glia are neural-crest
  derived and **not** mesenchymal stroma. Assigned stromal only as the closest of
  three available lineages. The 3-lineage vocabulary has no correct home for it.
- `CL:0000145` (professional antigen presenting cell) for UPMC `APC` — a
  functional grouping, not a marker-derived identity, and the cluster the
  evidence contradicts most confidently.
- `CL:0000738` (leukocyte) for UPMC `Naive immune cell` / CRC `immune cells` —
  generic terms standing in for labels that name no subset.

---

## 7. Files

| file | role |
|---|---|
| `data_preprocessing/celltype_registry.csv` | **the mapping** — cited, versioned, CL-grounded |
| `data_preprocessing/registry.py` | loader + structural/provenance validation |
| `data_preprocessing/markers.py` | marker normalisation + alias vocabulary |
| `data_preprocessing/lineage_evidence.py` | the falsifiability gate |
| `data_preprocessing/schema.py` | `CL_LINEAGE_ANCHOR`, `LINEAGE_MARKER_PANELS` |
| `data_preprocessing/datasets/*/lineage_evidence.csv` | per-cluster evidence (machine-readable) |
| `data_preprocessing/datasets/CRC_doublets/` | exclusion-bias variant cohort |
| `spatial_positional_encoding/run_verify.py` | `--perturb-map` sensitivity analysis |
| ~~`spatial_positional_encoding/src/validate_groups.py`~~ | **removed** — superseded by `lineage_evidence.py` |

### Commands

```bash
cd data_preprocessing
python run_ingest.py --dataset UPMC            # gate runs inside validation
python run_ingest.py --dataset CRC
python run_ingest.py --dataset CRC_doublets    # exclusion-bias variant

cd ../spatial_positional_encoding
python run_verify.py --dataset UPMC --taxonomy lineage --perturb-map
python run_verify.py --dataset CRC  --taxonomy lineage --perturb-map
```
