# Results

*Spatial features from in-situ proteomics imaging: what was measured, what held up,
and what did not.*

This report contains **results only**. For how the pipeline works — the code path,
the feature definitions, the test definitions — see
[PIPELINE_EXPLAINED.md](PIPELINE_EXPLAINED.md). For how to re-derive every number
here, see §9.

---

## 1. Provenance

Every number in this report comes from **one clean end-to-end run** on 2026-07-20:
all three cohorts re-ingested from raw, then every matrix recomputed by the current
code. Nothing is carried over from an earlier run.

| | |
|---|---|
| Cell-type registry | `celltype_registry.csv` v1.0.0 |
| Marker vocabulary | `markers.py` v1.0.0 |
| Seed | 1029 everywhere (ingest sampling, permutation null, split-half, RSF) |
| Permutations | 20 per sample (5 per perturbation scenario) |
| Determinism | verified — re-running any matrix with the same seed reproduces it byte-for-byte |

> **These numbers supersede the previously committed CSVs.** One check
> (§7.2) disagrees materially with the version committed before this run, and the
> provenance of that older file could not be established. That is the reason the
> whole pipeline was re-run from raw rather than partially refreshed.

---

## 2. Cohorts

| | UPMC | CRC | CRC_doublets |
|---|---|---|---|
| role | primary | **external validation** | measurement control, *not* a third cohort |
| tissue | head & neck (HNSCC) | colorectal, advanced stage | — same raw file as CRC — |
| source | this project's baseline cohort | Schürch et al. 2020 | |
| samples | 308 | 140 | 140 |
| cells | 2,061,102 | 240,554 | 244,504 |
| patients | 81 | 35 | 35 |
| native cell types | 16 | 29 (25 mapped) | 29 (27 mapped) |
| markers | 39 | 56 | 56 |
| survival | **yes** — 103 events, 33.4% | **none** | none |

Three further folders exist under `data_preprocessing/datasets/` — `keren_tnbc`,
`hubmap_intestine_codex`, `Ferguson` — with adapter configs but **no `processed/`
directory**. They were never successfully ingested (`Ferguson`'s adapter is still
an unedited `_template` copy) and contribute nothing to any number in this report.
`python run_pipeline.py --list` is the authoritative inventory.

CRC_doublets is CRC with its two ambiguous doublet clusters *mapped* instead of
dropped. It exists to test whether excluding interface cells manufactures the
result (§6.2). It is a variant of the same tissue and must never be reported as
independent evidence.

**Cohort labelling note.** UPMC is the head & neck cohort. Schürch et al. 2020 is
the *source of the CRC data* and, separately, the *literature source for the node
marker shortlist*. Three legacy documents describe UPMC itself as a "colorectal
CODEX atlas (Schürch 2020)"; that is a documentation error, corrected here and
flagged in §8.

---

## 3. What is being tested

Two feature blocks, both computed per sample, both tested **without any outcome
data**:

| block | what it is | features |
|---|---|---|
| **Enrichment** | abundance-corrected spatial organisation, read off a Delaunay neighbour graph. 0 = a random arrangement. | `kl_mean`, `kl_tumor`, `self_enrich`, `immune_tumor`, `stroma_tumor` |
| **Node markers** | the mean of one protein marker over *only* the cells of the type it is biologically read in | `tumor_ki67`, `tumor_mac_pdl1`, `cd8_granzymeb`, `cd4_pd1`, `cd4_icos`, `cd4_foxp3`, `tcell_cd45ro`, `apc_mac_hladr` |

Both are tested against the same baseline — **composition**, the fraction of each
cell type in the sample. That is the boring explanation ("you only need to know
what cells are there, not where they are or what they're expressing"), and a
feature only counts if it beats it.

Three label-free tests, pass bars in brackets: **signal** vs a within-sample
permutation null [`|z| ≥ 2`], **specificity** beyond composition [`≥ 0.5`], and
**stability** on independent half-samples [`≥ 0.5`]. All three must pass for
`STRONG`.

---

## 4. Result 1 — the enrichment block

`verify_<taxonomy>.csv`. Verdict is **STRONG on all 5 features, in all 4
cohort × taxonomy combinations.**

**Lineage taxonomy** (3 categories — portable, directly comparable across cohorts):

| feature | UPMC z | CRC z | UPMC spec / stab | CRC spec / stab |
|---|---|---|---|---|
| `kl_mean` | +1397.5 | +81.3 | 0.77 / 0.98 | 0.94 / 0.96 |
| `kl_tumor` | +2141.1 | +92.7 | 0.93 / 0.99 | 0.99 / 0.96 |
| `self_enrich` | +31.4 | +7.6 | 0.55 / 0.98 | 1.00 / 0.62 |
| `immune_tumor` | −51.6 | −9.4 | 1.00 / 0.99 | 0.87 / 0.88 |
| `stroma_tumor` | −52.2 | −4.6 | 0.83 / 0.97 | 0.85 / 0.81 |

**Native taxonomy** (16 types UPMC / 25 CRC — finer, per-dataset):

| feature | UPMC z | CRC z | UPMC spec / stab | CRC spec / stab |
|---|---|---|---|---|
| `kl_mean` | +100.1 | +6.4 | 0.75 / 0.93 | 1.00 / 0.64 |
| `kl_tumor` | +1250.5 | +57.4 | 0.58 / 0.99 | 1.00 / 0.89 |
| `self_enrich` | +16.6 | +6.5 | 0.67 / 0.79 | 1.00 / 0.65 |
| `immune_tumor` | −51.6 | −9.4 | 0.69 / 0.99 | 1.00 / 0.88 |
| `stroma_tumor` | −52.2 | −4.6 | 0.60 / 0.97 | 1.00 / 0.81 |

Reading these:

- **Signs are biologically coherent and consistent across cohorts.** `kl_mean` and
  `kl_tumor` positive (neighbourhoods differ from the bulk mix; tumour sits in its
  own niche); `immune_tumor` and `stroma_tumor` negative on both cohorts (immune
  and stromal cells contact tumour *less* than chance — exclusion, not
  infiltration).
- **`self_enrich` roughly quadruples at native resolution** (UPMC 0.94 → 2.21,
  CRC 0.78 → 2.40). Finer labels reveal finer self-segregation, exactly as they
  should. This is a sanity check the pipeline passes.
- **Magnitudes differ by two orders of magnitude between cohorts** (UPMC z in the
  thousands, CRC in the tens) and this is expected, not alarming: UPMC has ~15× the
  cells per sample, and z grows with sample size. Both are decisively past the bar.
- **CRC's `stroma_tumor` (z = −4.6, stability 0.81) is the weakest cell in the
  table** and is the one to watch if signal ever becomes marginal.

---

## 5. Result 2 — the node marker block

`verify_nodes_<taxonomy>.csv`. Support is essentially complete (UPMC 308/308 samples
for every feature; CRC ≥ 139/140).

`cond_z` here answers a different question from the enrichment block's `null_z`:
permuting cell-type labels makes the conditioning set a size-matched random draw,
so the null is roughly the bulk mean of that marker. **cond_z therefore asks
whether the marker is actually enriched in the cell type the feature is named for
— the premise celltype-conditioning rests on, and one that had never been tested.**
It is not a spatial claim; node features use no graph.

| feature | marker | UPMC cond_z | CRC cond_z | UPMC spec / stab | CRC spec / stab |
|---|---|---|---|---|---|
| `tumor_ki67` | Ki67 | **+20.7** | **+5.7** | 0.89 / 1.00 | 0.92 / 0.97 |
| `tumor_mac_pdl1` | PDL1 | **−7.6** ⚠ | +0.4 ⚠ | 0.99 / 1.00 | 1.00 / 0.98 |
| `cd8_granzymeb` | GranzymeB | **+23.3** | +0.3 ⚠ | 0.99 / 0.98 | 0.98 / 0.86 |
| `cd4_pd1` | PD1 | **+9.0** | **+5.4** | 0.96 / 0.97 | 1.00 / 0.96 |
| `cd4_icos` | ICOS | **+25.7** | **+7.3** | 1.00 / 1.00 | 1.00 / 0.97 |
| `cd4_foxp3` | FoxP3 | **+4.3** | **+5.9** | 0.99 / 0.99 | 0.89 / 0.89 |
| `tcell_cd45ro` | CD45RO | **+48.4** | **+14.0** | 0.97 / 0.99 | 0.99 / 0.96 |
| `apc_mac_hladr` | HLA-DR | **+40.5** | **+3.4** | 0.96 / 0.99 | 1.00 / 0.94 |

**Verdicts: UPMC 7/8 STRONG + 1 CONTRADICTED; CRC 6/8 STRONG + 2 weak.**

### 5.1 What holds

- **Composition-specificity is 0.89–1.00 on every feature, on both cohorts.** The
  cell-type mix explains almost none of the marker state. Sensible in hindsight —
  *how much* of a cell type is present and *what those cells express* are close to
  orthogonal axes — and it is exactly the claim the block needs.
- **Stability is 0.86–1.00 throughout.** The conditioned means are reproducible on
  independent halves of the same tissue.
- **`cd4_foxp3` and `tcell_cd45ro` — the two markers behind the project's best
  legacy survival configuration — are STRONG on both cohorts.** CRC has no survival
  and cannot corroborate a C-index, but it independently corroborates that these
  two features measure a real, celltype-enriched, reproducible, non-abundance
  quantity. That is the part of the claim external validation *can* reach without
  outcome data, and it now has it.

### 5.2 Two features fail, and the failures are informative

- **`tumor_mac_pdl1` fails on both cohorts.** On UPMC `cond_z = −7.6`: PDL1 is
  *depleted* in tumour+macrophage cells relative to a size-matched random draw —
  the opposite of the feature's premise. On CRC it is flat (+0.4). Its
  composition-specificity (0.99/1.00) and stability (1.00/0.98) are near-perfect,
  so it is a precise, reproducible measurement of *something* — just not of what
  its name asserts. Most likely "tumour + macrophage" is simply too wide a set for
  a broadly-expressed marker to be enriched against everything else.
  **Recommendation: do not use without re-deriving the conditioning.**
- **`cd8_granzymeb` is cohort-specific.** Strongly enriched on UPMC (+23.3), absent
  on CRC (+0.3, only 35% of samples clearing |z| > 2). A feature that behaves this
  differently across two cohorts is not portable, whatever it does on either alone.

Neither failure is visible from a survival result. Neither would have been found
without this matrix.

### 5.3 A caveat on cross-cohort comparison

`real_mean` is **not** comparable between cohorts: UPMC arrives already
arcsinh-normalised from source, CRC is arcsinh-normalised at ingest with
cofactor 5. The three test metrics are the comparable quantities; the raw means
are not.

---

## 6. Result 3 — robustness of the enrichment result

### 6.1 The cell-type → lineage mapping does not carry the result

Three of the five enrichment features are *defined* in terms of the lineage mapping,
so `--perturb-map` re-runs the whole matrix with every contested assignment dropped
and then flipped to whatever the marker evidence predicted.

Contested rows (registry `evidence_override` — where marker evidence disagreed with
the declared lineage and a human accepted it anyway):

| cohort | contested clusters |
|---|---|
| UPMC | `APC`, `Naive immune cell`, `Tumor (Podo+)` |
| CRC | `plasma cells`, `CD11b+ monocytes`, `CD4+ T cells GATA3+`, `CD163+ macrophages` |

**Result: no verdict changed in any scenario, on either cohort** — including
`evidence-all`, which flips every contested cluster simultaneously. Magnitudes do
move (flipping UPMC's `Tumor (Podo+)` to stromal lifts `kl_mean`'s
composition-specificity from 0.77 to 0.92; CRC's `flip:plasma cells→tumour` pushes
`stroma_tumor` from z −5.5 to −8.4), but every conclusion survives.

### 6.2 Excluding the doublets does not manufacture the result

Production CRC drops 4 native types (6.9% of cells), two of them doublets sitting
at the immune–tumour interface — precisely where three of the five features
measure. `CRC_doublets` maps them instead.

| feature | CRC z | CRC_doublets z |
|---|---|---|
| `kl_mean` | +81.3 | +92.4 |
| `kl_tumor` | +92.7 | +101.4 |
| `self_enrich` | +7.6 | +8.5 |
| `immune_tumor` | −9.4 | −8.6 |
| `stroma_tumor` | −4.6 | −4.4 |

All verdicts identical, node-feature verdicts identical too. The exclusion is not
load-bearing.

---

## 7. Result 4 — survival (UPMC only)

Survival is **not** how the feature blocks are validated — §4 and §5 are, and they
need no outcome data. This section answers the separate, narrower question: *do the
validated features also predict outcome on the one cohort that has any?*

Protocol: RandomSurvivalForest (100 trees, `random_state=1029`), GroupKFold by
`patient_id`, 10 folds, 308 samples / 103 events / 81 patients. All feature sets
are scored **on the same folds**, so deltas are paired.

### 7.1 The result depends on how strong the baseline is

| feature set | C (lineage baseline, 3 cats) | Δ vs base | C (native baseline, 16 cats) | Δ vs base |
|---|---|---|---|---|
| Composition baseline | 0.605 | — | **0.635** | — |
| Baseline + Enrichment | 0.626 | +0.021 | **0.652** | +0.017 |
| Baseline + Nodes | **0.655** | +0.051 | 0.634 | **−0.001** |
| Baseline + Enrichment + Nodes | 0.622 | +0.017 | **0.659** | +0.023 |
| Baseline + Noise (control) | 0.592 | −0.013 | 0.629 | −0.006 |
| Enrichment alone | 0.586 | −0.019 | 0.574 | −0.061 |
| Nodes alone | 0.603 | −0.002 | 0.603 | −0.032 |

**The single most important line in this table is `Baseline + Nodes`.** Against the
weak 3-category baseline it looks like the best block in the project (+0.051).
Against the strong 16-category baseline it adds **nothing** (−0.001). The node
block's apparent survival gain is largely information that the finer cell-type
composition already contains. The enrichment block, by contrast, adds a small but
consistent amount under *both* baselines (+0.021 / +0.017) and beats the noise
control under both.

### 7.2 None of it is statistically established

**Every delta in the table is smaller than its own fold-to-fold standard
deviation** — e.g. `Baseline + Nodes` at lineage resolution is +0.051 with an SD of
±0.121 across the 10 folds. With 81 patients and 103 events, this endpoint cannot
resolve effects of this size. These are directional observations, not findings.

A second check points the same way. The optional separability test
(classifying `survival_status` from the spatial block, GroupKFold by patient) gives:

| feature set | AUC |
|---|---|
| Composition baseline | 0.646 |
| Baseline + spatial block | 0.703 |
| Baseline + noise | 0.675 |
| Spatial block alone | 0.633 |

That is a *pass* — the block beats both baseline and noise. But the version of this
same file committed before this run reported the exact opposite (baseline 0.603,
+spatial 0.492, +noise 0.628 — a clear fail). The current numbers reproduce exactly
across repeated runs on the current code and data; the provenance of the older file
could not be established. **An endpoint whose sign flips between code states is not
a basis for any claim**, and it is treated here as further evidence that the
outcome analysis on this cohort is underpowered rather than as a result in either
direction.

### 7.3 Not comparable to the legacy C = 0.733

The retired pipeline reported C = 0.733 for `Celltype + Enrichment + Markers(2)`.
That number is **not** comparable to this table: different cohort framing
(exclude-normal, 307 QC-filtered samples vs 308 here), a hand-selected 2-marker
subset rather than the full 8, and a different feature-construction path. It is
preserved in
[discarded/legacy_pipeline/RESULT_REPORT_legacy.md](../discarded/legacy_pipeline/RESULT_REPORT_legacy.md)
together with the code that produced it. Do not quote the two side by side.

---

## 8. What is established, and what is not

**Established, on two independent cohorts, without any outcome data:**

1. The **enrichment block** carries real, composition-independent, reproducible
   spatial signal — 5/5 STRONG at both taxonomies on both cohorts.
2. That conclusion is **robust** to every contested cell-type assignment (§6.1) and
   to the doublet exclusion (§6.2).
3. **6 of 8 node markers** measure a real, celltype-enriched, reproducible,
   non-abundance quantity on both cohorts — including both markers behind the
   project's best legacy configuration.

**Not established:**

4. **That either block predicts survival.** Every delta is within fold noise (§7.2),
   and one companion check reversed sign between code states. The honest statement
   is that this cohort is underpowered for the question, not that the answer is no.
5. **`tumor_mac_pdl1`** — premise contradicted on UPMC, flat on CRC.
6. **`cd8_granzymeb`** — works on UPMC, not on CRC. Not portable.

**Corrected documentation errors:**

7. Three legacy documents described the UPMC head & neck cohort as a "colorectal
   CODEX atlas (Schürch et al. 2020)". UPMC and CRC are different cohorts; conflating
   them would collapse the external-validation claim entirely, since CRC would then
   validate nothing.
8. The previous `RESULT_REPORT.md` contained two unresolved git merge conflicts.
   Resolved in the archived copy.

---

## 9. Limitations

- **Two cohorts, one of them without outcome data.** External validation covers the
  *measurement*, not the *prognosis*.
- **CRC survival exists but was not obtained.** Schürch et al. 2020 Supplementary
  Table S1 carries OS days and vital status for these 35 patients. Merging it in is
  the single highest-value next step: it would make CRC a genuine external outcome
  test. It will still be underpowered at 35 patients.
- **The registry is cited and ontology-grounded but not ontology-*verified*.** Every
  row carries `verified=no` pending manual confirmation against the Cell Ontology
  (procedure in [CELLTYPE_MAPPING.md](CELLTYPE_MAPPING.md)). `schema.CL_LINEAGE_ANCHOR`
  is a curated lookup of the terms this project uses, not a traversal of the OBO
  ancestry — it catches inconsistency, it does not prove a term is the right one.
- **The permutation null uses 20 draws.** Enough to separate z ≈ 5 from z ≈ 0;
  the third decimal of a z in the thousands is not meaningful.
- **CRC's `groups` (CLR/DII) was not used.** It is the cohort's only outcome-adjacent
  label, but it reaches the sample parquets and not `manifest.parquet`, so the
  separability check cannot read it. It is also partly circular for this feature set
  — CLR vs DII is *defined by* immune infiltration architecture, which is close to
  what three of the five enrichment features measure.

---

## 10. Reproducing every number

```bash
# Phase 1 — ingest all three cohorts from raw
cd data_preprocessing
python run_ingest.py --dataset UPMC
python run_ingest.py --dataset CRC
python run_ingest.py --dataset CRC_doublets

# Phase 2 — enrichment block (§4, §6)
cd ../spatial_positional_encoding
python run_verify.py --dataset UPMC --taxonomy lineage --perturb-map
python run_verify.py --dataset UPMC --taxonomy native
python run_verify.py --dataset CRC  --taxonomy lineage --perturb-map
python run_verify.py --dataset CRC  --taxonomy native
python run_verify.py --dataset CRC_doublets --taxonomy lineage
python run_verify.py --dataset UPMC --taxonomy lineage --label-col survival_status

# Phase 2b — node marker block (§5)
python run_verify_nodes.py --dataset UPMC --taxonomy lineage
python run_verify_nodes.py --dataset UPMC --taxonomy native
python run_verify_nodes.py --dataset CRC  --taxonomy lineage
python run_verify_nodes.py --dataset CRC  --taxonomy native
python run_verify_nodes.py --dataset CRC_doublets --taxonomy lineage

# Phase 3 — optional survival, cohorts that have it (§7)
python run_survival.py --dataset UPMC --taxonomy lineage
python run_survival.py --dataset UPMC --taxonomy native
python run_survival.py --dataset CRC              # exits cleanly: no survival
```

All outputs land in `data_preprocessing/datasets/<NAME>/processed/verification/`.

## 11. Companion documents

| document | contents |
|---|---|
| [PIPELINE_EXPLAINED.md](PIPELINE_EXPLAINED.md) | how the pipeline works, end to end, plus a step-by-step manual verification guide |
| [VERIFICATION.md](VERIFICATION.md) | the verification design and metric definitions |
| [CELLTYPE_MAPPING.md](CELLTYPE_MAPPING.md) | the registry, its grounding, and the confirmation procedure |
| [LITERATURE_REVIEW.md](LITERATURE_REVIEW.md) | cohort provenance and verified citations |
| [ENRICHMENT_FEATURES_REPORT.md](ENRICHMENT_FEATURES_REPORT.md), [KL_MEAN_REPORT.md](KL_MEAN_REPORT.md), [NODE_FEATURES_REPORT.md](NODE_FEATURES_REPORT.md) | legacy-era deep dives on individual feature families — design rationale still current, **numbers superseded by this report** |
| [discarded/legacy_pipeline/](../discarded/legacy_pipeline/) | the retired pipeline and its results, kept for traceability |
