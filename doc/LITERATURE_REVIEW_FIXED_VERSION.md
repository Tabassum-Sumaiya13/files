# Literature Review — Published Work Backing Each Pipeline Component

**Verification status.** Every DOI below was resolved against the **Crossref API**
and confirmed to be a real `journal-article` whose author, year, venue, volume,
issue, pages and title match what is cited here. Citations were *not* written from
memory. Where no strong canonical source exists, that is stated explicitly rather
than filled with a plausible-looking reference.

**Confidence key:**
`VERIFIED` = DOI resolved, metadata matches.
`GAP` = no canonical citation found; see the note.

---

## 1. CODEX / multiplexed antibody imaging (the imaging platform)

CODEX (CO-Detection by indEXing) uses DNA-barcoded antibodies revealed over
iterative imaging cycles, demonstrated at 28 markers in mouse spleen and
scalable to 56+ markers on human tissue. This is the instrument that produced
both cohorts' data — it is why each cell has ~40–56 marker intensities *and* an
`(X, Y)` position, which is the precondition for any spatial analysis.

- `VERIFIED` **Goltsev Y, Samusik N, Kennedy-Darling J, et al.** *Deep Profiling of
  Mouse Splenic Architecture with CODEX Multiplexed Imaging.* **Cell** 174(4):968–981.e15 (2018).
  DOI: [10.1016/j.cell.2018.07.010](https://doi.org/10.1016/j.cell.2018.07.010)
  — the original CODEX technology paper; 28-antibody panel, introduced indexed
  niche (i-niche) analysis using Delaunay neighbor composition (see §3, §10).
- `VERIFIED` **Black S, Phillips D, Hickey JW, et al.** *CODEX multiplexed tissue
  imaging with DNA-conjugated antibodies.* **Nature Protocols** 16(8):3802–3835 (2021).
  DOI: [10.1038/s41596-021-00556-8](https://doi.org/10.1038/s41596-021-00556-8)
  — step-by-step protocol; cite for methods detail.

## 2. The CRC cohort (this pipeline's second dataset)

The exact dataset ingested as `datasets/CRC/`: iterative CODEX on FFPE tissue
microarrays, **140 tissue regions from 35 advanced-stage CRC patients, 56 markers**,
split into the **CLR** (Crohn's-like reaction) vs **DII** (diffuse inflammatory
infiltration) patient groups. Those numbers match the ingest exactly (140 samples,
35 patients, 56 markers), confirming provenance.

- `VERIFIED` **Schürch CM, Bhate SS, Barlow GL, et al.** *Coordinated Cellular
  Neighborhoods Orchestrate Antitumoral Immunity at the Colorectal Cancer Invasive
  Front.* **Cell** 182(5):1341–1359.e19 (2020).
  DOI: [10.1016/j.cell.2020.07.005](https://doi.org/10.1016/j.cell.2020.07.005)

> ⚠️ **Citation trap.** Many sources (and the automated search) report this paper as
> DOI `10.1016/j.cell.2020.10.021`. Crossref shows that DOI is **Cell 183(3):838** — a
> **one-page erratum**, not the article. Cite `10.1016/j.cell.2020.07.005`
> (Cell 182(5):1341–1359.e19). Optionally cite the erratum separately.

## 3. Cellular neighborhoods / spatial niches

The idea that a cell's *local neighborhood composition* — the cell types among its
spatial neighbors — is itself a meaningful unit of tissue biology. Two independent
formulations anchor this concept:

- **Goltsev et al. 2018** — defined "indexed niches" (i-niches) using
  **first-tier Delaunay neighbors** (not kNN), then clustered the relative cell-type
  frequency in each cell's neighbor ring into **100 i-niche types**. This is the
  origin of neighbor-composition analysis and the direct conceptual ancestor of this
  pipeline's neighbor-count matrix M.
- **Schürch et al. 2020** — clustered kNN windows across 140 CRC tissue regions into
  **nine conserved cellular neighborhoods (CNs)**: T cell enriched, Bulk tumor,
  Tumor/immune mix, Macrophage enriched, TLS (follicle), Tumor boundary,
  Tumor/immune mix (variant), Smooth muscle, Granulocyte enriched. CN-5 (TLS) was
  the main distinguishing feature between CLR and DII groups.

- `VERIFIED` **Goltsev et al. 2018** (as above) — origin of i-niche / Delaunay
  neighbor composition.
- `VERIFIED` **Schürch et al. 2020** (as above) — origin of the 9-CN framework in
  tumors.

## 4. Permutation-based neighborhood enrichment (the `null_z` test)

The core validation method of this pipeline: hold the spatial graph fixed, randomly
**shuffle cell-type labels**, recompute, and score the real value against that null.
This is a long-established test in multiplexed tissue imaging — it is exactly what
`run_verify.py` does to produce `null_z_median`.

- `VERIFIED` **Schapiro D, Jackson HW, Raghuraman S, et al.** *histoCAT: analysis of
  cell phenotypes and interactions in multiplex image cytometry data.*
  **Nature Methods** 14(9):873–876 (2017).
  DOI: [10.1038/nmeth.4391](https://doi.org/10.1038/nmeth.4391)
  — the origin of the shuffle-labels neighborhood/interaction permutation test.
- `VERIFIED` **Palla G, Spitzer H, Klein M, et al.** *Squidpy: a scalable framework
  for spatial omics analysis.* **Nature Methods** 19(2):171–178 (2022).
  DOI: [10.1038/s41592-021-01358-2](https://doi.org/10.1038/s41592-021-01358-2)
  — modern implementation (`nhood_enrichment`): permutation z-score on a fixed graph.

**This is the strongest methodological precedent in the whole pipeline** — the
permutation null is not a bespoke invention, it is the field-standard test.

## 5. Delaunay triangulation as the cell-neighbor graph

The pipeline defines "neighbors" by Delaunay triangulation of cell centroids —
parameter-free (no `k` or radius to justify), with mean degree ≈ 6 by Euler's
formula. Delaunay is an established option for spatial-omics neighbor graphs,
offered alongside kNN in the standard toolkits.

- `VERIFIED` **Palla et al. 2022 (Squidpy)** — supports Delaunay spatial graph construction.
- `VERIFIED` **Dries R, Zhu Q, Dong R, et al.** *Giotto: a toolbox for integrative
  analysis and visualization of spatial expression data.* **Genome Biology** 22(1):78 (2021).
  DOI: [10.1186/s13059-021-02286-2](https://doi.org/10.1186/s13059-021-02286-2)
  — pairs Delaunay neighbor graphs with permutation-based proximity enrichment.

> `GAP` There is **no single canonical paper** for "Delaunay graphs for tissue cell
> neighbors." The geometry predates biology (Delaunay, 1934). Cite the *implementations*
> above, and justify the choice on its parameter-free property, not on a founding paper.

## 6. arcsinh marker normalization (cofactor 5)

Marker intensities span orders of magnitude and can go negative after background
subtraction, so `arcsinh(x / 5)` is the standard compression in mass/fluorescence
cytometry. This is `APPLY_ARCSINH` / `ARCSINH_COFACTOR = 5.0` in the adapter configs.

- `VERIFIED` **Bendall SC, Simonds EF, Qiu P, et al.** *Single-Cell Mass Cytometry of
  Differential Immune and Drug Responses Across a Human Hematopoietic Continuum.*
  **Science** 332(6030):687–696 (2011).
  DOI: [10.1126/science.1198704](https://doi.org/10.1126/science.1198704)
  — the landmark CyTOF paper that introduced the arcsinh transformation and
  cofactor-5 convention (Supplementary Fig. S2).

## 7. Random Survival Forests and the concordance index

Used only in the *optional, per-dataset* survival downstream (`run_survival.py`):
RSF handles right-censored outcomes; Harrell's C-index measures how well predicted
risk ranks actual event order (0.5 = coin flip).

- `VERIFIED` **Ishwaran H, Kogalur UB, Blackstone EH, Lauer MS.** *Random survival
  forests.* **The Annals of Applied Statistics** 2(3):841–860 (2008).
  DOI: [10.1214/08-AOAS169](https://doi.org/10.1214/08-AOAS169)
- `VERIFIED` **Harrell FE, Califf RM, Pryor DB, Lee KL, Rosati RA.** *Evaluating the
  Yield of Medical Tests.* **JAMA** 247(18):2543–2546 (1982).
  DOI: [10.1001/jama.1982.03320430047030](https://doi.org/10.1001/jama.1982.03320430047030)
  — origin of the C-index.
- `VERIFIED` **Harrell FE Jr, Lee KL, Mark DB.** *Multivariable prognostic models…*
  **Statistics in Medicine** 15(4):361–387 (1996).
  DOI: [10.1002/(SICI)1097-0258(19960229)15:4<361::AID-SIM168>3.0.CO;2-4](https://doi.org/10.1002/(SICI)1097-0258(19960229)15:4%3C361::AID-SIM168%3E3.0.CO;2-4)
  — the standard reference formalizing C-index for censored data.

## 8. Patient-grouped cross-validation (leakage control)

Both cohorts have **multiple tissue samples per patient** (UPMC: 308 samples / 81
patients; CRC: 140 / 35). Splitting a patient's samples across train and test folds
leaks information and inflates scores — hence `StratifiedGroupKFold(groups=patient_id)`.

- `GAP` There is **no canonical origin paper** for grouped CV; it is a standard
  technique (scikit-learn `GroupKFold`). Cite the software and a demonstration of
  the harm instead:
- `VERIFIED` **Rosenblatt M, Tejavibulya L, Jiang R, Noble S, Scheinost D.** *Data
  leakage inflates prediction performance in connectome-based machine learning models.*
  **Nature Communications** 15(1):1829 (2024).
  DOI: [10.1038/s41467-024-46150-w](https://doi.org/10.1038/s41467-024-46150-w)
  — empirically shows correlated within-subject samples split across folds inflate
  apparent performance, and that grouping by subject is required.

## 9. Spatial TME architecture predicts outcome (why this pipeline exists)

The premise: *spatial arrangement carries prognostic information beyond cell-type
composition.* This is the motivation for testing features against a composition
baseline rather than assuming spatial features are informative.

- `VERIFIED` **Jackson HW, Fischer JR, Zanotelli VRT, et al.** *The single-cell
  pathology landscape of breast cancer.* **Nature** 578(7796):615–620 (2020).
  DOI: [10.1038/s41586-019-1876-x](https://doi.org/10.1038/s41586-019-1876-x)
  — cell communities (spatially contiguous groups) predictive of outcome; used
  Shannon entropy and KL divergence for spatial heterogeneity quantification.
- `VERIFIED` **Danenberg E, Bardwell H, Zanotelli VRT, et al.** *Breast tumor
  microenvironment structures are associated with genomic features and clinical
  outcome.* **Nature Genetics** 54(5):660–669 (2022).
  DOI: [10.1038/s41588-022-01041-y](https://doi.org/10.1038/s41588-022-01041-y)
  — unsupervised discovery of 10 recurrent TME structures; "suppressed expansion"
  (Treg + dysfunctional T cells) independently prognostic on top of cell composition.
- `VERIFIED` **Schürch et al. 2020** — found PD-1+CD4+ T cell frequency within
  CN-9 (granulocyte-enriched neighborhood) significantly associated with overall
  survival in **DII patients only** (n=18, 13 deaths; Cox p=0.006). CLR patients
  (4 deaths) had insufficient events for survival analysis. Neither overall
  PD-1+CD4+ frequency nor CN-9 abundance alone was prognostic — only the
  cell-type-within-neighborhood combination.

### Marker shortlist provenance

The 8 functional markers (Ki67, PD-1, ICOS, GranzymeB, CD45RO, PD-L1, FoxP3, HLA-DR)
encode celltype-conditioned functional states in `spatial_positional_encoding/src/node_features.py`
(via Cell Ontology terms from `data_preprocessing/celltype_registry.csv`). These are
standard immune functional markers used broadly across the spatial proteomics literature.
Their selection for this project was **inspired by**:

- **Schürch 2020** — PD-1+ CD4 T cells as top survival hit; ICOS, GranzymeB, CD45RO
  as key functional markers
- **Jackson 2020** — functional-marker means (Ki67, PD-L1, GranzymeB) prognostic in
  breast cancer IMC
- **Danenberg 2022** — FoxP3 in Tregs as part of the "suppressed expansion" TME
  structure

The current implementation (`node_features.py`) resolves cell-type conditioning
per-cohort through Cell Ontology terms, replacing the original hard-coded version
(`discarded/legacy_pipeline/src/marker_states.py`).

## 10. KL divergence / abundance-corrected enrichment

The pipeline's distinctive readout: divide observed neighbor fractions by the
global composition (`E[i][j] = P[i][j] / p_j`, `1.0 = chance`) and summarise with
KL divergence (`kl_mean`, `kl_tumor`) and log2 block enrichments.

> `GAP` **No canonical paper was found for this exact formulation** (per-sample
> KL of neighbor distribution vs. global composition, collapsed to a single scalar).
> The `kl_mean` ratio is a weaker version of the established Squidpy/histoCAT
> permutation z-score — it is mean-only and discards the variance term, causing
> ~43% artifact in cell types with <10 cells (KL_MEAN_REPORT.md §6). A properly
> calibrated z-score (e.g., Squidpy's `nhood_enrichment` or the Analytical NES
> closed form, arXiv 2506.18692) would not have this limitation.

**Closest conceptual anchors (none identical):**

1. **Goltsev et al. 2018** — i-niche analysis using Delaunay neighbor composition;
   the identical starting point as this pipeline's neighbor-count matrix M.
2. **Squidpy `nhood_enrichment` (Palla 2022)** — permutation z-score on a fixed
   graph; same conceptual goal (abundance-corrected neighborhood enrichment) with
   variance calibration.
3. **histoCAT (Schapiro 2017)** — origin of the shuffle-labels permutation test
   within a single image.
4. **Jackson et al. 2020** — used Shannon entropy and KL divergence for intratumor
   heterogeneity in breast cancer IMC; closest precedent for information-theoretic
   divergence on spatial omics distributions.
5. **Schrom E, et al.** *Spatial Patterning Analysis of Cellular Ensembles (SPACE)
   finds complex spatial organization at the cell and tissue levels.* **PNAS** 122(6) (2025).
   DOI: [10.1073/pnas.2412146122](https://doi.org/10.1073/pnas.2412146122)
   — uses KL divergence on co-occurrence distributions at the group-comparison
   level (transMI: pairwise KL between specimens → network modularity → Z-score).
   Technically different from our per-sample `kl_mean` (which is intra-specimen
   neighbor vs. composition), but supports the broader use of KL in spatial omics.

**What is genuinely novel in this project (not the ratio itself):**

- **Collapsing** the full cell-type-pair enrichment matrix into a single parsimonious
  scalar while preserving composition-independent signal
- A **dataset-agnostic verification framework** (3 survival-independent tests:
  `null_z`, `1-R²`, `stability_r`) that proves spatial signal is real and
  composition-independent on any cohort without requiring survival data
- **Cross-cohort portability** via Cell Ontology-grounded cell-type registry
  (`data_preprocessing/celltype_registry.csv`)

---

## Outstanding gap you must fill

`GAP` **The UPMC head-and-neck cohort's own source paper was not identified.** This
review covered the CRC dataset (Schürch et al.) but the UPMC cohort — the project's
baseline, referenced internally as the "spatsurv" baseline — has a source
publication that this search did not surface. **You need to supply that citation**;
it should be cited wherever UPMC results are reported.

## Summary table

| # | Component | Status | Primary citation |
|---|---|---|---|
| 1 | CODEX imaging | VERIFIED | Goltsev 2018 *Cell*; Black 2021 *Nat Protoc* |
| 2 | CRC dataset | VERIFIED | Schürch 2020 *Cell* 182(5) |
| 3 | Cellular neighborhoods | VERIFIED | Goltsev 2018 (i-niches); Schürch 2020 (9 CNs) |
| 4 | Permutation null test | VERIFIED | Schapiro 2017 *Nat Methods*; Palla 2022 *Nat Methods* |
| 5 | Delaunay neighbor graph | VERIFIED (impl.) | Palla 2022; Dries 2021 *Genome Biol* |
| 6 | arcsinh cofactor 5 | VERIFIED | Bendall 2011 *Science* |
| 7 | RSF + C-index | VERIFIED | Ishwaran 2008 *AOAS*; Harrell 1982/1996 |
| 8 | Grouped CV | GAP + support | Rosenblatt 2024 *Nat Commun* |
| 9 | Spatial TME → outcome | VERIFIED | Jackson 2020 *Nature*; Danenberg 2022 *Nat Genet*; Schürch 2020 (DII-only) |
| 10 | KL / abundance-corrected | **GAP** | Goltsev 2018 (method ancestor); Squidpy/histoCAT (permutation z-score); Jackson 2020 (KL in spatial); Schrom 2025 *PNAS* (supporting) |
| — | UPMC cohort source | **GAP** | *you must supply* |
