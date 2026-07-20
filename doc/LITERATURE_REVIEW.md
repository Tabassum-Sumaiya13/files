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
iterative imaging cycles, allowing ~50+ protein markers to be measured on a single
intact tissue section. This is the instrument that produced both cohorts' data —
it is why each cell has ~40–56 marker intensities *and* an `(X, Y)` position, which
is the precondition for any spatial analysis.

- `VERIFIED` **Goltsev Y, Samusik N, Kennedy-Darling J, et al.** *Deep Profiling of
  Mouse Splenic Architecture with CODEX Multiplexed Imaging.* **Cell** 174(4):968–981.e15 (2018).
  DOI: [10.1016/j.cell.2018.07.010](https://doi.org/10.1016/j.cell.2018.07.010)
  — the original CODEX technology paper.
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
k nearest spatial neighbors — is itself a meaningful unit of tissue biology.
Schürch et al. clustered these kNN windows into nine conserved "cellular
neighborhoods" (CNs). This is the direct conceptual ancestor of this pipeline's
neighbor-count matrix `M`, from which all 5 enrichment features are derived.

- `VERIFIED` **Schürch et al. 2020** (as above) — origin of the CN concept in tumors.
- `VERIFIED` **Goltsev et al. 2018** (as above) — introduces neighborhood-level
  analysis of tissue architecture.

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
  — the landmark CyTOF paper that established the arcsinh/cofactor-5 convention.

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
- `VERIFIED` **Danenberg E, Bardwell H, Zanotelli VRT, et al.** *Breast tumor
  microenvironment structures are associated with genomic features and clinical
  outcome.* **Nature Genetics** 54(5):660–669 (2022).
  DOI: [10.1038/s41588-022-01041-y](https://doi.org/10.1038/s41588-022-01041-y)
- `VERIFIED` **Schürch et al. 2020** — CN organization linked to CRC survival.

These three are also the source of the functional markers hard-coded in
`src/marker_states.py` (Ki67, PD-1, ICOS, GranzymeB, CD45RO, PD-L1, FoxP3, HLA-DR).

## 10. KL divergence / abundance-corrected enrichment

The pipeline's distinctive readout: divide observed neighbor fractions by the
global composition (`E = P[i][j] / p_j`, `0 = chance`) and summarise with
KL divergence (`kl_mean`, `kl_tumor`) and log2 block enrichments.

> `GAP` **No canonical paper was found for this exact formulation.** The closest
> verified anchor is:
- `VERIFIED` **Schrom E, et al.** *Spatial Patterning Analysis of Cellular Ensembles
  (SPACE) finds complex spatial organization at the cell and tissue levels.*
  **PNAS** 122(6) (2025).
  DOI: [10.1073/pnas.2412146122](https://doi.org/10.1073/pnas.2412146122)
  — uses KL divergence on co-occurrence distributions with randomized nulls that
  control for compositional abundance.

**Read this as an opportunity, not a weakness.** Components 1–9 are all standard and
well-cited; the abundance-corrected KL enrichment readout — and the
survival-independent verification matrix built on it — is where this project's
**novel contribution** sits. Frame it as: *standard permutation-null methodology
(histoCAT/squidpy) applied to a new abundance-corrected readout, validated
dataset-agnostically.*

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
| 3 | Cellular neighborhoods | VERIFIED | Schürch 2020; Goltsev 2018 |
| 4 | Permutation null test | VERIFIED | Schapiro 2017 *Nat Methods*; Palla 2022 *Nat Methods* |
| 5 | Delaunay neighbor graph | VERIFIED (impl.) | Palla 2022; Dries 2021 *Genome Biol* |
| 6 | arcsinh cofactor 5 | VERIFIED | Bendall 2011 *Science* |
| 7 | RSF + C-index | VERIFIED | Ishwaran 2008 *AOAS*; Harrell 1982/1996 |
| 8 | Grouped CV | GAP + support | Rosenblatt 2024 *Nat Commun* |
| 9 | Spatial TME → outcome | VERIFIED | Jackson 2020 *Nature*; Danenberg 2022 *Nat Genet* |
| 10 | KL / abundance-corrected | **GAP** | Schrom 2025 *PNAS* (closest anchor) |
| — | UPMC cohort source | **GAP** | *you must supply* |
