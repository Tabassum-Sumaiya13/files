# Enrichment Features Report — Abundance-Corrected Spatial Scalars

*Design rationale for the enrichment block. Cohort: **UPMC head & neck (HNSCC)**. Companion to [KL_MEAN_REPORT.md](KL_MEAN_REPORT.md) (deep dive on the single strongest feature).*

> **Two corrections to this document, 2026-07-20.**
>
> 1. **Cohort.** Earlier versions described this data as a "colorectal CODEX atlas".
>    It is the **UPMC head & neck** cohort. The colorectal (Schürch et al. 2020)
>    data is the separate **CRC** cohort used for external validation.
> 2. **Source and numbers are superseded.** `src/enrichment_features.py` has been
>    retired to [../discarded/legacy_pipeline/](../discarded/legacy_pipeline/) —
>    it was hardwired to UPMC's 16 integer cluster ids and prebuilt edge files. The
>    same five scalars are now computed by
>    [`src/spatial_features.py`](../spatial_positional_encoding/src/spatial_features.py)
>    on any cohort. Current numbers: [RESULT_REPORT.md](RESULT_REPORT.md). The
>    **feature definitions and reasoning below are still current.***

---

## 1. TL;DR

- **What they are:** 5 scalar features per sample that measure spatial organisation *after removing what cell abundance alone explains*. The base paper flattens a 16×16 neighbour matrix into **256** columns; these collapse the same matrix into **5** abundance-corrected numbers.
- **The key move:** divide observed neighbourhood fractions by the global mix — `E[i][j] = P[i][j] / p_j`. `E = 1` is exactly chance; the celltype baseline already knows `p`, so this deletes the redundant part and keeps only the spatial signal.
- **Headline finding:** the 5-scalar block recovers **almost all** of the 256-block's survival performance (**C = 0.716 vs 0.718**) using **5 features instead of 256**, and cleanly beats its width-matched noise control (+0.028 [+0.005, +0.053]).
- **Caveat:** the 5 features are **highly redundant with each other** (r up to 0.80), and within the block **only `kl_mean` reliably adds signal** — the other four are near-zero or hurt. Alone (without celltype composition) the block is *worse* than baseline (0.580).

---

## 2. Why these features exist — the confound they remove

The base paper's readout is `P[i][j]` = "of type i's neighbours, what fraction are type j", flattened to 256 columns. Under a **random** spatial arrangement, every row of `P` collapses to the global composition:

```
P[i][j]  ->  p_j        (p = global proportion of each celltype)
```

The celltype-proportion baseline **already hands the model `p`**. So most of those 256 columns re-encode information the baseline has — which is why `Celltype + delaunay(256)` gains only a little over `Celltype` alone: it is largely the same information twice, spread across 256 noisy columns.

Dividing by `p` removes the shared part:

```
E[i][j] = P[i][j] / p_j       E = 1 -> exactly chance
                              E > 1 -> genuinely enriched
                              E < 1 -> genuinely depleted
```

**Every enrichment feature is 0 under the random-mix null, and none is computable from abundance alone.** 5 features against ~27 patient events is EPV ≈ 5 — the same order as the >10 rule of thumb, versus the 272-wide blocks that sit ~90× off it.

Graph = **Delaunay** (Experiment A winner, parameter-free).

---

## 3. The 5 features

| Feature | Definition | What it measures |
|---|---|---|
| `kl_mean` | mean over celltypes of KL(neighbourhood ‖ global mix), in bits | how organised the whole tissue is (0 = random soup) |
| `kl_tumor` | KL(pooled-tumour neighbourhood ‖ global mix) | how distinctive the tumour niche is |
| `self_enrich` | mean log₂(P[i][i] / p_i) | do celltypes clump with their own kind beyond chance |
| `immune_tumor` | log₂(P[immune→tumor] / p_tumor) | immune infiltration of tumour beyond chance |
| `stroma_tumor` | log₂(P[stroma→tumor] / p_tumor) | stromal/vessel–tumour interface beyond chance |

Celltype groups: **IMMUNE** = APC, B, CD4 T, CD8 T, Granulocyte, Macrophage, Naive immune. **TUMOR** = the 6 tumour variants. **STROMA** = Lymph vessel, Stromal/Fibroblast, Vessel.

Numerical guards: a half-count floor (0.5) keeps `log2(0)` finite; rows with zero neighbours are skipped; every feature returns 0 for an empty sample.

---

## 4. Characteristics — distribution and redundancy

### 4.1 Distribution across 307 samples

| Feature | mean | std | min | max | reads as |
|---|---|---|---|---|---|
| kl_mean | 0.552 | 0.181 | 0.121 | 1.260 | organised (>0) |
| kl_tumor | 0.240 | 0.186 | 0.008 | 1.117 | tumour niche mildly distinctive |
| self_enrich | **2.212** | 0.345 | 1.232 | 3.123 | strong self-clustering (~4.6× chance) |
| immune_tumor | **−0.767** | 0.464 | −2.438 | −0.030 | immune **depleted** around tumour |
| stroma_tumor | **−1.225** | 0.662 | −3.151 | 0.151 | stroma strongly **depleted** around tumour |

**Biological read-through:** `self_enrich` is large and positive — cells overwhelmingly sit with their own kind. `immune_tumor` and `stroma_tumor` are **negative in nearly every sample** — immune and stromal cells are *excluded* from the tumour compartment (log₂ < 0 = below chance). This is the classic "immune-excluded / cold tumour" architecture, quantified.

### 4.2 They are highly correlated with each other

```
              kl_mean  kl_tumor  self_enrich  immune_tumor  stroma_tumor
kl_mean          1.00      0.44         0.80         -0.65         -0.66
kl_tumor         0.44      1.00         0.36         -0.77         -0.36
self_enrich      0.80      0.36         1.00         -0.48         -0.47
immune_tumor    -0.65     -0.77        -0.48          1.00          0.56
stroma_tumor    -0.66     -0.36        -0.47          0.56          1.00
```

`kl_mean` ↔ `self_enrich` at **r = 0.80** (both dominated by self-clustering — the KL matrix is almost entirely diagonal). `kl_tumor` ↔ `immune_tumor` at **r = −0.77**. So the block has **far fewer than 5 independent dimensions** — a distinctive tumour niche *is* an immune-excluded one. This is why adding all five buys little over the best one.

---

## 5. Findings — does the block predict survival?

Protocol: RandomSurvivalForest, patient-grouped StratifiedGroupKFold, 20 seeds; `dC vs NOISE` = paired gain vs a same-width random block (the honest, dilution-cancelling test). Cohort: 307 samples, 81 patients, 103 events.

### 5.1 The block vs the 256-column readout (Experiment C)

| Feature set | Feats | C-index | ΔC vs NOISE | Verdict |
|---|---|---|---|---|
| Celltype + delaunay(256) | 272 | 0.718 | +0.062 [+0.030, +0.099] | REAL SIGNAL |
| **Celltype + Enrichment** | **21** | **0.716** | +0.028 [+0.005, +0.053] | **REAL SIGNAL** |
| Celltype Proportions (baseline) | 16 | 0.690 | — | baseline |
| Enrichment alone | 5 | 0.580 | — | **worse than baseline** |

**Two conclusions:**

1. **Efficiency.** The 5-scalar block reaches 0.716 vs the 256-block's 0.718 — it captures essentially all the recoverable spatial signal with **~2% of the columns**, and beats its own noise control cleanly. As a *compact readout*, abundance-corrected enrichment works.
2. **Complementary, not standalone.** `Enrichment alone` = 0.580, well below baseline. Spatial organisation without knowing *which* cells are present is not enough — the features only help *on top of* composition.

### 5.2 Within the block, one feature does the work

From the ablation in [KL_MEAN_REPORT.md](KL_MEAN_REPORT.md) (exclude-normal cohort, 27 events):

- Each feature alone, added to celltype: **`kl_mean` +0.021 [+0.003, +0.038]** is the only one whose CI clears 0. `immune_tumor`, `kl_tumor`, `stroma_tumor` ≈ +0.006–0.008 (CIs cross 0). `self_enrich` −0.006 (hurts).
- **Leave-one-out** from the full 5: dropping `kl_mean` costs −0.009; dropping *any other* feature **improves** the model (+0.001 to +0.005). The top 16 of 31 subsets all contain `kl_mean`; the best subset without it ranks 17th.

So the 5-block's usable signal is concentrated in **`kl_mean`**, with the rest largely redundant (§4.2 explains why).

### 5.3 Redundant on top of the full 256

`Celltype + delaunay(256) + kl_mean` = 0.7188 vs `Celltype + delaunay(256)` = 0.7186 → paired ΔC = **+0.0003 [−0.023, +0.016]**, i.e. zero. This is expected: every enrichment scalar is a *deterministic function* of the same 16×16 matrix the 256 columns come from. There is no new information — only a question of whether the forest can derive the summary itself. It can. **Pick one:** the 256 for the top number, the enrichment scalars for parsimony/interpretability.

---

## 6. Per-feature verdicts

| Feature | Detects real structure | Adds to model | Biology | Verdict |
|---|---|---|---|---|
| `kl_mean` | PASS (15× noise floor, see KL report) | **+0.021 — the workhorse** | tracks tissue type & recurrence | **keep** |
| `kl_tumor` | partial | ≈0 | tumour-niche distinctiveness | redundant with immune_tumor |
| `self_enrich` | PASS but ≈ kl_mean (r 0.80) | −0.006 (hurts) | self-clustering | drop (subsumed) |
| `immune_tumor` | weak | ≈0 | immune exclusion (good biology) | drop — best story, no signal |
| `stroma_tumor` | weak | ≈0 | desmoplastic barrier | drop |

The irony worth noting: `immune_tumor` had the strongest literature backing (immune desert/excluded/inflamed phenotypes) and contributes nothing here; `kl_mean`, invented from a statistical argument, is the one that works.

---

## 7. Limitations

1. **Redundant with the 256** — adds exactly zero on top of the full block.
2. **Mean-only enrichment** — this is a ratio, not a permutation z-score (unlike Squidpy `nhood_enrichment` / histoCAT). It discards the variance term, so 2× enrichment from 10,000 edges and from 12 edges look identical. Rare celltypes inject artifact (43% below 10 cells for `kl_mean`).
3. **High internal correlation** — effectively 1–2 independent dimensions, not 5.
4. **Standalone-weak** — only meaningful added to celltype composition.
5. **Overall survival** — `kl_mean` is null for OS in this cohort (p = 0.80) but significant for *recurrence* (p = 0.005–0.010); the enrichment signal is better chased through recurrence than OS (see KL report §9, §14).

---

## 8. Bottom line

> **The enrichment features are an efficient, honest, abundance-corrected readout: 5 scalars recover nearly all of the 256-block's survival signal (0.716 vs 0.718) and beat noise. But they are highly redundant with each other and with the 256, and the usable signal is concentrated almost entirely in `kl_mean`.**

**Recommendation:** use the enrichment block when **parsimony/interpretability** matters (21 features, high EPV, clean noise-beating verdict); use the full 256 for the **best raw number**; never stack both (redundant). If reducing the block, keep **`kl_mean`** and consider dropping the other four.

```bash
python spatial_positional_encoding/run_survival.py --experiment C   # enrichment(5) vs delaunay(256)
```

*Related: [KL_MEAN_REPORT.md](KL_MEAN_REPORT.md) (single-feature deep dive), [NODE_FEATURES_REPORT.md](spatial_positional_encoding/NODE_FEATURES_REPORT.md) (functional-marker node features).*
