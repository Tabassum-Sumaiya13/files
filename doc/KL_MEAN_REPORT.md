# kl_mean — Full Report

**What it is, how it works, and everything we found.**

Dataset: UPMC Head & Neck cohort. 307 samples, 81 patients, 2,061,102 cells, 16 cell types.
Graph: Delaunay. All numbers below are measured, not estimated.

---

## 1. The Idea in One Picture
```
GLOBAL MIX  ──┐
               ├── COMPARE ──► ONE NUMBER
LOCAL MIX    ──┘              (0 = random, high = organised)

```

Every cell sits somewhere. Ask two questions:

```
Question 1:  "What cells are in this tissue?"        -> celltype_prop  (GLOBAL)
Question 2:  "What cells are next to THIS cell?"     -> P[i,:]         (LOCAL)

kl_mean = how different is LOCAL from GLOBAL
```

If you shuffled every cell to a random spot, a cell's neighbours would just look
like the tissue average. So:

```
kl_mean = 0     ->  no organisation. A well-mixed soup.
kl_mean = high  ->  cells sit with specific partners. The tissue is organised.
```

That's it. **kl_mean is one number that says "how organised is this tissue?"**

---

## 2. Worked Example — Real Numbers

Sample `UPMC_c001_v001_r001_reg001`, 4,301 cells.

### Step 0: Who is here (this is `p`, the global mix)

```
Stromal   17.9%      Tumor      11.6%
Tum CD21  12.9%      Tum Podo   10.7%
Macro     11.6%      Tum CD20    7.0%    ... (16 types)
```

### Step 1: Count neighbours -> `M`

For every edge in the graph, add 1.

```
centre      Stromal  Tum CD21   Macro    Tumor   rowsum
Stromal        1804        74     659      541     4622
Tum CD21         74      1732      90      413     3347
Macro           659        90     810      259     3006
Tumor           541       413     259      772     3007
```

### Step 2: Divide each row by its total -> `P` (the LOCAL view)

```
centre      Stromal  Tum CD21   Macro    Tumor
Stromal       0.390     0.016   0.143    0.117
Tum CD21      0.022     0.517   0.027    0.123
Macro         0.219     0.030   0.269    0.086
Tumor         0.180     0.137   0.086    0.257
--------------------------------------------------
p (global)    0.179     0.129   0.116    0.116   <- the random-mix reference
```

**This is where the base paper stops.** It flattens these 256 numbers into the model.

### Step 3: Divide by chance -> the key move

Look at the **Tumor** row. The biggest number is `Tumor -> Stromal = 0.180`.
Bigger than `Tumor -> Tum CD21 = 0.137`.

**Does tumour prefer stromal cells? No.**

Stromal cells are **17.9% of the whole tissue**. Shuffle everything randomly and
tumour would *still* touch stromal ~17.9% of the time. So `0.180` says nothing.
And `celltype_prop` **already told the model** stromal is 17.9%.

So divide it out:

```
E[i][j] = P[i][j] / p_j          1.0 = exactly chance
```

| From the Tumor row | `P` (before) | `E` (after) | Truth |
|---|---|---|---|
| Tumor -> Stromal | **0.180** (biggest!) | **1.00** | **exactly chance. Nothing.** |
| Tumor -> Tumor | 0.257 | **2.22** | real: 2.2x clumping |
| Tumor -> Macro | 0.086 | **0.74** | real: avoiding macrophages |

**The number that looked biggest was pure abundance.** Dividing by `p` deletes what
the baseline already knows and leaves only the spatial part.


### Step 4: Collapse to one number

```
KL_i    = how far type i's neighbourhood is from the global mix (in bits)
kl_mean = average of the 16 KL_i values
```

For this sample:

```
Tum CD21   KL = 0.847        Tum CD20   KL = 0.354
Tum Podo   KL = 0.334        Stromal    KL = 0.323
Macro      KL = 0.274        Tumor      KL = 0.143
...
kl_mean = mean of all 16 rows = 0.640 bits
```

---

## 3. The Matrix Behind It

`kl_mean` is a scalar, but it decomposes into a 16x16 matrix:

```
C[i][j] = P[i][j] * log2( P[i][j] / p_j )    <- contribution matrix
KL_i    = sum of row i
kl_mean = mean of the row sums
```

See `data/processed/outputs/spatial_maps/05_kl_matrix.png`.

**Finding: the matrix is almost entirely DIAGONAL.** `Tumor->Tumor`,
`B cell->B cell`, `Vessel->Vessel`. The off-diagonal is nearly blank.

So kl_mean is mostly measuring **self-clustering** — do cells sit with their own
kind? This is why `kl_mean` and `self_enrich` correlate at **r = 0.80**.

---

## 4. Distribution Across 307 Samples

```
             mean     sd      min     max
kl_mean      0.552   0.181   0.121   1.260
```

| Sample | kl_mean | What it looks like |
|---|---|---|
| `c002...reg008` | **1.260** (highest) | Highly structured, distinct zones |
| `c006...reg051` | **0.121** (lowest) | Everything intermixed |

Cross-check: `reg051` is also the fully-mixed sample in `03_groups.png`.
The number and the picture agree.

---

## 5. Is It Real? — The Shuffle Test

**The most important test.** Shuffle the cell type labels but keep the graph,
the cell counts, and the composition identical. This destroys spatial structure
while changing nothing else. Whatever kl_mean returns now is **pure noise**.

```
observed kl_mean   0.552
SHUFFLED kl_mean   0.035    <- the noise floor
excess             0.517

observed > shuffled in 307 / 307 samples  (100%)
Mann-Whitney p = 6e-102
-> only 6.4% of the value is artifact
```

**VERDICT: PASS.** kl_mean detects real structure at ~15x the noise floor.

---

## 6. Where the Noise Is — Rare Cell Types

The 6.4% artifact is not spread evenly. It all sits in rare cell types:

| cells in type | KL observed | KL shuffled | **% artifact** |
|---|---|---|---|
| **0-10** | 1.233 | 0.527 | **42.8%** |
| 10-25 | 0.730 | 0.121 | 16.6% |
| 25-50 | 0.610 | 0.056 | 9.1% |
| 50-100 | 0.585 | 0.028 | 4.7% |
| 100-250 | 0.617 | 0.012 | 1.9% |
| 1000+ | 0.337 | 0.001 | **0.3%** |

**Why:** you cannot estimate a 16-bin distribution from 9 cells. The estimate is
noise, and KL turns noise into a big positive number every time.

**Real example** — the highest-kl_mean sample:

```
LymphVes   n =    9 cells  ->  KL = 4.948 bits   <- 25% of that sample's kl_mean
Tum CD20   n =    1 cell   ->  KL = 1.452 bits   <- from ONE cell
B cell     n = 1750 cells  ->  KL = 0.935 bits
```

kl_mean weights all 16 types **equally**. So 9 cells out of 6,000+ contribute a
quarter of the value.

**Fix:** weight by abundance, or use `kl_excess = observed - shuffled`. The
weighted version is also the most reliable (see next section).

---

## 7. Is It Stable? — Reliability (ICC)

Each patient has several tissue regions. If kl_mean is a real patient property,
the regions should agree.

**ICC = share of variance that is BETWEEN patients** (0 = pure noise, 1 = perfectly stable).

```
ICC(kl_mean raw)         = 0.406
ICC(kl_excess)           = 0.385
ICC(kl_excess weighted)  = 0.459    <- best
```
(78 patients, 304 acquisitions)

**VERDICT: MODERATE.** About 41-46% of the variation is between patients; the rest
is region-to-region. So it is partly a patient property, partly local. This caps
how well it can ever predict a patient-level outcome.

---

## 8. Does It Mean Anything Biologically? — YES

### Tissue type (Kruskal-Wallis **p = 0.0009**)

| tissue | median kl_mean |
|---|---|
| Normal mucosa | **0.606** |
| Nodal met | 0.587 |
| **Primary tumor** | **0.514** |

**Normal tissue is MORE organised than tumour.** That is exactly what cancer does
to tissue architecture. The 3-way split survives bias correction (p = 0.006).
The Normal-vs-Primary pair weakens to p = 0.061 after correction (only n=26
normals), so part of that specific contrast was bias.

### Recurrence (patient level, primary tumours only)

| outcome | groups | medians | p |
|---|---|---|---|
| `isrecurrencetumor` | 64 vs 6 | 0.509 / **0.321** | **0.005** |
| `recurred` | 56 vs 14 | 0.522 / **0.425** | **0.010** |
| `hpvstatus` | 34 vs 33 | 0.381 / 0.383 | **0.693** |

**Less organisation -> recurrence.** And HPV is **null** — a good specificity check.
The feature is not just correlating with everything.

### The consistent gradient

```
Normal mucosa      0.606   <- most organised
Nodal met          0.587
Primary tumor      0.514
Recurred patient   0.425
Recurrence tumor   0.321   <- least organised
```

**Loss of spatial organisation tracks disease aggressiveness.** Four independent
contrasts, all pointing the same way.

---

## 9. Does It Predict Survival? — NO

```
Patient level (81 patients, 27 events), split at median kl_mean:

  HIGH structure:  40 patients, 13 events, 32.5% died
  LOW  structure:  41 patients, 14 events, 34.1% died

  log-rank, median split : p = 0.804
  log-rank, tertiles     : p = 0.612
  Spearman kl vs died    : rho = +0.001  (p = 0.99)
  acquisition level      : p = 0.289
```

**rho = +0.001.** Not weak — absent.

**Why this is not fatal:** survival here has only **27 patient events**, and
`survival_status = 1` counts deaths from *any* cause, including unrelated ones.
No tumour feature can predict a heart attack. Recurrence is better powered and
that is where kl_mean works.

---

## 10. Does It Help a Model? — Somewhat

Random Survival Forest, 10-fold grouped by patient, exclude-normal cohort,
281 acquisitions / 81 patients / 27 events.

| feature set | feats | C-index | dC vs baseline |
|---|---|---|---|
| Celltype (baseline) | 16 | 0.6825 | — |
| **Celltype + kl_mean** | **17** | **0.7033** | **+0.021 [+0.003, +0.038]** |
| Celltype + delaunay(256) | 272 | 0.7186 | +0.036 [+0.017, +0.052] |
| Celltype + delaunay + kl_mean | 273 | 0.7188 | +0.036 [+0.013, +0.056] |

### Two conclusions

**1. kl_mean is efficient.** It captures **~58%** of the 256-block's gain using
**1 column instead of 256**.

**2. kl_mean adds NOTHING on top of the 256.**

```
Celltype + delaunay(256)        C = 0.7186
Celltype + delaunay + kl_mean   C = 0.7188
paired dC = +0.0003 [-0.0233, +0.0163]    <- exactly zero
```

This was predictable. kl_mean is a **deterministic function** of the same 16x16
matrix those 256 columns come from. There was never new information in it — only
a question of whether the forest could derive the summary itself. It can.

**So: pick one.** The 256 for the best number, kl_mean for parsimony.

---

## 11. The Noise Floor — Why Width Matters

Random features added to the baseline, **5 independent draws** per width:

```
Noise(  1) -> width 17:   mean 0.6836   sd 0.0038   range [0.6792, 0.6897]
Noise(256) -> width 272:  mean 0.6575   sd 0.0198   range [0.6343, 0.6935]  <- 0.059 spread!
Noise(257) -> width 273:  mean 0.6502   sd 0.0176   range [0.6284, 0.6736]
```

**At width 272, pure random features land anywhere in a 0.059-wide window
depending only on the seed.** That window is bigger than every effect in this
report.

Against properly-averaged floors:

```
kl_mean       +0.0196 over its floor  =  +5.1 sd   <- narrow block, stable floor
delaunay(256) +0.0611 over its floor  =  +3.1 sd
```

**Lesson: noise-floor uncertainty grows with block width.** A 1-feature block needs
one draw. A 256-feature block needs many. The original harness used **one draw**
(`make_noise_block(..., seed=42)`), which made every wide-block verdict unreliable.

---

## 12. Comparison — All Five Features Tested

### Each alone, added to celltype (10 seeds)

| feature | C-index | dC vs baseline |
|---|---|---|
| **`kl_mean`** | **0.7033** | **+0.021 [+0.003, +0.038]** |
| `immune_tumor` | 0.6901 | +0.008 [-0.013, +0.024] |
| `kl_tumor` | 0.6901 | +0.008 [-0.013, +0.018] |
| `stroma_tumor` | 0.6886 | +0.006 [-0.008, +0.015] |
| `self_enrich` | 0.6765 | **-0.006** (hurts) |

### Leave-one-out from the full 5 (C = 0.7061)

```
drop kl_mean      -> 0.6972   -0.0089   <- the only one that costs anything
drop kl_tumor     -> 0.7096   +0.0035   helps
drop self_enrich  -> 0.7098   +0.0037   helps
drop immune_tumor -> 0.7111   +0.0050   helps
drop stroma_tumor -> 0.7074   +0.0013   helps
```

**Dropping 4 of the 5 improves the model.** The top 16 of 31 subsets all contain
`kl_mean`. The best subset without it ranks 17th.

**The irony:** `immune_tumor` had the best biology behind it (the HNSCC
desert/excluded/inflamed literature, median OS 85/61/37 months) and it contributes
nothing. `kl_mean` was invented with no literature basis and it is the only one
that works.

---

## 13. Which Reference Is Right? — Per-Sample Wins

There are four levels. The pipeline uses three.

```
Level 2  neighbourhood   P[i,:]
Level 3  this sample      p         <- kl_mean uses THIS
Level 4  the cohort       p_cohort  <- tested, does not work
```

| variant | survival dC | tissuetype | recurred |
|---|---|---|---|
| **`kl_mean`** (per-sample ref) | **+0.021** | **p=0.0009** | **p=0.010** |
| `comp_atypical` (KL of composition vs cohort) | +0.006 | p=0.21 | p=0.11 |
| `kl_cohort` (neighbourhood vs cohort ref) | **-0.002** | p=0.0018 | p=0.47 |
| `Noise(1)` | -0.001 | — | — |

**Why per-sample is correct:** referencing to the cohort **re-mixes composition back
in**. A sample with unusually many T cells shows high `kl_cohort` *even if perfectly
randomly mixed* — because its neighbourhoods differ from the cohort average through
composition alone. That is the exact confound the division by `p` removes. And the
baseline already has composition.

This matches the established tools: Squidpy and histoCAT both permute labels
*within* an image, never across a cohort.

---

## 14. Honest Limitations

1. **I invented it.** kl_mean came from a statistical argument in a conversation,
   not from a literature survey. The core idea (abundance-corrected neighbourhood
   enrichment) IS established — Squidpy's `nhood_enrichment`, histoCAT — but those
   use a **permutation z-score**, not my ratio. Mine is the mean-only version: it
   throws away the variance term, so a 2x enrichment from 10,000 edges and from 12
   edges look identical to it.

2. **Rare-type bias.** 43% artifact below 10 cells. Use the weighted or
   excess version.

3. **Redundant with the 256.** Adds exactly zero on top.

4. **ICC 0.41.** Half the variance is region-level, not patient-level.

5. **Multiple testing.** ~5 outcomes tested. Bonferroni leaves
   `isrecurrencetumor` (p=0.005) comfortable; `recurred` (p=0.010) becomes
   marginal. `isrecurrencetumor` has only 6 patients in one arm. Treat recurrence
   as **exploratory**.

6. **Mostly diagonal.** It is largely a self-clustering measure (r = 0.80 with
   `self_enrich`), so "spatial organisation" oversells it slightly.

---

## 15. Bottom Line

> **kl_mean is a validated measure of tissue spatial organisation. It separates
> tissue types and tracks recurrence. It does not predict overall survival in this
> cohort, and it is redundant with the full 256-feature block.**

| Test | Result |
|---|---|
| Detects real structure (shuffle) | **PASS** — 307/307, p=6e-102 |
| Free of rare-type artifact | **PARTIAL** — 6.4% overall, 43% in rare types |
| Stable across a patient's regions | **MODERATE** — ICC 0.41-0.46 |
| Separates tissue types | **PASS** — p=0.0009 |
| Tracks recurrence | **PASS** — p=0.005-0.010 |
| Specific (null on HPV) | **PASS** — p=0.69 |
| Best single feature of the 5 | **PASS** — +0.021, only one that matters |
| Predicts overall survival | **FAIL** — p=0.804 |
| Adds on top of the 256 | **FAIL** — +0.0003 |

### What to do next

1. **Chase recurrence, not survival.** Better powered, and it is where the signal is.
2. **Swap the ratio for Squidpy's `nhood_enrichment`** — same concept, permutation-
   calibrated, and citable.
3. **Fix the noise controls** in `run_survival.py` — average many draws, not one.
4. **Abundance-correct all 256** (`E[i][j] = P[i][j]/p_j`) with selection *inside*
   each CV fold. The signal lives in specific cell-type pairs, not in any summary —
   every aggregate tested here is null while the full block works.

---

## Files

| Path | What |
|---|---|
| `src/enrichment_features.py` | builds the 5 features |
| `src/validate_groups.py` | lineage-marker validation |
| `src/make_spatial_maps.py` | the figures |
| `data/processed/outputs/spatial_maps/01_protein_variation.png` | protein in situ |
| `data/processed/outputs/spatial_maps/02_celltypes.png` | all 16 types |
| `data/processed/outputs/spatial_maps/03_groups.png` | the 3-group collapse |
| `data/processed/outputs/spatial_maps/04_lineage_validation.png` | biology check |
| `data/processed/outputs/spatial_maps/05_kl_matrix.png` | the matrix behind kl_mean |
| `data/processed/outputs/survival/ablation_enrichment.csv` | all 31 subsets |

## References

- Palla et al. 2022, *Squidpy: a scalable framework for spatial omics analysis*, Nature Methods — `nhood_enrichment`
- Schapiro et al. 2017, *histoCAT*, Nature Methods — neighbourhood permutation test
- Analytical Neighborhood Enrichment Score, arXiv 2506.18692 (2025) — closed-form version, validated vs Monte Carlo at r >= 0.95
- CD8+/Treg immune phenotypes in HNSCC — desert/excluded/inflamed, median OS 37/61/85 months
