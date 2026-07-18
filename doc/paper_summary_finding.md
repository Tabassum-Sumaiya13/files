# Spatial Survival Features — Combined Summary

*One summary of three reports: Enrichment features, `kl_mean` (deep dive), and Node/marker features. Spatial positional-encoding survival pipeline on a cancer CODEX atlas.*

---

## 1. The big picture in one diagram

The whole project asks one question in three ways:

```
   TISSUE  =  a graph of cells
              (each cell = a node, edges = who is next to whom)

              │
              ▼
   Turn the tissue into a short list of numbers per patient
              │
              ▼
   RandomSurvivalForest  →  predict who does worse

   BASELINE = "which cell types are present" (16 proportions)  → C ≈ 0.69
   QUESTION = does adding a new feature beat a RANDOM block of the same size?
```

That last line is the honest test used everywhere in these reports. It is called **`ΔC vs NOISE`**. Adding any columns can help by luck, so each new feature is compared to a same-size block of random numbers. Beating that is real signal.

> **One-line takeaway:** All three feature families add **real but small** signal *on top of* knowing the cell types. **None of them work on their own.** And in every case, **fewer features win** — a couple of good numbers beat a big pile.

---

## 2. What the three reports cover

| Report | Feature family | Size | Core idea |
|---|---|---|---|
| **Enrichment** | 5 spatial scalars | 5 | Squash a 16×16 neighbour table into 5 "how organised is the tissue" numbers, after removing what cell counts already explain |
| **kl_mean** | 1 scalar (deep dive) | 1 | The single strongest enrichment feature. "How different is a cell's local neighbourhood from the tissue average" |
| **Node features** | protein markers | 8 → 2 | Bring back the measured proteins (e.g. FoxP3, Ki67), each read only in the cell type where it means something |

*Note on data: the Enrichment and Node reports name a **colorectal** cohort (Schürch 2020); the kl_mean report names a **UPMC Head & Neck** cohort. The method is identical; the tissue label differs between the source documents. Event counts line up (~27 patient-level, ~103 acquisition-level).*

---

## 3. The one shared trick — "divide by chance"

Two of the three families (Enrichment + kl_mean) rest on the same move. It is worth seeing once.

```
P[i][j] = "of cell type i's neighbours, what fraction are type j"

Problem: if type j is 18% of the whole tissue, then even RANDOM mixing
         makes j show up ~18% of the time. That 18% is not organisation.
         And the baseline model ALREADY knows the 18%.

Fix:     E[i][j] = P[i][j] / p_j

         E = 1  → exactly chance (delete it, tells us nothing)
         E > 1  → genuinely clumped together
         E < 1  → genuinely avoiding each other
```

This deletes the part the baseline already has, and keeps only the spatial part.

---

## 4. Enrichment features (the 5 scalars)

| Feature | Plain meaning |
|---|---|
| `kl_mean` | how organised the whole tissue is (0 = random soup) |
| `kl_tumor` | how distinctive the tumour's neighbourhood is |
| `self_enrich` | do cells clump with their own kind |
| `immune_tumor` | do immune cells get into the tumour |
| `stroma_tumor` | do stromal cells sit at the tumour edge |

**Biology it revealed:** cells strongly sit with their own kind, and immune + stromal cells are **pushed out of the tumour** in almost every sample — the classic "cold / immune-excluded tumour."

**Main findings:**

| Feature set | Feats | C-index | vs NOISE | Verdict |
|---|---|---|---|---|
| Celltype + delaunay (full table) | 272 | 0.718 | +0.062 | real signal |
| **Celltype + Enrichment** | **21** | **0.716** | +0.028 | **real signal** |
| Celltype baseline | 16 | 0.690 | — | baseline |
| Enrichment alone | 5 | 0.580 | — | worse than baseline |

- **Efficiency win:** 5 numbers recover almost all of what the 256-column table gives (0.716 vs 0.718).
- **But standalone-weak:** on their own (0.580) they are worse than baseline.
- **Very redundant with each other** (correlations up to 0.80). Effectively 1–2 real dimensions, not 5. Inside the block, **only `kl_mean` reliably helps** — dropping any of the other four actually *improves* the model.

---

## 5. `kl_mean` — the one feature that does the work

This is the strongest of the 5, studied on its own.

```
kl_mean = 0     →  well-mixed soup, no organisation
kl_mean = high  →  cells sit with specific partners, tissue is organised
```

**Is it real? The shuffle test (the key check):**

```
observed kl_mean   0.552
SHUFFLED labels    0.035   ← noise floor
                   -----
real structure is ~15× the noise, in 307 / 307 samples
```

**PASS** — it detects genuine structure.

**What it means biologically (all point the same way):**

```
Normal mucosa   0.606   ← most organised
Nodal met       0.587
Primary tumor   0.514
Recurred        0.425
Recurrence tumor 0.321  ← least organised
```

**Losing spatial organisation tracks a more aggressive disease.** It separates tissue types (p = 0.0009) and tracks **recurrence** (p = 0.005–0.010).

**Where it stumbles:**

| Test | Result |
|---|---|
| Detects real structure | PASS |
| Separates tissue types | PASS |
| Tracks recurrence | PASS |
| Predicts overall survival | **FAIL** (p = 0.80) |
| Adds anything on top of the full 256 table | **FAIL** (+0.0003) |
| Free of rare-cell-type noise | PARTIAL (43% noise when a type has <10 cells) |

Two honest caveats: it does **not** predict overall survival here (survival has only 27 events and counts *any* death), and it adds **nothing** on top of the full neighbour table — because it is just a summary of that same table. **Pick one, not both.**

---

## 6. Node features (bring back the proteins)

Until this work, only the *cell type label* reached the model. The 39 measured proteins were dropped. The idea: add a few proteins back — but read each one **only in the cell type where it is meaningful**.

```
Ki67 in a tumour cell    →  the tumour is dividing fast   (meaningful)
Ki67 averaged over ALL   →  noise                         (meaningless)
```

Why not all 39? With so few patient deaths, hundreds of columns would drown the signal. So a short, biology-backed list was used (FoxP3, GranzymeB, PD-1, Ki67, etc.).

**Main findings:**

| Feature set | Feats | C-index | vs NOISE | Verdict |
|---|---|---|---|---|
| **Celltype + Enrichment + Markers(2)** | **23** | **0.733** | **+0.053** | **best result** |
| Celltype + Markers (full 8) | 24 | 0.716 | +0.042 | real signal |
| Celltype baseline | 16 | 0.690 | — | baseline |
| Markers alone (8) | 8 | 0.648 | — | worse than baseline |

- **Less is more, again:** cutting 8 markers down to **just 2** (`cd4_foxp3` + `tcell_cd45ro`) *raised* the score and tightened the result.
- Those two markers alone (Treg suppression + T-cell memory) nearly matched the whole 16-way cell-type mix — the immune functional state really matters for outcome.
- **Complementary, not a substitute:** markers alone (0.648) lose to the baseline. Knowing *which* cells are present beats knowing their *state* — but combining both wins.

---

## 7. The three lessons that repeat in every report

```
1. COMPLEMENTARY, NOT STANDALONE
   Every family loses to the baseline when used alone.
   They only help ON TOP OF the cell-type composition.

2. LESS IS MORE
   Enrichment: usable signal is basically 1 feature (kl_mean).
   Markers:    8 → 2 markers scored HIGHER.
   Big blocks dilute the signal across too many columns.

3. REDUNDANT WITH THE BIG TABLE
   The full 256-column neighbour table already contains the
   enrichment/kl_mean information. Adding them on top = ~zero.
   Choose: the 256 for the top raw number, OR the small blocks
   for a simple, readable, robust model.
```

---

## 8. Bottom line and recommendation

| Goal | Use this | Score |
|---|---|---|
| **Best single number** | Celltype + Enrichment + Markers(2) | **C = 0.733** |
| Best raw spatial number | Celltype + full 256 table | C = 0.718 |
| Simple + robust | Celltype + kl_mean (1 feature) or + Markers(2) | ~0.703–0.718 |

**Recommended setup:** the two-marker node block plus the enrichment block on top of cell-type composition.

```bash
python run_survival.py --experiment C --with-markers --marker-keep cd4_foxp3,tcell_cd45ro
```

**Suggested next steps (from the reports):**
- Chase **recurrence**, not overall survival — it is better powered and it is where the signal actually lives.
- Swap the home-made ratio for a permutation-calibrated tool (Squidpy `nhood_enrichment`).
- Fix the noise controls — average many random draws, not one, because wide blocks are unstable.
- Abundance-correct the full 256 table and pick features *inside* each cross-validation fold — the signal sits in specific cell-type pairs, not in any single summary.