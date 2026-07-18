# Data Flow: Input → Feature → Model → Output

---

## Common Input Files (Used by All Features)

### 1. marker_names.csv
```
┌────────────────┐
│    marker      │
├────────────────┤
│ CD31           │
│ CD57           │
│ CD4            │
│ CD45RA         │
│ CD8            │
│ ...            │
│ (40 markers)   │
└────────────────┘
```

### 2. sample_metadata.csv (Survival Labels)
```
┌───────────────────────────────┬────────────┬──────────────────┬─────────────────┬─────────────────┐
│ acquisition_id                │ patient_id │ survival_status  │ survival_day    │ tissuetype      │
├───────────────────────────────┼────────────┼──────────────────┼─────────────────┼─────────────────┤
│ UPMC_c001_v001_r001_reg001    │ 1          │ 1                │ 365             │ Primary tumor   │
│ UPMC_c001_v001_r001_reg002    │ 1          │ 1                │ 365             │ Primary tumor   │
│ UPMC_c002_v001_r001_reg001    │ 2          │ 0                │ 1825            │ Normal mucosa   │
│ UPMC_c003_v001_r001_reg001    │ 3          │ 1                │ 548             │ Primary tumor   │
│ ...                           │ ...        │ ...              │ ...             │ ...             │
└───────────────────────────────┴────────────┴──────────────────┴─────────────────┴─────────────────┘
```

### 3. labeled_arcsinh_norm_data.pkl (Cell Expression Data)
```
┌──────────────────────────────┬─────────┬─────────┬─────────┬─────────┬─────────┬─────────────────────┐
│ sample_id                    │ cell_id │ CD31    │ CD57    │ CD4     │ ...     │ cluster_label       │
├──────────────────────────────┼─────────┼─────────┼─────────┼─────────┼─────────┼─────────────────────┤
│ UPMC_c001_v001_r001_reg001   │ 1       │ 0.523   │ 0.102   │ 2.145   │ ...     │ CD4 T cell          │
│ UPMC_c001_v001_r001_reg001   │ 2       │ 0.087   │ 0.034   │ 0.056   │ ...     │ Tumor (Ki67+)       │
│ UPMC_c001_v001_r001_reg001   │ 3       │ 1.234   │ 0.876   │ 0.123   │ ...     │ Macrophage          │
│ ...                          │ ...     │ ...     │ ...     │ ...     │ ...     │ ...                 │
│ (2 million rows)             │         │         │         │         │         │                     │
└──────────────────────────────┴─────────┴─────────┴─────────┴─────────┴─────────┴─────────────────────┘
```

---

## Feature 1: Biomarker Region (40 features)- protein in all tissue ( + empty space)

### **Input File: `expression_biomarkers/{sample}_cell_info.csv`**
```
┌─────────────┬─────────┬─────────┬─────────┬─────────┬─────────┐
│             │ CD31    │ CD57    │ CD4     │ CD8     │ ...     │
├─────────────┼─────────┼─────────┼─────────┼─────────┼─────────┤
│ sum         │ 52341.2 │ 8234.5  │ 12456.7 │ 9876.3  │ ...     │
│ average     │ 25.42   │ 4.00    │ 6.05    │ 4.80    │ ...     │  ← Uses this row
│ max         │ 255.0   │ 198.3   │ 245.6   │ 201.2   │ ...     │
└─────────────┴─────────┴─────────┴─────────┴─────────┴─────────┘
```
```
Input File (raw):           arcsinh()           Feature Matrix:
average = 25.42      →    arcsinh(25.42)    →     3.24
average = 4.00       →    arcsinh(4.00)     →     2.08
average = 6.05       →    arcsinh(6.05)     →     2.51

```
```
Sample: UPMC_c001_v001_r001_reg001

Step 1: Read file expression_biomarkers/UPMC_c001_v001_r001_reg001_cell_info.csv

Step 2: Get 'average' row values:
        CD31=25.42, CD57=4.00, CD4=6.05, CD8=4.80, ...

Step 3: Apply arcsinh() to each value:
        arcsinh(25.42) = 3.93  → stored as CD31 feature
        arcsinh(4.00)  = 2.09  → stored as CD57 feature
        arcsinh(6.05)  = 2.50  → stored as CD4 feature
        arcsinh(4.80)  = 2.27  → stored as CD8 feature

Step 4: Merge with survival labels from sample_metadata.csvunsupervised methods to
cluster cells rather than supervised cell type assignment when
constructing the neighborhood matrix and non-cell-based K-
functions for biomarker
```
### Feature Matrix (Input to Model)
```
┌───────────────────────────────┬─────────┬─────────┬─────────┬─────────┬──────────────────┬─────────────────┐
│ acquisition_id                │ CD31    │ CD57    │ CD4     │ ...     │ survival_status  │ survival_day    │
├───────────────────────────────┼─────────┼─────────┼─────────┼─────────┼──────────────────┼─────────────────┤
│ UPMC_c001_v001_r001_reg001    │ 3.24    │ 2.08    │ 2.51    │ ...     │ 1                │ 365             │
│ UPMC_c002_v001_r001_reg001    │ 2.87    │ 1.92    │ 2.78    │ ...     │ 0                │ 1825            │
│ UPMC_c003_v001_r001_reg001    │ 3.56    │ 2.34    │ 2.12    │ ...     │ 1                │ 548             │
│ ...                           │ ...     │ ...     │ ...     │ ...     │ ...              │ ...             │
└───────────────────────────────┴─────────┴─────────┴─────────┴─────────┴──────────────────┴─────────────────┘
Shape: (~307 samples × 40 features)
```

---

## Feature 2: Biomarker Cell (40 features)- protein only inside cell

### Computation: Group by sample → Mean of each protein
```python
#** From labeled_arcsinh_norm_data.pkl**
cell_biomarkers.groupby('acquisition_id').mean(numeric_only=True)
```
```
┌──────────────────────────────┬─────────┬─────────┬─────────┬─────────┬─────────┐
│ sample_id                    │ cell_id │ CD31    │ CD57    │ CD4     │ ...     │
├──────────────────────────────┼─────────┼─────────┼─────────┼─────────┼─────────┤
│ UPMC_c001_v001_r001_reg001   │ 1       │ 0.523   │ 0.102   │ 2.145   │ ...     │
│ UPMC_c001_v001_r001_reg001   │ 2       │ 0.087   │ 0.034   │ 0.056   │ ...     │
│ UPMC_c001_v001_r001_reg001   │ 3       │ 1.234   │ 0.876   │ 0.123   │ ...     │
│ UPMC_c001_v001_r001_reg001   │ 4       │ 1.567   │ 0.189   │ 2.456   │ ...     │
│ UPMC_c001_v001_r001_reg001   │ 5       │ 1.045   │ 0.034   │ 1.987   │ ...     │
│ ... (5000 cells in sample 1) │ ...     │ ...     │ ...     │ ...     │ ...     │
├──────────────────────────────┼─────────┼─────────┼─────────┼─────────┼─────────┤
│ UPMC_c002_v001_r001_reg001   │ 1       │ 0.654   │ 0.234   │ 1.876   │ ...     │
│ UPMC_c002_v001_r001_reg001   │ 2       │ 0.789   │ 0.456   │ 0.987   │ ...     │
│ ... (6000 cells in sample 2) │ ...     │ ...     │ ...     │ ...     │ ...     │
└──────────────────────────────┴─────────┴─────────┴─────────┴─────────┴─────────┘
```
```
Sample: UPMC_c001_v001_r001_reg001 (has 5000 cells)

Step 1: Filter all cells belonging to this sample
        → 5000 rows

Step 2: Calculate MEAN of each protein column across all 5000 cells:

        CD31 values: [0.523, 0.087, 1.234, 1.567, 1.045, ...]  (5000 values)
        CD31 mean = (0.523 + 0.087 + 1.234 + 1.567 + 1.045 + ...) / 5000
        CD31 mean = 0.892  ← This goes in the feature matrix!

        CD57 values: [0.102, 0.034, 0.876, 0.189, 0.034, ...]  (5000 values)
        CD57 mean = (0.102 + 0.034 + 0.876 + 0.189 + 0.034 + ...) / 5000
        CD57 mean = 0.234  ← This goes in the feature matrix!

        (Repeat for all 40 proteins)

Step 3: Merge with survival labels from sample_metadata.csv
```

### Feature Matrix (Input to Model)
```
┌───────────────────────────────┬─────────┬─────────┬─────────┬─────────┬──────────────────┬─────────────────┐
│ acquisition_id                │ CD31    │ CD57    │ CD4     │ ...     │ survival_status  │ survival_day    │
├───────────────────────────────┼─────────┼─────────┼─────────┼─────────┼──────────────────┼─────────────────┤
│ UPMC_c001_v001_r001_reg001    │ 0.892   │ 0.234   │ 1.456   │ ...     │ 1                │ 365             │
│ UPMC_c002_v001_r001_reg001    │ 0.756   │ 0.312   │ 1.234   │ ...     │ 0                │ 1825            │
│ UPMC_c003_v001_r001_reg001    │ 0.945   │ 0.189   │ 1.567   │ ...     │ 1                │ 548             │
│ ...                           │ ...     │ ...     │ ...     │ ...     │ ...              │ ...             │
└───────────────────────────────┴─────────┴─────────┴─────────┴─────────┴──────────────────┴─────────────────┘
Shape: (~307 samples × 40 features)
```

---

## Feature 3: Cell Type Proportion (16 features)

### Computation: Count cells per type → Normalize to proportions
```python
# From labeled_arcsinh_norm_data.pkl
expr.groupby(['acquisition_id', 'cluster_label']).size() / total_cells
```
```
┌──────────────────────────────┬─────────┬─────────────────────┐
│ sample_id                    │ cell_id │ cluster_label       │
├──────────────────────────────┼─────────┼─────────────────────┤
│ UPMC_c001_v001_r001_reg001   │ 1       │ CD4 T cell          │
│ UPMC_c001_v001_r001_reg001   │ 2       │ Tumor (Ki67+)       │
│ UPMC_c001_v001_r001_reg001   │ 3       │ Tumor (Ki67+)       │
│ UPMC_c001_v001_r001_reg001   │ 4       │ CD4 T cell          │
│ UPMC_c001_v001_r001_reg001   │ 5       │ Macrophage          │
│ UPMC_c001_v001_r001_reg001   │ 6       │ Tumor (Ki67+)       │
│ UPMC_c001_v001_r001_reg001   │ 7       │ CD8 T cell          │
│ UPMC_c001_v001_r001_reg001   │ 8       │ Tumor (Ki67+)       │
│ UPMC_c001_v001_r001_reg001   │ 9       │ B cell              │
│ UPMC_c001_v001_r001_reg001   │ 10      │ Macrophage          │
│ ... (5000 cells total)       │ ...     │ ...                 │
└──────────────────────────────┴─────────┴─────────────────────┘
```
```
Sample: UPMC_c001_v001_r001_reg001 (has 5000 cells total)

Step 1: COUNT cells of each type in this sample

        Cell Type         │ Count
        ──────────────────┼───────
        CD4 T cell        │  600
        CD8 T cell        │  400
        Macrophage        │  750
        Tumor (Ki67+)     │ 1750
        B cell            │  300
        Stromal           │  500
        ... (16 types)    │  ...
        ──────────────────┼───────
        TOTAL             │ 5000

Step 2: DIVIDE by total cells to get PROPORTION

        Cell Type         │ Count │ Proportion
        ──────────────────┼───────┼────────────
        CD4 T cell        │  600  │ 600/5000  = 0.12
        CD8 T cell        │  400  │ 400/5000  = 0.08
        Macrophage        │  750  │ 750/5000  = 0.15
        Tumor (Ki67+)     │ 1750  │ 1750/5000 = 0.35
        B cell            │  300  │ 300/5000  = 0.06
        ... (16 types)    │  ...  │ ...
        ──────────────────┴───────┴────────────
        SUM of proportions = 1.00 (always!)

Step 3: Store in feature matrix row
```

### Feature Matrix (Input to Model)
```
┌───────────────────────────────┬───────────┬───────────┬─────────────┬───────────┬──────────────────┬─────────────────┐
│ acquisition_id                │ CD4_Tcell │ CD8_Tcell │ Macrophage  │ Tumor     │ survival_status  │ survival_day    │
├───────────────────────────────┼───────────┼───────────┼─────────────┼───────────┼──────────────────┼─────────────────┤
│ UPMC_c001_v001_r001_reg001    │ 0.12      │ 0.08      │ 0.15        │ 0.35      │ 1                │ 365             │
│ UPMC_c002_v001_r001_reg001    │ 0.18      │ 0.14      │ 0.10        │ 0.22      │ 0                │ 1825            │
│ UPMC_c003_v001_r001_reg001    │ 0.05      │ 0.03      │ 0.08        │ 0.52      │ 1                │ 548             │
│ ...                           │ ...       │ ...       │ ...         │ ...       │ ...              │ ...             │
└───────────────────────────────┴───────────┴───────────┴─────────────┴───────────┴──────────────────┴─────────────────┘
Shape: (~307 samples × 16 features)
```

---

## Feature 4: Neighborhood Matrix (256 features)

### Input File: `celltype_neighborhoods/k10/{sample}.npy`
```
Pre-computed 16×16 numpy array per sample:
Each sample has a pre-computed 16×16 matrix counting neighbor relationships:

File: celltype_neighborhoods/k10/UPMC_c001_v001_r001_reg001.npy

For each cell, we found its 10 nearest neighbors (k=10), then counted:
"How many neighbors of each type does each cell type have?"

RAW COUNTS MATRIX:
                    │ Tumor │ CD4T │ CD8T │ Macro │ Bcell │ ...  │ ROW SUM
────────────────────┼───────┼──────┼──────┼───────┼───────┼──────┼─────────
Tumor (1750 cells)  │ 7350  │ 2100 │ 1400 │ 2625  │  875  │ ...  │ 17500
CD4 T cell (600)    │  600  │ 1500 │ 1080 │  480  │  720  │ ...  │  6000
CD8 T cell (400)    │  320  │  600 │ 1200 │  400  │  240  │ ...  │  4000
Macrophage (750)    │ 1500  │  750 │  900 │ 2625  │  600  │ ...  │  7500
...                 │  ...  │  ... │  ... │  ...  │  ...  │ ...  │   ...
```
row normalization
```
NORMALIZED MATRIX (proportions):
                    │ Tumor │ CD4T │ CD8T │ Macro │ Bcell │ ...  │ ROW SUM
────────────────────┼───────┼──────┼──────┼───────┼───────┼──────┼─────────
Tumor               │ 0.42  │ 0.12 │ 0.08 │ 0.15  │ 0.05  │ ...  │  1.00
CD4 T cell          │ 0.10  │ 0.25 │ 0.18 │ 0.08  │ 0.12  │ ...  │  1.00
CD8 T cell          │ 0.08  │ 0.15 │ 0.30 │ 0.10  │ 0.06  │ ...  │  1.00
Macrophage          │ 0.20  │ 0.10 │ 0.12 │ 0.35  │ 0.08  │ ...  │  1.00
...                 │  ...  │  ... │  ... │  ...  │  ...  │ ...  │   ...

Calculation example:
Tumor→Tumor = 7350 / 17500 = 0.42
Tumor→CD4T  = 2100 / 17500 = 0.12
```

```
Sample: UPMC_c001_v001_r001_reg001

INPUT (.npy file):
┌─────────────────────────────────────────┐
│  Raw 16×16 count matrix                 │
│  [7350, 2100, 1400, ...]                │
│  [600,  1500, 1080, ...]                │
│  ...                                    │
└─────────────────────────────────────────┘
            │
            ▼ Row normalize
            │
┌─────────────────────────────────────────┐
│  Normalized 16×16 matrix                │
│  [0.42, 0.12, 0.08, ...]                │
│  [0.10, 0.25, 0.18, ...]                │
│  ...                                    │
└─────────────────────────────────────────┘
            │
            ▼ Flatten
            │
┌─────────────────────────────────────────────────────────────────────┐
│  Feature vector (256 values)                                        │
│  [0.42, 0.12, 0.08, 0.15, ..., 0.10, 0.25, 0.18, ..., 0.08, ...]   │
│   T→T   T→CD4 T→CD8 T→Mac      CD4→T CD4→CD4 CD4→CD8                │
└─────────────────────────────────────────────────────────────────────┘
            │
            ▼ Add to DataFrame row
            │
OUTPUT (Feature Matrix):
┌─────────────────────────────┬───────┬────────┬────────┬─────┬─────────────────┐
│ acquisition_id              │ T→T   │ T→CD4  │ T→CD8  │ ... │ survival_status │
├─────────────────────────────┼───────┼────────┼────────┼─────┼─────────────────┤
│ UPMC_c001_v001_r001_reg001  │ 0.42  │ 0.12   │ 0.08   │ ... │ 1               │
└─────────────────────────────┴───────┴────────┴────────┴─────┴─────────────────┘
```

### Feature Matrix (Input to Model)
```
┌───────────────────────────────┬────────────────┬────────────────┬────────────────┬─────┬──────────────────┬─────────────────┐
│ acquisition_id                │ Tumor→Tumor    │ Tumor→CD4T     │ Tumor→CD8T     │ ... │ survival_status  │ survival_day    │
├───────────────────────────────┼────────────────┼────────────────┼────────────────┼─────┼──────────────────┼─────────────────┤
│ UPMC_c001_v001_r001_reg001    │ 0.42           │ 0.12           │ 0.08           │ ... │ 1                │ 365             │
│ UPMC_c002_v001_r001_reg001    │ 0.28           │ 0.18           │ 0.15           │ ... │ 0                │ 1825            │
│ UPMC_c003_v001_r001_reg001    │ 0.55           │ 0.05           │ 0.03           │ ... │ 1                │ 548             │
│ ...                           │ ...            │ ...            │ ...            │ ... │ ...              │ ...             │
└───────────────────────────────┴────────────────┴────────────────┴────────────────┴─────┴──────────────────┴─────────────────┘
Shape: (~307 samples × 256 features)
```

---

## Feature 5: Ripley K Function (80 features)

### Input File: `k_fns_norm_by_uw_qc_labeled.csv`

pre computed 
```
RAW DATA                          R SCRIPT                      PYTHON
──────────────────────────────────────────────────────────────────────────

cell_locations_and_labels.csv     compute_k_fns.R               rsf_prediction.py
(X, Y coordinates)                      │                            │
        │                               │                            │
        ▼                               ▼                            ▼
┌─────────────────┐             ┌──────────────────┐         ┌────────────────┐
│ Cell positions  │   ──────►   │ Calculate K(r)   │  ──────►│ Load CSV       │
│ + protein expr  │             │ for each marker  │         │ + merge labels │
└─────────────────┘             │ at r=30,80       │         │ + train RSF    │
                                └──────────────────┘         └────────────────┘
                                        │
                                        ▼
                                k_fns_norm_by_uw_qc_labeled.csv
                                (307 samples × 80 features)

```                               
```
┌───────────────────────────────┬──────────────┬──────────────┬──────────────┬──────────────┬─────────┐
│ acquisition_id                │ CD31_r30     │ CD31_r80     │ CD57_r30     │ CD57_r80     │ ...     │
├───────────────────────────────┼──────────────┼──────────────┼──────────────┼──────────────┼─────────┤
│ UPMC_c001_v001_r001_reg001    │ 1.234        │ 2.567        │ 0.876        │ 1.234        │ ...     │
│ UPMC_c002_v001_r001_reg001    │ 0.987        │ 1.876        │ 1.123        │ 1.567        │ ...     │
│ UPMC_c003_v001_r001_reg001    │ 1.456        │ 2.890        │ 0.654        │ 0.987        │ ...     │
│ ...                           │ ...          │ ...          │ ...          │ ...          │ ...     │
└───────────────────────────────┴──────────────┴──────────────┴──────────────┴──────────────┴─────────┘
(40 markers × 2 radii = 80 features)
```
```
CD31_r30 = 1.234 (clustered):
"Cells with high CD31 (endothelial marker) are clustered together"


CD57_r30 = 0.876 (dispersed):
"Cells with high CD57 (NK cell marker) are spread out"



Result: 
- Value > 1 → CLUSTERED (more neighbors than random)
- Value = 1 → RANDOM
- Value < 1 → DISPERSED (fewer neighbors than random)
```

### Feature Matrix (Input to Model)
```
┌───────────────────────────────┬──────────────┬──────────────┬──────────────┬─────┬──────────────────┬─────────────────┐
│ acquisition_id                │ CD31_r30     │ CD31_r80     │ CD57_r30     │ ... │ survival_status  │ survival_day    │
├───────────────────────────────┼──────────────┼──────────────┼──────────────┼─────┼──────────────────┼─────────────────┤
│ UPMC_c001_v001_r001_reg001    │ 1.234        │ 2.567        │ 0.876        │ ... │ 1                │ 365             │
│ UPMC_c002_v001_r001_reg001    │ 0.987        │ 1.876        │ 1.123        │ ... │ 0                │ 1825            │
│ UPMC_c003_v001_r001_reg001    │ 1.456        │ 2.890        │ 0.654        │ ... │ 1                │ 548             │
│ ...                           │ ...          │ ...          │ ...          │ ... │ ...              │ ...             │
└───────────────────────────────┴──────────────┴──────────────┴──────────────┴─────┴──────────────────┴─────────────────┘
Shape: (~307 samples × 80 features)
```

---
IMPORTANT: The SAME 10 folds are used for all features!
(Same patients in same folds - fair comparison)
```
┌─────────────────────────────────────────────────────────────────────────────┐
│                                                                             │
│  RUN 1: biomarker_region (40 features)                                      │
│  ┌─────────────────────────────────────────────────────────────────────┐   │
│  │ Fold 0 → C=0.655 │ Fold 1 → C=0.648 │ ... │ Fold 9 → C=0.661       │   │
│  └─────────────────────────────────────────────────────────────────────┘   │
│  Mean C-Index = 0.655                                                       │
│                                                                             │
│  RUN 2: biomarker_cell (40 features)                                        │
│  ------similar for other feature                                                        │
└─────────────────────────────────────────────────────────────────────────────┘
```
```
┌─────────────────────────────────────────────────────────────────────────┐
│                                                                         │
│  FOLD 0:                                                                │
│                                                                         │
│  ┌─────────────────┐      ┌─────────────────┐      ┌─────────────────┐ │
│  │ Train Data      │      │ RSF Model       │      │ Test Data       │ │
│  │ (276 samples)   │──────► (100 trees)     │──────► (31 samples)    │ │
│  │                 │ fit()│                 │score()│                 │ │
│  └─────────────────┘      └─────────────────┘      └─────────────────┘ │
│                                  │                         │            │
│                                  │                         ▼            │
│                                  │              ┌─────────────────────┐ │
│                                  │              │ rsf.score() returns │ │
│                                  │              │ C-Index = 0.712     │ │
│                                  │              └─────────────────────┘ │
│                                  │                         │            │
│                                  ▼                         ▼            │
│                           cv_concord[0] = 0.712                        │
│                                                                         │
└─────────────────────────────────────────────────────────────────────────┘
```
## Output Files (After Training)
![alt text](image-5.png)
![alt text](image-6.png)

### Directory: `rsf_risk_scores/`
```
rsf_risk_scores/
├── avg_biomarker_region_concordance.npy    ← biomarker_region results
├── avg_biomarker_cell_concordance.npy      ← biomarker_cell results
├── celltype_proportions_concordance.npy    ← celltype_prop results
├── neighbor_mat_concordance.npy            ← neighbor_mat results
└── ripley_concordance.npy                  ← ripley results
```

### 1. `{feature}_concordance.npy`--- c-index
```python
# Shape: (10,) - One C-index per CV fold
array([0.712, 0.698, 0.725, 0.701, 0.718, 0.695, 0.732, 0.708, 0.699, 0.715])

# Mean C-index = 0.710

{feature}_concordance.npy = [0.712, 0.698, 0.725, ...]

┌─────────────────────────────────────────────────────────┐
│  10 values = C-Index for each cross-validation fold     │
│                                                         │
│  Each value = "What % of patient pairs were correctly   │
│                ranked by risk in this test fold?"       │
│                                                         │
│  Mean = 0.710 → Overall model performance               │
│  Std  = 0.012 → Model is consistent across folds        │
└─────────────────────────────────────────────────────────┘
```

### 2. `{feature}_pred_risks.npy`
```python
# Shape: (307,) - Risk score for each sample
array([2.456, 1.234, 3.567, 0.987, 2.123, ...])

# Higher value = Higher risk of death
```

### 3. `{feature}_feature_imp.npy`
Permutation Importance Method
```
ORIGINAL MODEL: C-Index = 0.712 (on test set)

After shuffling Feature 0 (Tumor→Tumor):
C-Index = 0.698  (drops by 0.014)

→ Feature 0 importance = 0.014 (how much performance decreased)
```
```python
# Shape: (10, num_features, 2) - [mean, std] importance per fold
# Example for neighbor_mat: (10, 256, 2)

# Fold 0, Feature 0 (Tumor→Tumor): mean=0.012, std=0.003
# Fold 0, Feature 1 (Tumor→CD4T): mean=0.008, std=0.002
# ...
```
```┌─────────────────────────────────────────────────────────────────────┐
│                                                                     │
│  {feature}_feature_imp.npy                                          │
│                                                                     │
│  • Shape: (10, num_features, 2)                                    │
│  • [fold, feature, 0] = mean importance (performance drop)        │
│  • [fold, feature, 1] = std importance (variability)               │
│                                                                     │
│  Calculated by:                                                     │
│  1. Train RSF on fold's training data                              │
│  2. For each feature: shuffle 15 times, measure C-index drop       │
│  3. Store mean and std of drops per feature per fold               │
│                                                                     │
│  Higher values = More important features                           │
│                                                                     │
└─────────────────────────────────────────────────────────────────────┘
```

### 4. `{feature}_split_info.npy`
tracks which samples were used for testing in each cross-validation fold. This allows you to:
```python
# Shape: (10, 307) - Which samples in test set per fold
# 1 = test sample, 0 = train sample

array([[0, 0, 1, 1, 0, 0, 0, 1, ...],  # Fold 0
       [1, 1, 0, 0, 1, 0, 0, 0, ...],  # Fold 1
       ...])
```

### 5. `estimators/{feature}_cv{i}_rsf.pkl`
```python
# Saved RandomSurvivalForest model for each fold
# Can be loaded with pickle.load() for later prediction
```



---