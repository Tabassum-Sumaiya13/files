# Spatial Proteomics Cancer Survival Prediction Pipeline


**Paper:** "Deriving spatial features from in situ proteomics imaging to enhance cancer survival analysis"  
**Authors:** Dayao et al.  
**Dataset:** UPMC Head & Neck Cancer Cohort


## Core question?
>Can we predict how long a cancer patient will survive based on the spatial arrangement of cells in their tumor?



## What is Head & Neck Squamous Cell Carcinoma (HNSCC)?

- **Location:** Mouth, throat, larynx, pharynx
- **Incidence:** ~65,000 new cases/year in the US
- **5-year survival:** 50-65% (varies by stage and HPV status)
- **HPV connection:** HPV-positive tumors have better prognosis


##  Data Source

| Property | Value |
|----------|-------|
| Institution | University of Pittsburgh Medical Center (UPMC) |
| Cancer Type | Head & Neck Squamous Cell Carcinoma (HNSCC) |
| Patients | ~80 unique patients |
| Samples | ~400 tissue regions (multiple per patient) |
| Cells | ~2,061,102 individual cells |
| Proteins | 40 biomarkers measured |
| Cell Types | 16 distinct cell types |
| Technology | Multiplexed Ion Beam Imaging (MIBI) |

<div style="page-break-after: always;"></div>

## Pipeline
```
┌──────────────────────────────────────────────────────────────────┐
│  TISSUE IMAGES → Individual Cells with Markers & Locations      │
└──────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌──────────────────────────────────────────────────────────────────┐
│  FEATURE EXTRACTION (5 different approaches)                     │
│                                                                  │
│  1. Biomarker Region: Average protein intensity per tissue area  │
│  2. Biomarker Cell: Average protein intensity within cells       │
│  3. Celltype Proportion: % of each cell type in the tissue       │
│  4. Neighborhood Matrix: Who neighbors whom? (cell interactions) │
│  5. Ripley's K: Are certain cells clustered or dispersed?        │
└──────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌──────────────────────────────────────────────────────────────────┐
│  RANDOM SURVIVAL FOREST (RSF)                                    │
│  - Machine learning model for survival prediction                │
│  - 10-fold cross-validation (split by patient)                   │
│  - Predicts: Who will survive longer?                            │
└──────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌──────────────────────────────────────────────────────────────────┐
│  OUTPUT: Concordance Index (0.5 = random, 1.0 = perfect)         │
│  + Feature Importance (which markers/cells matter most)          │
└──────────────────────────────────────────────────────────────────┘
```
<div style="page-break-after: always;"></div>

### **step 1: Raw data collection**
```
Multiplexed Tissue Images
         │
         ▼
┌─────────────────────────────────────────────────┐
│  1 Sample = 1 tissue region from 1 patient      │
│  • ~2 million cells per sample                  │
│  • Each cell has:                               │
│      - X-Y coordinates (location)               │
│      - Protein intensity (40 markers)           │
│  • Patient records: survival time, status       │
└─────────────────────────────────────────────────┘
         │
         ▼
   Files: qc_acq_ids_labeled.csv (sample IDs)
          sample_metadata.csv (patient info)
          cell_locations_and_labels.csv (cell data)

```
### **step 2: Normalize Data**
```
Raw protein intensities
         │
         ▼
   arcsinh() normalization
         │
         ▼
   Reduces extreme values, makes data comparable

```
<div style="page-break-after: always;"></div>

### **step 3: Cell type labeling**

```
Normalized data
         │
         ▼
   Clustering algorithm
         │
         ▼
   Each cell gets a type label (16 cell types)
         │
         ▼
   File: labeled_arcsinh_norm_data.pkl

```

### **step 4: Compute spatial feature & extract traditional proteomics feature**

```
Cell locations + Cell types
         │
 FEATURE EXTRACTION (5 different approaches)                      
                                                                 
 1. Biomarker Region: Average protein intensity per tissue area  
 2. Biomarker Cell: Average protein intensity within cells       
 3. Celltype Proportion: % of each cell type in the tissue       
 4. Neighborhood Matrix: Who neighbors whom? (cell interactions) 
 5. Ripley's K: Are certain cells clustered or dispersed?  

```
<div style="page-break-after: always;"></div>

### **step 5: train & evaluate**

```
Each feature set (separately)
         │
         ▼
┌─────────────────────────────────────┐
│  10-Fold Cross Validation           │
│  (Grouped by patient - no leakage)  │
└─────────────────────────────────────┘
         │
         ▼
┌─────────────────────────────────────┐
│  Random Survival Forest (RSF)       │
│  • 100 trees                        │
│  • Predicts: risk score (survival)  │
└─────────────────────────────────────┘
         │
         ▼
┌─────────────────────────────────────┐
│  Calculate:                         │
│  • C-Index (prediction accuracy)    │
│  • Feature importance (permutation) │
└─────────────────────────────────────┘
         │
         ▼
   Output: rsf_risk_scores/

```
<div style="page-break-after: always;"></div>

## file summary 
**inputs--**

```
dataset_info.tar.gz
├── labeled_arcsinh_norm_data.pkl  ──────► biomarker_cell (protein per cell type)
├── cell_locations_and_labels.csv  ──────► celltype_prop (cell type counts)
└── sample_metadata.csv            ──────► SURVIVAL DATA (target variable)

biomarker_expr_summary.tar.gz
└── expression_biomarkers/         ──────► biomarker_region (tissue-level protein)

all_k_neighborhood_mats.tar.gz
└── k1-k20/{acq_id}.npy            ──────► neighbor_mat (cell interactions)

k_fns_norm_by_uw_qc_labeled.csv    ──────► ripley (clustering patterns)

patwa_comparisons.tar.gz           ──────► biomarker_threshold, biomarker_interaction
denvar.tar.gz                      ──────► DenVar

```
**outputs--**
```
rsf_risk_scores/               → Risk predictions + feature importance

```

<div style="page-break-after: always;"></div>


## Raw Data Files
### **File 1: `cell_locations_and_labels.csv` (2,061,102 rows)**
Every single cell's position and type (2 million+ rows):

                                
• X, Y coordinates (WHERE the cell is)                            
• Cell type label (WHAT type of cell)                             
                                                                   
 Used for: neighbor_mat, celltype_prop, ripley  


![alt text](image-8.png)

### **File 2: `labeled_arcsinh_norm_data.pkl` (2,061,102 rows)**

Protein levels for every cell (normalized)
 • 2 million+ rows (one per cell)                                  
 • 40 protein levels per cell (HOW MUCH protein)                   
 • Already normalized with arcsinh()                               
                                                                    
 Used for: biomarker_cell, biomarker_region  


![alt text](image-3.png)


### **File 3: `sample_metadata.csv` (~400 rows)**

Patient survival info + clinical data
• ~400 rows (one per tissue sample)                               
• Survival status (0=alive, 1=dead)                               
• Survival days (how long)                                        
• Clinical info (HPV, cancer site, etc.)                          
                                                                   
Used for: SURVIVAL LABELS (target for prediction)   

![alt text](image-4.png)

<div style="page-break-after: always;"></div>

## Data collection method

**Multiplexed Ion Beam Imaging (MIBI)**

MIBI is a technology that can measure **40+ proteins simultaneously** at **single-cell resolution** while preserving **spatial information**.

```
MULTIPLEXED IMAGING (MIBI/IMC):
- 40+ markers simultaneously
- Metal-tagged antibodies
- Mass spectrometry detection
- Single-cell resolution

```
**Initial Data generation (MIBI)**

```
Step 1: TISSUE BIOPSY
        Patient tumor → Surgical removal → Tissue sample
        
Step 2: TISSUE SECTIONING  
        Tissue → Thin slices (4-5 μm) → Mounted on slides
        
Step 3: ANTIBODY LABELING
        40 different antibodies, each tagged with a unique metal isotope
        Each antibody binds to one specific protein
        
Step 4: MASS SPECTROMETRY IMAGING
        Ion beam scans tissue → Metals vaporized → Mass detected
        Result: Pixel-by-pixel protein intensity map
        
Step 5: CELL SEGMENTATION
        Image analysis identifies individual cells
        Each cell gets X-Y coordinates + protein profile
        
Step 6: CELL TYPE CLASSIFICATION
        Clustering based on protein markers → Cell type labels

```
<div style="page-break-after: always;"></div>

### **The 40 Protein Biomarkers**

| Category | Markers | Biological Role |
|----------|---------|-----------------|
| **T cell** | CD3e, CD4, CD8, CD45RA, CD45RO | T cell identification & activation |
| **B cell** | CD20, CD21, CD38 | B cell identification |
| **Myeloid** | CD11b, CD11c, CD14, CD68 | Macrophages, dendritic cells |
| **NK cell** | CD56, CD57, CD16 | Natural killer cells |
| **Activation** | CD69, CD134, ICOS, HLA-DR | Immune activation status |
| **Checkpoint** | PD1, PDL1, CD152 (CTLA4) | Immune checkpoints |
| **Regulatory** | FoxP3 | Regulatory T cells |
| **Cytotoxic** | GranzymeB | Cytotoxic function |
| **Proliferation** | Ki67 | Cell division |
| **Tumor** | PanCK, p16, Podoplanin | Tumor cell markers |
| **Stromal** | Vimentin, aSMA, CollagenIV | Stromal/fibroblast markers |
| **Other** | CD31, CD34, CD47, CD49f, CD117, CD15, TMEM16A | Various functions |

<div style="page-break-after: always;"></div>

### **The 16 Cell Type**s

Cell types were determined by clustering cells based on protein expression:

```
IMMUNE CELLS:
├── CD4 T cell          (helper T cells - coordinate immune response)
├── CD8 T cell          (cytotoxic T cells - kill tumor cells)
├── B cell              (antibody production)
├── Naive immune cell   (not yet activated)
├── Macrophage          (phagocytosis, antigen presentation)
├── APC                 (Antigen Presenting Cell)
├── Granulocyte         (innate immunity)
└── Regulatory/exhausted (suppressive immune cells)

TUMOR CELLS:
├── Tumor (Ki67+)       (actively proliferating tumor)
├── Tumor (CD21+)       (tumor expressing CD21)
├── Tumor (CD20+)       (tumor expressing CD20)
├── Tumor (Podo+)       (tumor expressing Podoplanin)
└── Squamous epithelium (potentially pre-cancerous)

STROMAL CELLS:
├── Stromal / Fibroblast (structural support)
├── Endothelial         (blood vessel lining)
└── Smooth muscle       (vessel walls)

```

<div style="page-break-after: always;"></div>


##  Preprocessing Steps

### **1. arcsinh Normalization**

**Problem:** Raw protein intensity values have extreme outliers and non-normal distribution.

**Solution:** Apply arcsinh transformation (inverse hyperbolic sine).
```
Before normalization:
Values: [0, 0.1, 0.5, 1, 10, 100, 1000]

After arcsinh:
Values: [0, 0.10, 0.48, 0.88, 3.0, 5.3, 7.6]

Effect: Compresses large values, preserves small differences

```

**Why needed:**

Protein intensity values are highly skewed (most cells low, few very high)
Log transformation is common but can't handle zeros
Arcsinh handles zeros and reduces skewness


### **2. Cell Type Clustering**

Cells are grouped into 16 types using k-nearest neighbors (kNN) classification based on protein expression profiles:
```
16 clusters identified:
├── Based on known marker combinations
├── CD3+CD4+ → CD4 T cell
├── CD3+CD8+ → CD8 T cell
├── CD20+ → B cell
├── PanCK+ → Tumor/Epithelial
└── etc.

┌─────────────────────────────────────────────────────────────┐
│  CLUSTERING ALGORITHM                                        │
│                                                              │
│  Input: 40 protein values per cell                          │
│                                                              │
│  Example cell: CD4=high, CD8=low, CD45=high → "CD4 T cell"  │
│  Example cell: PanCK=high, CD45=low → "Tumor"               │
│  Example cell: CD68=high, CD11c=high → "Macrophage"         │
│                                                              │
│  Result: 16 cell type categories                            │
└─────────────────────────────────────────────────────────────┘

```

<div style="page-break-after: always;"></div>

## Each feature explained:

```
RAW DATA (2 million cells)
          │
          ├──────────────────────────────────────────────────────────┐
          │                                                          │
          ▼                                                          ▼
┌─────────────────────────┐                           ┌─────────────────────────┐
│  SPATIAL FEATURES       │                           │  NON-SPATIAL FEATURES   │
│  (THIS PAPER'S METHOD)  │                           │  (TRADITIONAL)          │
│                         │                           │                         │
│  1. neighbor_mat        │                           │  3. biomarker_region    │
│  2. ripley              │                           │  4. biomarker_cell      │
│                         │                           │  5. celltype_prop       │
└─────────────────────────┘                           └─────────────────────────┘
          │                                                          │
          └──────────────────────────┬───────────────────────────────┘
                                     │
                                     ▼
                          ┌─────────────────────────┐
                          │  SURVIVAL PREDICTION    │
                          │  (Random Survival       │
                          │   Forest)               │
                          └─────────────────────────┘

```

<div style="page-break-after: always;"></div>

### **Feature 1:Neighborhood Matrix (256 feature)**
> "For each cell type, what types of cells are its neighbors?"
```
Uses X,Y coordinates to find each cell's 10 nearest neighbors
Answers: "Which cell types are physically next to each other?"

Example insight:
"In patients who died quickly, Macrophages were next to Tumors"
"In patients who survived, T cells were next to Tumors"

```

example 16*16 matrix
```
              Neighbors →
            Tumor  Macro  Tcell  Bcell  ...
Center ↓    
Tumor       0.40   0.20   0.15   0.05   ...  (Tumors near Tumors)
Macrophage  0.35   0.25   0.20   0.05   ...  (Macros near Tumors)
T cell      0.20   0.15   0.30   0.10   ...  

```

### **feature 2: Ripley k function (80 feature)**

> "Are cells of the same type clustered together, spread apart, or randomly distributed?"
```
Uses X,Y coordinates to measure clustering at different distances
Answers: "Are high-expressing cells clustered together?"

Example insight:
"High CD21 expression is spatially clustered in bad-prognosis patients"
```
![alt text](image-6.png)

### **feature 3: Biomarker cell (40 feature)**
Protein intensity per cell type

### **feature 4: Biomarker region (40 feature)**
Protein intensity across entire tissue

### **Feature 5: cell type proportion ( 16 feature)**
Proportion of each cell type
<div style="page-break-after: always;"></div>

![alt text](image-5.png)



###  Comparison Methods

### **DenVar (Density Variation)
**
From another research paper — groups samples by how protein density varies spatially:
> Groups samples by how protein density varies across the tissue

```
File: denvar/{protein}_DenVar_clusters.csv

┌──────────┬─────────┬──────────────────────────────┐
│ Unnamed  │ cluster │           sample             │
├──────────┼─────────┼──────────────────────────────┤
│   100    │    1    │ UPMC_c004_v001_r001_reg009   │
│   106    │    2    │ UPMC_c004_v001_r001_reg015   │
│   141    │    2    │ UPMC_c004_v001_r001_reg060   │
└──────────┴─────────┴──────────────────────────────┘

Cluster = 1, 2, or 3 (which density group the sample belongs to)
One file per protein (40 proteins = 40 files)

```

### **Patwa et al. Methods**

**Method A: biomarker_threshold**
>For each sample, what % of cells are "positive" for each protein?
```
File: comparisons/biomarker_frac_positivity_qc_labeled.csv

Shape: (307 samples × 40 proteins)

Sample               │  CD31  │  CD57  │  CD4   │ ... │
─────────────────────┼────────┼────────┼────────┼─────┤
patient_001_region_A │  0.52  │  0.50  │  0.75  │ ... │  ← 52% of cells are CD31+
patient_002_region_A │  0.65  │  0.26  │  0.74  │ ... │
patient_003_region_A │  0.71  │  0.03  │  0.77  │ ... │

Value = fraction (0 to 1) of cells positive for that protein
```

**Method B: biomarker_interaction_counts**
>Count cells that are positive for BOTH protein A AND protein B
```
File: comparisons/interaction_biomarker_features.csv

Shape: (307 samples × 780 protein pairs)

Sample               │ CD31-CD31 │ CD31-CD57 │ CD31-CD4 │ ... │
─────────────────────┼───────────┼───────────┼──────────┼─────┤
patient_001_region_A │   4038    │   5656    │   9306   │ ... │
patient_002_region_A │   8856    │   6929    │  18448   │ ... │
patient_003_region_A │  11253    │   1063    │  23039   │ ... │

Value = COUNT of cells that are positive for BOTH proteins
40 proteins × 40 proteins = 780 unique pairs (including self-pairs)

```

<div style="page-break-after: always;"></div>

## **Survival analysis**
```
survival_status: 0 = Alive/Censored, 1 = Died
survival_day:    Number of days from diagnosis to death/last follow-up

```
## **train & evaluate**
```
┌────────────────────┐     ┌─────────────────┐     
│  Extract Feature   │ ──► │ + Survival      │ 
│  (different data)  │     │   Labels        │     
└────────────────────┘     └─────────────────┘    
         │
         ▼
┌─────────────────────────────────────┐
│  10-Fold Cross Validation           │
│  (Grouped by patient - no leakage)  │
└─────────────────────────────────────┘
         │
         ▼
┌─────────────────────────────────────┐
│  Random Survival Forest (RSF)       │
│  • 100 trees                        │
│  • Predicts: risk score (survival)  │
└─────────────────────────────────────┘
         │
         ▼
┌─────────────────────────────────────┐
│  Calculate:                         │
│  • C-Index (prediction accuracy)    │
│  • Feature importance (permutation) │
└─────────────────────────────────────┘
         │
         ▼
   Output: rsf_risk_scores/
```

<div style="page-break-after: always;"></div>

## **RSF Algorithm Overview**

```
RANDOM SURVIVAL FOREST
━━━━━━━━━━━━━━━━━━━━━━

SINGLE DECISION TREE:
                    Root
                     │
            ┌────────┴────────┐
      CD8_prop > 0.1?    CD8_prop ≤ 0.1?
            │                  │
       Lower Risk         Higher Risk

RANDOM FOREST:
  Tree 1: Uses features {A, C, F} → Risk₁
  Tree 2: Uses features {B, D, G} → Risk₂
  Tree 3: Uses features {A, E, H} → Risk₃
  ...
  Tree 100: Uses features {C, F, I} → Risk₁₀₀
  
  Final Risk = Average(Risk₁, Risk₂, ..., Risk₁₀₀)
```


**Data flow:**
```
Raw cell data → Feature extraction → Merge with labels → Train RSF → Evaluate
```

<div style="page-break-after: always;"></div>

## Results Interpretation

### C-Index Results

| Feature | C-Index | What It Means |
|---------|---------|---------------|
| celltype_prop | 0.705 | Cell composition alone is highly predictive |
| neighbor_mat | 0.704 | Spatial arrangement predicts survival well |
| biomarker_cell | 0.662 | Protein levels per cell = moderate |
| biomarker_region | 0.655 | Bulk protein levels = moderate |
| ripley | 0.528 | Clustering pattern alone = poor predictor |

<div style="page-break-after: always;"></div>

### What C-Index Means

```
C-Index = 1.0  → Perfect prediction (never happens)
C-Index = 0.7+ → Good (clinically useful)
C-Index = 0.5  → Random guess (useless)
```

### Key Findings

| Finding | Meaning |
|---------|---------|
| **neighbor_mat (0.704) ≈ celltype_prop (0.705)** | WHERE cells are located is as important as WHAT cells are present |
| **Combined (0.730) > Individual** | Combining features improves prediction |
| **Spatial > Traditional proteomics** | 0.704 > 0.655 means location beats protein levels |
| **Ripley K failed (0.528)** | Simple clustering statistics not enough |

<div style="page-break-after: always;"></div>

### Top Predictive Features (from neighbor_mat)

| Feature | Importance | Biological Meaning |
|---------|------------|-------------------|
| Tumor (Podo+) → Tumor (CD21+) | 0.012 | Tumor-tumor interactions matter |
| Tumor (Ki67+) → B cell | 0.009 | B cells near proliferating tumor = good? |
| Macrophage → Stromal | 0.007 | Immune-stromal interface important |
| Tumor (CD21+) → CD8 T cell | 0.005 | T cells infiltrating tumor = good |




