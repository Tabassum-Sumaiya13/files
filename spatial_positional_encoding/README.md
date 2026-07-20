# spatial_positional_encoding/ — Phase 2: feature verification

> **Name is historical.** Laplacian positional encoding was removed from the
> pipeline. This package now builds and verifies two feature blocks — spatial
> **enrichment** and celltype-conditioned **node markers** — on any cohort that
> `data_preprocessing/` has ingested. The retired PE code is in
> [`discarded/legacy_pipeline/`](../discarded/legacy_pipeline/).

Everything here reads the **canonical** cohort layout, so the same code runs on
UPMC, CRC, or any future dataset with no per-cohort branching:

```
data_preprocessing/datasets/<NAME>/processed/
    samples/sample_<acquisition_id>.parquet   X, Y, cluster_label, lineage, <markers>
    manifest.parquet                          patient_id, survival (NaN if none)
    marker_columns.txt
```

## Layout

```
spatial_positional_encoding/
├── run_verify.py         Phase 2  — enrichment block matrix (+ map perturbation)
├── run_verify_nodes.py   Phase 2b — node marker block matrix
├── run_survival.py       Phase 3  — OPTIONAL, cohorts with survival only
├── run_pipeline.py       convenience dispatcher / dataset lister
└── src/
    ├── cohort.py           load an ingested cohort, stream its samples
    ├── spatial_features.py Delaunay graph + the 5 enrichment scalars
    └── node_features.py    celltype-conditioned marker features
```

That is the whole active surface. If a module is not in this list it was retired —
see [`discarded/legacy_pipeline/README.md`](../discarded/legacy_pipeline/README.md)
for what each retired file was and what replaced it.

## Running it

```bash
# what is ingested, and which cohorts have survival
python run_pipeline.py --list

# enrichment block — the 5 spatial scalars
python run_verify.py --dataset CRC --taxonomy lineage
python run_verify.py --dataset CRC --taxonomy lineage --perturb-map

# node marker block — the 8 celltype-conditioned markers
python run_verify_nodes.py --dataset CRC

# optional survival downstream (exits cleanly on a survival-less cohort)
python run_survival.py --dataset UPMC
```

Outputs land in `<cohort>/processed/verification/`:

| file | contents |
|---|---|
| `verify_<taxonomy>.csv` | enrichment matrix — signal / specificity / stability / verdict |
| `verify_nodes_<taxonomy>.csv` | node marker matrix, plus per-feature support |
| `node_conditioning.csv` | which native cell types each node feature resolved to |
| `perturbation_<taxonomy>.csv` | the matrix re-run under each contested lineage assignment |
| `separability_<taxonomy>_<label>.csv` | optional classification check against a manifest label |
| `survival_<taxonomy>.csv` | optional RSF C-index comparison |

## The two blocks

**Enrichment** (`src/spatial_features.py`) — a parameter-free Delaunay graph
(mean degree ≈ 6) is built from X/Y, collapsed to a K×K neighbour-count matrix, and
read out as 5 abundance-corrected scalars: `kl_mean`, `kl_tumor`, `self_enrich`,
`immune_tumor`, `stroma_tumor`. Each is **0 under a random arrangement**, so a
non-zero value means organisation rather than abundance.

**Node markers** (`src/node_features.py`) — the mean of one protein marker over
*only* the cells of the type it is biologically read in. Conditioning is resolved
per cohort from `celltype_registry.csv` **by Cell Ontology term**, so `cd4_pd1`
becomes `['CD4 T cell']` on UPMC and `['CD4+ T cells', 'CD4+ T cells CD45RO+',
'CD4+ T cells GATA3+', 'Tregs']` on CRC with no hand mapping. Samples with no cells
of the conditioning type are reported as unsupported, never zero-filled.

## Taxonomy

`--taxonomy lineage` uses the 3-way `immune / tumour / stromal` column present in
every cohort — portable and directly comparable. `--taxonomy native` uses each
dataset's own `cluster_label` (16 for UPMC, 25 for CRC), which is finer but not
comparable across cohorts. For node features the flag only changes the
*composition baseline*; conditioning is always by CL term.

## Requirements

`pip install -r requirements.txt`. `run_survival.py` additionally needs
`scikit-survival`, which requires **scikit-learn ≤ 1.8** — the script pins
system-site packages ahead of user-site to satisfy that; see the note in
`discarded/legacy_pipeline/README.md` if the import fails.

## Results

Current numbers: [`doc/RESULT_REPORT.md`](../doc/RESULT_REPORT.md).
How it all works end to end, plus a manual verification walkthrough:
[`doc/PIPELINE_EXPLAINED.md`](../doc/PIPELINE_EXPLAINED.md).
