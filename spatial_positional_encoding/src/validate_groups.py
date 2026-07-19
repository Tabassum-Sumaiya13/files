"""Biologically validate the IMMUNE / TUMOR / STROMA grouping against the
protein data. Does each celltype actually express its lineage's marker?

Canonical lineage markers (textbook, not chosen from this data):
  CD45       pan-leukocyte      -> every immune cell must be CD45-high
  PanCK      pan-cytokeratin    -> every epithelial/tumour cell must be PanCK-high
  Vimentin   mesenchymal        -> stromal/fibroblast
  aSMA       smooth muscle      -> vessel walls, myofibroblasts
  CD31       endothelial        -> blood vessels
  Podoplanin lymphatic endothel -> lymph vessels
"""
import numpy as np
import pandas as pd
from pathlib import Path

RAW = Path(r"d:/Desktop/FYDP/FYDP 3/files/spatial_positional_encoding/data/raw")
OUT = Path(r"C:/Users/User/AppData/Local/Temp/claude/d--Desktop-FYDP-FYDP-3-files/08a1fd2b-c50e-46b3-ad94-e7abe111ce07/scratchpad")

LINEAGE = ["CD45", "PanCK", "Vimentin", "aSMA", "CD31", "Podoplanin",
           "CD3e", "CD20", "CD21", "CD68"]
COLS = LINEAGE + ["cluster_label"]

print("Reading 2M cells in chunks (only lineage markers)...")
sums, counts = {}, {}
for chunk in pd.read_csv(RAW / "labeled_arcsinh_norm_data.csv",
                         usecols=COLS, chunksize=250_000):
    g = chunk.groupby("cluster_label")
    s = g[LINEAGE].sum()
    n = g.size()
    for lab in s.index:
        sums[lab] = sums.get(lab, 0) + s.loc[lab]
        counts[lab] = counts.get(lab, 0) + n.loc[lab]

mean = pd.DataFrame({lab: sums[lab] / counts[lab] for lab in sums}).T
mean.index.name = "cluster_label"
n_cells = pd.Series(counts, name="n_cells")
print(f"  {int(n_cells.sum()):,} cells across {len(mean)} celltypes\n")

# z-score each marker ACROSS celltypes: "is this type high for this marker
# relative to the other types?"
z = (mean - mean.mean()) / mean.std()

MY_GROUP = {
    "APC": "IMMUNE", "B cell": "IMMUNE", "CD4 T cell": "IMMUNE",
    "CD8 T cell": "IMMUNE", "Granulocyte": "IMMUNE", "Macrophage": "IMMUNE",
    "Naive immune cell": "IMMUNE",
    "Tumor": "TUMOR", "Tumor (CD15+)": "TUMOR", "Tumor (CD20+)": "TUMOR",
    "Tumor (CD21+)": "TUMOR", "Tumor (Ki67+)": "TUMOR", "Tumor (Podo+)": "TUMOR",
    "Lymph vessel": "STROMA", "Stromal / Fibroblast": "STROMA", "Vessel": "STROMA",
}
EXPECT = {"IMMUNE": "CD45", "TUMOR": "PanCK", "STROMA": None}

print("=" * 96)
print("PER-CELLTYPE LINEAGE MARKER EXPRESSION (z-scored across celltypes)")
print("=" * 96)
hdr = f"{'celltype':22s}{'n':>9s}{'group':>8s}" + "".join(f"{m[:8]:>10s}" for m in LINEAGE[:6])
print(hdr)
print("-" * 96)
order = sorted(mean.index, key=lambda x: (MY_GROUP.get(x, "ZZ"), x))
for lab in order:
    row = "".join(f"{z.loc[lab, m]:10.2f}" for m in LINEAGE[:6])
    print(f"{lab:22s}{int(n_cells[lab]):9,d}{MY_GROUP.get(lab,'?'):>8s}" + row)

print("\n" + "=" * 96)
print("VERDICT — does each type express its lineage marker above the cross-type mean?")
print("=" * 96)
for lab in order:
    g = MY_GROUP.get(lab, "?")
    marker = EXPECT.get(g)
    if marker is None:
        # stroma: any of the structural markers
        best = max(["Vimentin", "aSMA", "CD31", "Podoplanin"],
                   key=lambda m: z.loc[lab, m])
        ok = z.loc[lab, best] > 0
        detail = f"best structural = {best} (z={z.loc[lab, best]:+.2f})"
    else:
        ok = z.loc[lab, marker] > 0
        detail = f"{marker} z={z.loc[lab, marker]:+.2f}"
    # cross-check: is it high for the WRONG lineage?
    conflict = ""
    if g == "TUMOR" and z.loc[lab, "CD45"] > z.loc[lab, "PanCK"]:
        conflict = "  <-- CD45 > PanCK: looks IMMUNE, not tumour!"
    if g == "IMMUNE" and z.loc[lab, "PanCK"] > z.loc[lab, "CD45"]:
        conflict = "  <-- PanCK > CD45: looks EPITHELIAL, not immune!"
    print(f"  {'PASS' if ok else 'FAIL'}  {lab:22s} {g:7s} {detail}{conflict}")

mean.to_csv(OUT / "lineage_mean.csv")
z.to_csv(OUT / "lineage_z.csv")
n_cells.to_csv(OUT / "celltype_counts.csv")
print(f"\nsaved lineage_mean.csv / lineage_z.csv / celltype_counts.csv")
