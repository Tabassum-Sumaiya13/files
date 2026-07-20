"""Spatial maps: protein variation, celltypes, and the 3-group collapse.

Colors come from the dataviz reference palette (validated: 3-slot categorical
passes all-pairs CVD in light mode; magenta carries a contrast WARN -> relief is
the always-visible labelled legend).
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, to_rgb, rgb_to_hsv, hsv_to_rgb
from matplotlib.lines import Line2D
from pathlib import Path

ROOT = Path(r"d:/Desktop/FYDP/FYDP 3/files/spatial_positional_encoding")
OUT = ROOT / "data/processed/outputs/spatial_maps"
OUT.mkdir(parents=True, exist_ok=True)
UM_PER_PX = 0.377

# --- palette (from dataviz references/palette.md) ---
SURFACE, INK, INK2, MUTED = "#fcfcfb", "#0b0b0b", "#52514e", "#898781"
BLUE_RAMP = ["#cde2fb", "#b7d3f6", "#9ec5f4", "#86b6ef", "#6da7ec", "#5598e7",
             "#3987e5", "#2a78d6", "#256abf", "#1c5cab", "#184f95", "#104281", "#0d366b"]
SEQ = LinearSegmentedColormap.from_list("blue_seq", BLUE_RAMP)
GROUP_COLOR = {"TUMOR": "#2a78d6", "IMMUNE": "#008300", "STROMA": "#e87ba4"}

GROUPS = {
    "TUMOR":  ["Tumor", "Tumor (CD15+)", "Tumor (CD20+)", "Tumor (CD21+)",
               "Tumor (Ki67+)", "Tumor (Podo+)"],
    "IMMUNE": ["APC", "B cell", "CD4 T cell", "CD8 T cell", "Granulocyte",
               "Macrophage", "Naive immune cell"],
    "STROMA": ["Lymph vessel", "Stromal / Fibroblast", "Vessel"],
}
LABEL_OF_ID = ["APC", "B cell", "CD4 T cell", "CD8 T cell", "Granulocyte",
               "Lymph vessel", "Macrophage", "Naive immune cell",
               "Stromal / Fibroblast", "Tumor", "Tumor (CD15+)", "Tumor (CD20+)",
               "Tumor (CD21+)", "Tumor (Ki67+)", "Tumor (Podo+)", "Vessel"]


def shades(base_hex, n):
    """n lightness steps of ONE hue -> composite encoding (hue=group, value=type).

    Interpolates tint(white) -> base -> shade(black) so the steps actually
    separate; an HSV value-only sweep left adjacent tumour types identical.
    """
    base = np.array(to_rgb(base_hex))
    light = base + (1.0 - base) * 0.66     # tint toward white
    dark = base * 0.38                     # shade toward black
    ramp = LinearSegmentedColormap.from_list("grp", [light, base, dark])
    return [ramp(x)[:3] for x in np.linspace(0.04, 0.96, n)]


CT_COLOR = {}
for g, members in GROUPS.items():
    for lab, c in zip(members, shades(GROUP_COLOR[g], len(members))):
        CT_COLOR[lab] = c

# --- pick 3 samples spanning immune_tumor (the feature the maps should show) ---
enr = pd.read_parquet(ROOT / "data/processed/neighbor_features/enrichment.parquet")
srt = enr.sort_values("immune_tumor")
SAMPLES = [srt.index[0], srt.index[len(srt) // 2], srt.index[-1]]
TAGS = ["most excluded", "median", "least excluded"]
print("Samples chosen by immune_tumor:")
for s, t in zip(SAMPLES, TAGS):
    print(f"  {s:32s} immune_tumor = {enr.loc[s,'immune_tumor']:+.3f}  ({t})")

# --- load coords + labels ---
loc = pd.read_csv(ROOT / "data/raw/cell_locations_and_labels.csv",
                  usecols=["ACQUISITION_ID", "CELL_ID", "X", "Y", "CLUSTER_ID"])
loc = loc[loc.ACQUISITION_ID.isin(SAMPLES)]

# --- load protein for those samples only ---
MARKERS = ["PanCK", "CD45", "Vimentin"]
MARKER_DESC = {"PanCK": "tumour / epithelial", "CD45": "immune (pan-leukocyte)",
               "Vimentin": "stromal / mesenchymal"}
chunks = []
for ch in pd.read_csv(ROOT / "data/raw/labeled_arcsinh_norm_data.csv",
                      usecols=MARKERS + ["sample_id", "cell_id"], chunksize=250_000):
    ch = ch[ch.sample_id.isin(SAMPLES)]
    if len(ch):
        chunks.append(ch)
expr = pd.concat(chunks)
print(f"\nexpression rows for these samples: {len(expr):,}")

# MERGE on (sample, cell) — the two files are not row-aligned, so positional
# indexing would silently paint the wrong cell's protein onto a coordinate.
cells = loc.merge(expr, left_on=["ACQUISITION_ID", "CELL_ID"],
                  right_on=["sample_id", "cell_id"], how="inner")
print(f"merged coords+protein: {len(cells):,} cells "
      f"(loc={len(loc):,}, expr={len(expr):,})")
assert len(cells) > 0, "merge produced no rows — check key dtypes"
for sid in SAMPLES:
    a, b = (loc.ACQUISITION_ID == sid).sum(), (cells.ACQUISITION_ID == sid).sum()
    print(f"  {sid.replace('UPMC_',''):26s} loc={a:6,d}  merged={b:6,d}"
          f"{'  <-- LOSS' if b < a else ''}")


def frame(ax, title, sub=None):
    ax.set_facecolor(SURFACE)
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_color("#e1e0d9"); sp.set_linewidth(0.8)
    ax.set_title(title, fontsize=10, color=INK, pad=4)
    # xlabel, not free text: it is inside the layout engine so it cannot clip.
    if sub:
        ax.set_xlabel(sub, fontsize=8, color=MUTED, labelpad=6)


def scalebar(ax, xs, ys):
    L = 200 / UM_PER_PX  # 200 um
    x0 = xs.min() + (xs.max() - xs.min()) * 0.05
    y0 = ys.min() + (ys.max() - ys.min()) * 0.05
    ax.plot([x0, x0 + L], [y0, y0], color=INK2, lw=1.6, solid_capstyle="butt")
    ax.text(x0 + L / 2, y0 + (ys.max() - ys.min()) * 0.02, "200 µm",
            ha="center", va="bottom", fontsize=6.5, color=INK2)


# ---------------------------------------------------------------- FIGURE 1
fig, axes = plt.subplots(len(SAMPLES), len(MARKERS),
                         figsize=(4.1 * len(MARKERS), 4.0 * len(SAMPLES)))
fig.patch.set_facecolor(SURFACE)
for r, (sid, tag) in enumerate(zip(SAMPLES, TAGS)):
    s = cells[cells.ACQUISITION_ID == sid]
    for c, m in enumerate(MARKERS):
        ax = axes[r, c]
        v = s[m].values
        lo, hi = np.percentile(v, [1, 99])
        sc = ax.scatter(s["X"].values, s["Y"].values, c=v, cmap=SEQ,
                        vmin=lo, vmax=hi, s=1.4, linewidths=0, rasterized=True)
        ax.invert_yaxis(); ax.set_aspect("equal")
        frame(ax, f"{m}  ·  {MARKER_DESC[m]}" if r == 0 else m)
        if c == 0:
            ax.set_ylabel(f"{sid.replace('UPMC_','')}\n{tag}", fontsize=8, color=INK2)
        if c == len(MARKERS) - 1:
            cb = fig.colorbar(sc, ax=ax, fraction=0.042, pad=0.02)
            cb.set_label("arcsinh intensity", fontsize=7, color=INK2)
            cb.ax.tick_params(labelsize=6, colors=MUTED)
            cb.outline.set_visible(False)
        if r == len(SAMPLES) - 1 and c == 0:
            scalebar(ax, s["X"].values, s["Y"].values)
fig.suptitle("Protein variation in situ — raw cell coordinates, one dot per cell",
             fontsize=14, color=INK, y=1.005)
fig.text(0.5, 0.986, "Sequential single-hue ramp: light = low, dark = high. "
         "Scaled 1st–99th percentile per panel.", ha="center", fontsize=9, color=MUTED)
fig.tight_layout(rect=[0, 0, 1, 0.978])
fig.savefig(OUT / "01_protein_variation.png", dpi=170, facecolor=SURFACE,
            bbox_inches="tight")
plt.close(fig)
print("  -> 01_protein_variation.png")

# ---------------------------------------------------------------- FIGURE 2
fig, axes = plt.subplots(1, len(SAMPLES), figsize=(5.0 * len(SAMPLES), 5.4))
fig.patch.set_facecolor(SURFACE)
for ax, sid, tag in zip(np.atleast_1d(axes), SAMPLES, TAGS):
    s = loc[loc.ACQUISITION_ID == sid]
    cols = [CT_COLOR[LABEL_OF_ID[i]] for i in s["CLUSTER_ID"].values]
    ax.scatter(s["X"], s["Y"], c=cols, s=1.6, linewidths=0, rasterized=True)
    ax.invert_yaxis(); ax.set_aspect("equal")
    frame(ax, f"{sid.replace('UPMC_','')}", f"{tag}  ·  {len(s):,} cells")
handles = []
for gi, (g, members) in enumerate(GROUPS.items()):
    if gi:
        handles.append(Line2D([], [], ls="", marker="", label=" "))
    handles.append(Line2D([], [], ls="", marker="", label=f"{g}"))
    for lab in members:
        handles.append(Line2D([], [], marker="o", ls="", markersize=7,
                              markerfacecolor=CT_COLOR[lab],
                              markeredgecolor="none", label=f"   {lab}"))
leg = fig.legend(handles=handles, loc="center left", bbox_to_anchor=(1.005, 0.5),
                 frameon=False, fontsize=8.5, labelcolor=INK2, ncol=1,
                 handletextpad=0.6, labelspacing=0.55)
fig.suptitle("All 16 cell types — hue = biological group, shade = type within group",
             fontsize=14, color=INK, y=1.02)
fig.text(0.5, 0.965, "Composite encoding: 16 categories exceed any CVD-safe "
         "categorical palette, so hue carries the group and lightness the type. "
         "Read the group reliably; exact type only within a group.",
         ha="center", fontsize=9, color=MUTED)
fig.tight_layout(rect=[0, 0, 0.84, 0.955])
fig.savefig(OUT / "02_celltypes.png", dpi=170, facecolor=SURFACE,
            bbox_inches="tight")
plt.close(fig)
print("  -> 02_celltypes.png")

# ---------------------------------------------------------------- FIGURE 3
GROUP_OF = {lab: g for g, mem in GROUPS.items() for lab in mem}
fig, axes = plt.subplots(1, len(SAMPLES), figsize=(5.0 * len(SAMPLES), 5.4))
fig.patch.set_facecolor(SURFACE)
for ax, sid, tag in zip(np.atleast_1d(axes), SAMPLES, TAGS):
    s = loc[loc.ACQUISITION_ID == sid]
    g = [GROUP_OF[LABEL_OF_ID[i]] for i in s["CLUSTER_ID"].values]
    ax.scatter(s["X"], s["Y"], c=[GROUP_COLOR[x] for x in g], s=1.6,
               linewidths=0, rasterized=True)
    ax.invert_yaxis(); ax.set_aspect("equal")
    it = enr.loc[sid, "immune_tumor"]
    frame(ax, f"{sid.replace('UPMC_','')}",
          f"{tag}  ·  immune_tumor = {it:+.2f}  ({2**it:.2f}× chance)")
handles = [Line2D([], [], marker="o", ls="", markersize=9,
                  markerfacecolor=GROUP_COLOR[g], markeredgecolor="none", label=g)
           for g in ["TUMOR", "IMMUNE", "STROMA"]]
fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, -0.02),
           frameon=False, fontsize=10, labelcolor=INK2, ncol=3)
fig.suptitle("The 3-group collapse — what immune_tumor actually measures",
             fontsize=14, color=INK, y=1.015)
fig.text(0.5, 0.962, "Left: immune (green) walled out of tumour (blue) nests.  "
         "Right: the same three groups fully intermixed.",
         ha="center", fontsize=8.5, color=MUTED)
fig.tight_layout(rect=[0, 0.05, 1, 0.952])
fig.savefig(OUT / "03_groups.png", dpi=170, facecolor=SURFACE, bbox_inches="tight")
plt.close(fig)
print("  -> 03_groups.png")

# ---------------------------------------------------------------- FIGURE 4
# Biological validation. z-scores are polar (above/below the cross-type mean),
# so this is a DIVERGING scale: two hues + neutral gray midpoint, never a rainbow.
SCRATCH = Path(r"C:/Users/User/AppData/Local/Temp/claude"
               r"/d--Desktop-FYDP-FYDP-3-files/08a1fd2b-c50e-46b3-ad94-e7abe111ce07"
               r"/scratchpad")
z = pd.read_csv(SCRATCH / "lineage_z.csv", index_col=0)
nc = pd.read_csv(SCRATCH / "celltype_counts.csv", index_col=0).iloc[:, 0]

DIV = LinearSegmentedColormap.from_list(
    "red_gray_blue", ["#d03b3b", "#e89a9a", "#f0efec", "#86b6ef", "#0d366b"])

MARKS = ["CD45", "PanCK", "Vimentin", "aSMA", "CD31", "Podoplanin"]
EXPECTED = {"IMMUNE": "CD45", "TUMOR": "PanCK"}
order, rowgroup = [], []
for g in ["TUMOR", "IMMUNE", "STROMA"]:
    for lab in GROUPS[g]:
        order.append(lab); rowgroup.append(g)

Z = z.loc[order, MARKS].values
fig, ax = plt.subplots(figsize=(9.2, 7.4))
fig.patch.set_facecolor(SURFACE)
im = ax.imshow(Z, cmap=DIV, vmin=-2.6, vmax=2.6, aspect="auto")
ax.set_xticks(range(len(MARKS)))
ax.set_xticklabels(MARKS, fontsize=9, color=INK2)
ax.set_yticks(range(len(order)))
ax.set_yticklabels([f"{l}  ({nc[l]:,})" for l in order], fontsize=8.5, color=INK2)
ax.tick_params(length=0)
for sp in ax.spines.values():
    sp.set_visible(False)

for i in range(len(order)):
    for j in range(len(MARKS)):
        v = Z[i, j]
        ax.text(j, i, f"{v:+.2f}", ha="center", va="center", fontsize=7.5,
                color="#ffffff" if abs(v) > 1.5 else INK)

# ring the cell each row is REQUIRED to be high in
for i, (lab, g) in enumerate(zip(order, rowgroup)):
    m = EXPECTED.get(g)
    if m:
        j = MARKS.index(m)
        ok = z.loc[lab, m] > 0
        ax.add_patch(plt.Rectangle((j - .5, i - .5), 1, 1, fill=False,
                                   edgecolor="#0ca30c" if ok else "#d03b3b", lw=2.4))

# group separators
for b in [len(GROUPS["TUMOR"]) - .5,
          len(GROUPS["TUMOR"]) + len(GROUPS["IMMUNE"]) - .5]:
    ax.axhline(b, color=INK, lw=1.4)
for g, y in [("TUMOR", 2.5), ("IMMUNE", 9), ("STROMA", 14)]:
    ax.text(-2.45, y, g, fontsize=10, color=GROUP_COLOR[g], rotation=90,
            va="center", ha="center", weight="bold")

cb = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
cb.set_label("z-score across the 16 celltypes  (0 = cross-type average)",
             fontsize=8, color=INK2)
cb.ax.tick_params(labelsize=7, colors=MUTED); cb.outline.set_visible(False)

ax.set_title("Biological validation — does each celltype express its lineage marker?",
             fontsize=13, color=INK, pad=26)
ax.text(0.5, 1.035, "Green ring = passes its required marker · red ring = fails.  "
        "Immune must be CD45-high; tumour must be PanCK-high.",
        transform=ax.transAxes, ha="center", fontsize=8.5, color=MUTED)
fig.tight_layout()
fig.savefig(OUT / "04_lineage_validation.png", dpi=170, facecolor=SURFACE,
            bbox_inches="tight")
plt.close(fig)
print("  -> 04_lineage_validation.png")

print(f"\nSaved to {OUT}")
