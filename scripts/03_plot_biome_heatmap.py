#!/usr/bin/env python3
"""
Panel B - biome distribution of the 43 cargo-rich plasmids (MGnify Branchwater,
containment >= 0.90). Rows = plasmids (sorted by cross-biome breadth), columns = biome
categories, cell = number of near-complete detections (log colour scale). The two
plasmids selected for the deep-dive have their labels shown in bold.

Input : P2_plasmid_biome_matrix.csv  (same directory)
Output: Figure1B_biome_heatmap.pdf   (vector, Arial, editable text)

Run:  python plot_biome_heatmap.py
Deps: pandas, numpy, matplotlib
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

# ---------------------------------------------------------------- config
HERE = os.path.dirname(os.path.abspath(__file__))
INFILE  = os.path.join(HERE, "P2_plasmid_biome_matrix_bioproject.csv")
OUTFILE = os.path.join(HERE, "Figure1B_biome_heatmap.pdf")

# Arial to match the iTOL tree panel; Liberation Sans is a metric-compatible
# fallback on Linux, DejaVu as last resort. Keep text as editable TrueType in the PDF.
matplotlib.rcParams["font.family"]   = ["Arial", "Liberation Sans", "DejaVu Sans"]
matplotlib.rcParams["pdf.fonttype"]  = 42
matplotlib.rcParams["ps.fonttype"]   = 42
matplotlib.rcParams["svg.fonttype"]  = "none"

# biome columns to plot, in display order, with pretty labels
BIOMES = ["soil","plant/rhizosphere","food/silage","host-associated","wastewater",
          "compost/digester","freshwater","sediment","marine","peat/wetland",
          "lichen/moss","other","unclassified"]
PRETTY = {b: b[0].upper()+b[1:] for b in BIOMES}

CMAP = "YlGnBu"          # sequential, colour-blind friendly (light = few, dark = many)
ZERO_COLOR = "#ffffff"   # absent (0 detections) drawn white
SHOW_GROUP_LABELS   = True   # annotate the "recovered elsewhere" vs "not recovered" split
DROP_EMPTY_BIOMES   = True   # hide biome columns with no detection in any plasmid
SHOW_ONLY_TRAVELLERS = True  # True: show only the 16 plasmids recovered elsewhere,
                             # collapsing the 27 "stay-at-home" plasmids into a caption note

# ---------------------------------------------------------------- load
df = pd.read_csv(INFILE)

# optionally show only plasmids recovered near-complete beyond their source genome
n_all = len(df)
if SHOW_ONLY_TRAVELLERS:
    df = df[df["total_near_complete"] > 0].reset_index(drop=True)
n_hidden = n_all - len(df)

# optionally drop biome columns that are empty across all displayed plasmids
biomes = list(BIOMES)
if DROP_EMPTY_BIOMES:
    biomes = [b for b in biomes if df[b].sum() > 0]

# keep the file's existing order (already sorted by breadth: broadest at top)
labels   = df["plasmid"].tolist()
selected = (df["selected"].astype(str).str.lower() == "yes").tolist()
mat = df[biomes].to_numpy(dtype=float)
n_rows, n_cols = mat.shape

# split index between plasmids recovered elsewhere (total>0) and the rest
n_travellers = int((df["total_near_complete"] > 0).sum())

# ---------------------------------------------------------------- figure
masked = np.ma.masked_where(mat <= 0, mat)          # zeros -> masked (white)
vmax = mat.max()
norm = LogNorm(vmin=1, vmax=vmax)
cmap = matplotlib.colormaps[CMAP].copy()
cmap.set_bad(ZERO_COLOR)

fig_h = 0.24 * n_rows + 2.2
fig_w = 0.52 * n_cols + 3.2
fig, ax = plt.subplots(figsize=(fig_w, fig_h))

mesh = ax.pcolormesh(masked, cmap=cmap, norm=norm,
                     edgecolors="#d9d9d9", linewidth=0.4)

# annotate non-zero cells with the count
log_vmax = np.log10(vmax)
for i in range(n_rows):
    for j in range(n_cols):
        v = mat[i, j]
        if v > 0:
            frac = np.log10(v) / log_vmax if log_vmax > 0 else 1.0
            ax.text(j + 0.5, i + 0.5, f"{int(v)}",
                    ha="center", va="center", fontsize=11,
                    color="white" if frac > 0.55 else "#222222")

# axes / ticks
ax.set_xticks(np.arange(n_cols) + 0.5)
ax.set_xticklabels([PRETTY[b] for b in biomes], rotation=40, ha="right", fontsize=14)
ax.set_yticks(np.arange(n_rows) + 0.5)
ax.set_yticklabels(labels, fontsize=13)
for tick, sel in zip(ax.get_yticklabels(), selected):
    if sel:
        tick.set_fontweight("bold")     # highlight selected plasmids
ax.invert_yaxis()                       # broadest plasmid at the top
ax.set_xlim(0, n_cols); ax.set_ylim(n_rows, 0)
ax.tick_params(length=0)
for s in ax.spines.values():
    s.set_visible(False)

# divider + group labels
if SHOW_GROUP_LABELS and 0 < n_travellers < n_rows:
    ax.axhline(n_travellers, color="#444444", linewidth=1.1)
    ax.text(n_cols + 0.15, n_travellers / 2, f"recovered elsewhere\n(n={n_travellers})",
            rotation=90, va="center", ha="left", fontsize=13, color="#444444")
    ax.text(n_cols + 0.15, (n_travellers + n_rows) / 2,
            f"not recovered\nbeyond source (n={n_rows - n_travellers})",
            rotation=90, va="center", ha="left", fontsize=13, color="#444444")

# colourbar
cbar = fig.colorbar(mesh, ax=ax, fraction=0.025, pad=0.14,
                    ticks=[1, 3, 10, 30, 100])
cbar.ax.set_yticklabels(["1", "3", "10", "30", "100"], fontsize=14)
cbar.set_label("Distinct BioProjects (n, log scale)", fontsize=15)
cbar.outline.set_visible(False)

ax.set_title("Biome distribution of cargo-rich plasmids\n(distinct BioProjects, containment ≥ 0.90)",
             fontsize=18, pad=8)

if SHOW_ONLY_TRAVELLERS and n_hidden:
    fig.text(0.01, 0.005,
             f"{n_hidden} of {n_all} plasmids were not recovered near-complete beyond their "
             f"source genome and are not shown.", fontsize=13, color="#444444")
fig.tight_layout()
fig.savefig(OUTFILE, bbox_inches="tight")
print("wrote", OUTFILE)
