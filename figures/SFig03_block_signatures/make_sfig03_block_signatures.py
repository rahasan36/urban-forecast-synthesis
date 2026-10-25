"""Supplementary Fig. 3. Boundary density of each forecast relative to the block median.

For each 100 x 100 km block in which every forecast maps at least 10 km2 of urban land, a forecast's boundary
density (boundary length per unit urban area) is divided by the median density of the nine forecasts in that block.
a-c  share of blocks where the forecast lies below (blue) or above (red) the block median, dark shading below
     0.8 or above 1.25 times the median; SSP1 2020 (a), SSP1 2100 (b) and SSP5 2100 (c).
d    distribution of the ratio to the block median by forecast, SSP1 2020 (solid) and 2100 (hatched).
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
FIG_DATA = os.path.join(REPO, "data", "figures")
NAME = "SFig03_block_signatures"

FORECASTS = ["Li2021", "Chen2020", "Chen2022", "Gao2021", "Hou2022", "Koomen2UP", "He2023", "Zhang2023", "Zhao2025"]
LABEL = {"Li2021": "Li 2021", "Chen2020": "Chen 2020", "Chen2022": "Chen 2022", "Gao2021": "Gao 2021",
         "Hou2022": "Hou 2022", "Koomen2UP": "2UP 2023", "He2023": "He 2023", "Zhang2023": "Zhang 2023",
         "Zhao2025": "FU3D 2025"}
CASES = [("SSP1", 2020), ("SSP1", 2100), ("SSP5", 2100)]
DARK_BLUE, LIGHT_BLUE, LIGHT_RED, DARK_RED = "#2B5C9E", "#A9C4E4", "#F2B3C6", "#B0244F"

plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
                     "font.size": 8, "pdf.fonttype": 42, "axes.linewidth": 0.6, "hatch.linewidth": 0.5})


def out_path(folder, ext):
    path = os.path.join(folder, f"{NAME}.{ext}")
    k = 1
    while os.path.exists(path):
        path = os.path.join(folder, f"{NAME}_r{k}.{ext}")
        k += 1
    return path


def density_ratios(blocks_dir, ssp, year):
    """Boundary density of each forecast divided by the block median, for blocks with >= 10 km2 in every forecast."""
    z = np.load(os.path.join(blocks_dir, f"blocks100km_{ssp}_{year}.npz"))
    area = z["area"].astype(float)
    boundary = z["bnd"].astype(float)
    eligible = area.min(0) >= 10
    density = boundary[:, eligible] / area[:, eligible]
    return density / np.median(density, 0)


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--blocks", default=os.path.join(FIG_DATA, "blocks100km"),
                    help="folder with blocks100km_<SSP>_<year>.npz")
    ap.add_argument("--out", default=HERE)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    ratios = {c: density_ratios(a.blocks, *c) for c in CASES}
    order = list(np.argsort([np.median(ratios[("SSP1", 2020)][k]) for k in range(9)]))

    fig = plt.figure(figsize=(7.2, 5.6))
    width = 0.235
    for i, c in enumerate(CASES):
        ax = fig.add_axes([0.11 + i * (width + 0.055), 0.56, width, 0.39])
        r = ratios[c]
        for row, k in enumerate(order):
            far_below = 100 * (r[k] < 0.8).mean()
            below = 100 * (r[k] < 1).mean()
            far_above = 100 * (r[k] > 1.25).mean()
            above = 100 * (r[k] > 1).mean()
            ax.barh(row, -below, color=LIGHT_BLUE, height=0.72, lw=0)
            ax.barh(row, -far_below, color=DARK_BLUE, height=0.72, lw=0)
            ax.barh(row, above, color=LIGHT_RED, height=0.72, lw=0)
            ax.barh(row, far_above, color=DARK_RED, height=0.72, lw=0)
        ax.axvline(0, color="black", lw=0.6)
        ax.set_xlim(-100, 100)
        ax.set_xticks([-100, -50, 0, 50, 100])
        ax.set_xticklabels(["100", "50", "0", "50", "100"])
        ax.set_yticks(np.arange(9))
        ax.set_yticklabels([LABEL[FORECASTS[k]] for k in order] if i == 0 else [])
        ax.set_ylim(-0.6, 8.6)
        ax.spines[["top", "right"]].set_visible(False)
        ax.tick_params(length=2)
        ax.set_xlabel("Blocks (%)")
        ax.text(0.5, 1.02, f"{c[0]}, {c[1]}", transform=ax.transAxes, ha="center", va="bottom")
        x_letter = 0.11 + i * (width + 0.055) - (0.1 if i == 0 else 0.035)
        fig.text(x_letter, 0.985, "abc"[i], fontsize=10, fontweight="bold", va="top")
    fig.legend(handles=[Patch(color=DARK_BLUE, label="Below 0.8 × block median"),
                        Patch(color=LIGHT_BLUE, label="Below block median"),
                        Patch(color=LIGHT_RED, label="Above block median"),
                        Patch(color=DARK_RED, label="Above 1.25 × block median")],
               loc="upper center", bbox_to_anchor=(0.55, 0.465), ncol=4, frameon=False, handlelength=1.0,
               columnspacing=1.0)

    axd = fig.add_axes([0.11, 0.09, 0.87, 0.28])
    pos = np.arange(9)
    for c, hatch, off in ((("SSP1", 2020), "", -0.18), (("SSP1", 2100), "////", 0.18)):
        data = [np.log2(ratios[c][k]) for k in order]
        bp = axd.boxplot(data, positions=pos + off, widths=0.3, patch_artist=True, showfliers=False, whis=(5, 95),
                         medianprops=dict(color="black", lw=0.7), boxprops=dict(lw=0.5),
                         whiskerprops=dict(lw=0.5), capprops=dict(lw=0.5))
        for box, k in zip(bp["boxes"], order):
            m = np.median(ratios[c][k])
            box.set_facecolor(DARK_RED if m > 1.05 else (DARK_BLUE if m < 0.95 else "#D9D9D9"))
            box.set_hatch(hatch)
            box.set_edgecolor("#333333")
            box.set_alpha(0.9)
    axd.axhline(0, color="black", lw=0.6, ls=(0, (2, 2)))
    axd.set_xticks(pos)
    axd.set_xticklabels([LABEL[FORECASTS[k]] for k in order])
    axd.set_xlim(-0.6, 8.6)
    axd.set_yticks([-2, -1, 0, 1])
    axd.set_yticklabels(["0.25×", "0.5×", "1×", "2×"])
    axd.set_ylim(-2.3, 1.4)
    axd.set_ylabel("Boundary density relative\nto block median")
    axd.spines[["top", "right"]].set_visible(False)
    axd.tick_params(length=2)
    axd.legend(handles=[Patch(facecolor="white", edgecolor="#333333", label="2020"),
                        Patch(facecolor="white", edgecolor="#333333", hatch="////", label="2100")],
               loc="lower right", frameon=False, ncol=2, handlelength=1.4, title="SSP1", title_fontsize=8)
    fig.text(0.0, 0.4, "d", fontsize=10, fontweight="bold", va="top")

    for ext in ("png", "pdf"):
        fig.savefig(out_path(a.out, ext), dpi=300, facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    main()
