"""Supplementary Fig. 5. Urban boundaries adjoining cropland and forest across scenarios.

Global length of urban boundary in contact with cropland (A-C) and forest (D-F) of the 2016 land-cover
baseline, in 2020, 2050 and 2100, for the nine forecasts. Columns show SSP1, SSP3 and SSP5.
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
FIG_DATA = os.path.join(REPO, "data", "figures")
NAME = "SFig05_contact_over_time"

FORECASTS = ["Li2021", "Chen2020", "Chen2022", "Gao2021", "Hou2022", "Koomen2UP", "He2023", "Zhang2023", "Zhao2025"]
LABELS = ["Li 2021", "Chen 2020", "Chen 2022", "Gao 2021", "Hou 2022", "2UP 2023", "He 2023", "Zhang 2023", "FU3D 2025"]
COLORS = ["#0077BB", "#EE7733", "#DDAA00", "#009988", "#CC3311", "#33BBEE", "#AA3377", "#777777", "#222222"]
COLOR = dict(zip(FORECASTS, COLORS))
PATHWAYS = [1, 3, 5]
YEARS = [2020, 2050, 2100]

plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
    "font.size": 13, "axes.labelsize": 13, "axes.titlesize": 14, "xtick.labelsize": 12, "ytick.labelsize": 12,
    "legend.fontsize": 12, "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": 0.7,
    "pdf.fonttype": 42, "ps.fonttype": 42, "savefig.facecolor": "white",
})


def out_path(folder, ext):
    path = os.path.join(folder, f"{NAME}.{ext}")
    k = 1
    while os.path.exists(path):
        path = os.path.join(folder, f"{NAME}_r{k}.{ext}")
        k += 1
    return path


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--landcover_metrics", default=os.path.join(FIG_DATA, "results", "landcover_metrics.csv"))
    ap.add_argument("--out", default=HERE)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    lc = pd.read_csv(a.landcover_metrics)
    lc = lc[lc.threshold.eq(0.4) & lc.region.eq("Global")]

    fig, axs = plt.subplots(2, 3, figsize=(11.6, 9.5))
    fig.subplots_adjust(left=0.11, right=0.975, top=0.92, bottom=0.16, wspace=0.28, hspace=0.40)
    for i, c in enumerate(["Cropland", "Forest"]):
        for j, s in enumerate(PATHWAYS):
            ax = axs[i, j]
            for p in FORECASTS:
                q = lc[(lc["product"] == p) & lc.landcover.eq(c) & lc.ssp.eq(s)].sort_values("year")
                ax.plot(q.year, q.contact_km / 1e6, color=COLOR[p], marker="o", ms=3, lw=1.8)
            ax.set(xticks=YEARS, ylim=(0, None))
            ax.grid(axis="y", color="#DDDDDD", lw=0.6)
            ax.set_axisbelow(True)
            ax.text(-0.075, 1.035, "ABCDEF"[i * 3 + j], transform=ax.transAxes, ha="right", va="bottom",
                    fontsize=16, fontweight="bold")
            ax.set_title(f"{c} · SSP{s}", pad=13)
            if j == 0:
                ax.set_ylabel("Contact length (million km)")
    fig.legend([Line2D([], [], color=COLOR[p], lw=2) for p in FORECASTS], LABELS, loc="lower center",
               bbox_to_anchor=(0.5, 0.012), ncol=3, frameon=False, fontsize=12)

    for ext in ("png", "pdf"):
        fig.savefig(out_path(a.out, ext), dpi=300)
    plt.close(fig)


if __name__ == "__main__":
    main()
