"""Divergent urban boundaries and high spatial uncertainty in global urban land forecasts

Author: Rakibul Ahasan
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
FIG_DATA = os.path.join(REPO, "data", "figures")
NAME = "SFig01_landcover_conversion"

FORECASTS = ["Li2021", "Chen2020", "Chen2022", "Gao2021", "Hou2022", "Koomen2UP", "He2023", "Zhang2023", "Zhao2025"]
LABELS = ["Li 2021", "Chen 2020", "Chen 2022", "Gao 2021", "Hou 2022", "2UP 2023", "He 2023", "Zhang 2023", "FU3D 2025"]
COLORS = ["#0077BB", "#EE7733", "#DDAA00", "#009988", "#CC3311", "#33BBEE", "#AA3377", "#777777", "#222222"]
COLOR = dict(zip(FORECASTS, COLORS))
CLASS_COLOR = {"Cropland": "#B87828", "Forest": "#228866", "Other nonurban": "#ADB7C0"}
PATHWAYS = [1, 3, 5]
YEARS = [2020, 2050, 2100]

plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
    "font.size": 14, "axes.labelsize": 14, "axes.titlesize": 16, "xtick.labelsize": 13, "ytick.labelsize": 13,
    "legend.fontsize": 13, "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": 0.8,
    "lines.linewidth": 2, "pdf.fonttype": 42, "ps.fonttype": 42, "savefig.facecolor": "white",
})


def out_path(folder, ext):
    path = os.path.join(folder, f"{NAME}.{ext}")
    k = 1
    while os.path.exists(path):
        path = os.path.join(folder, f"{NAME}_r{k}.{ext}")
        k += 1
    return path


def panel(ax, letter, title=None):
    ax.text(-0.055, 1.035, letter, transform=ax.transAxes, fontsize=19, fontweight="bold", ha="right", va="bottom")
    if title:
        ax.set_title(title, pad=14)


def ygrid(ax, axis="y"):
    ax.grid(axis=axis, color="#DDDDDD", lw=0.7)
    ax.set_axisbelow(True)


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--global_metrics", default=os.path.join(FIG_DATA, "results", "global_metrics.csv"))
    ap.add_argument("--landcover_changes", default=os.path.join(FIG_DATA, "results", "landcover_endpoint_changes.csv"))
    ap.add_argument("--out", default=HERE)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    g = pd.read_csv(a.global_metrics)
    g = g[g.threshold == 0.4]
    cl = pd.read_csv(a.landcover_changes)
    cl = cl[cl.threshold == 0.4]

    fig, axs = plt.subplots(2, 3, figsize=(14, 10.8), gridspec_kw={"height_ratios": [1, 1.4]})
    fig.subplots_adjust(left=0.11, right=0.975, bottom=0.19, top=0.935, wspace=0.27, hspace=0.44)
    base = g[g.year == 2020].set_index(["product", "ssp"]).urban_area_km2
    lim = max((r.urban_area_km2 - base.loc[(r.product, r.ssp)]) / 1e6 for r in g.itertuples()) * 1.07
    for j, s in enumerate(PATHWAYS):
        ax = axs[0, j]
        for i, p in enumerate(FORECASTS):
            q = g[(g["product"] == p) & (g.ssp == s)].sort_values("year")
            ax.plot(q.year, (q.urban_area_km2 - base.loc[p, s]) / 1e6, color=COLOR[p], marker="o", ms=4,
                    ls="--" if i >= 6 else "-")
        ax.set(xticks=YEARS, ylim=(0, [3, 4, lim][j]))
        ygrid(ax)
        panel(ax, "ABC"[j], f"SSP{s}")
        if j == 0:
            ax.set_ylabel("Additional urban area\nsince 2020 (million km²)")

        ax = axs[1, j]
        left = np.zeros(9)
        for c in ["Cropland", "Forest", "Other nonurban"]:
            values = []
            for p in FORECASTS:
                q = cl[(cl["product"] == p) & (cl.ssp == s) & (cl.region == "Global")]
                if c == "Other nonurban":
                    q = q[~q.landcover.isin(["Cropland", "Forest", "Baseline urban"])]
                else:
                    q = q[q.landcover == c]
                values.append(q.overlap_km2_change.sum() / 1e6)
            ax.barh(range(9), values, left=left, color=CLASS_COLOR[c], height=0.65)
            left += values
        ax.set_yticks(range(9), LABELS if j == 0 else [""] * 9)
        ax.invert_yaxis()
        ax.set_xlabel("Net nonurban overlap increase\n2020–2100 (million km²)")
        ygrid(ax, "x")
        panel(ax, "DEF"[j])
        ax.set_xlim(0, lim)

    axs[1, 1].legend([Patch(color=CLASS_COLOR[c]) for c in CLASS_COLOR], list(CLASS_COLOR), loc="upper left",
                     bbox_to_anchor=(2, 1.7), bbox_transform=axs[1, 1].transData, ncol=1, frameon=False,
                     fontsize=10.5, handlelength=1.4, handletextpad=0.5, borderaxespad=0, labelspacing=0.65)
    fig.legend([Line2D([], [], color=COLOR[p], lw=2.5, ls="--" if i >= 6 else "-") for i, p in enumerate(FORECASTS)],
               LABELS, loc="lower center", bbox_to_anchor=(0.5, 0.012), ncol=5, frameon=False,
               columnspacing=1.5, handlelength=2.1)

    for ext in ("png", "pdf"):
        fig.savefig(out_path(a.out, ext), dpi=300, pad_inches=0.1)
    plt.close(fig)


if __name__ == "__main__":
    main()
