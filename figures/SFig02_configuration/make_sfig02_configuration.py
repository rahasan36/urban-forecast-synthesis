"""Divergent urban boundaries and high spatial uncertainty in global urban land forecasts

Author: Rakibul Ahasan
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import pandas as pd
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
FIG_DATA = os.path.join(REPO, "data", "figures")
NAME = "SFig02_configuration"

FORECASTS = ["Li2021", "Chen2020", "Chen2022", "Gao2021", "Hou2022", "Koomen2UP", "He2023", "Zhang2023", "Zhao2025"]
LABELS = ["Li 2021", "Chen 2020", "Chen 2022", "Gao 2021", "Hou 2022", "2UP 2023", "He 2023", "Zhang 2023", "FU3D 2025"]
COLORS = ["#0077BB", "#EE7733", "#DDAA00", "#009988", "#CC3311", "#33BBEE", "#AA3377", "#777777", "#222222"]
COLOR = dict(zip(FORECASTS, COLORS))
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


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--configuration", default=os.path.join(FIG_DATA, "results", "configuration_metrics.csv"))
    ap.add_argument("--out", default=HERE)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    conf = pd.read_csv(a.configuration)

    fig, axs = plt.subplots(2, 3, figsize=(14, 9.8))
    fig.subplots_adjust(left=0.1, right=0.975, top=0.93, bottom=0.15, wspace=0.28, hspace=0.34)
    metrics = ["mean_patch_area_ha", "edge_density_m_per_ha"]
    for j, s in enumerate(PATHWAYS):
        for i, p in enumerate(FORECASTS):
            q = conf[(conf["product"] == p) & (conf.ssp == s) & (conf.region == "Global")].sort_values("year")
            for row, m in enumerate(metrics):
                axs[row, j].plot(q.year, q[m] / q[m].iloc[0], marker="o", ms=4, color=COLOR[p],
                                 ls="--" if i >= 6 else "-")
        for row in range(2):
            ax = axs[row, j]
            ax.set_xticks(YEARS)
            ax.set_yscale("log")
            ax.set_ylim((0.7, 8) if row == 0 else (0.9, 4.5))
            ax.set_yticks([1, 2, 4, 8] if row == 0 else [1, 2, 3, 4])
            ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda v, pos: f"{v:g}"))
            ax.yaxis.set_minor_formatter(mticker.NullFormatter())
            ax.axhline(1, color="#888", ls="--", lw=0.9)
            ax.grid(axis="y", color="#DDDDDD", lw=0.7)
            ax.set_axisbelow(True)
            ax.text(-0.055, 1.035, "ABCDEF"[row * 3 + j], transform=ax.transAxes, fontsize=19,
                    fontweight="bold", ha="right", va="bottom")
            if row == 0:
                ax.set_title(f"SSP{s}", pad=14)
    axs[0, 0].set_ylabel("Mean patch area / 2020 (×)")
    axs[1, 0].set_ylabel("Edge density / 2020 (×)")
    fig.legend([Line2D([], [], color=COLOR[p], lw=2.5, ls="--" if i >= 6 else "-") for i, p in enumerate(FORECASTS)],
               LABELS, loc="lower center", bbox_to_anchor=(0.5, 0.01), ncol=5, frameon=False,
               columnspacing=1.5, handlelength=2.1)

    for ext in ("png", "pdf"):
        fig.savefig(out_path(a.out, ext), dpi=300, pad_inches=0.1)
    plt.close(fig)


if __name__ == "__main__":
    main()
