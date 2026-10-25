"""Figure 1. Variation among the nine forecasts in urban boundary length relative to urban area.

a-c  global urban boundary length of each forecast in 2020, 2050 and 2100 under SSP1, SSP3 and SSP5.
d    boundary length per unit urban area under SSP1.
e    2100/2020 ratio of boundary length against the 2100/2020 ratio of urban area, all forecasts and pathways.
f-h  coefficient of variation of urban area and of boundary length across the nine forecasts, with the
     range over the nine leave-one-out ensembles.
"""
import argparse
import os
import statistics

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import pandas as pd
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
FIG_DATA = os.path.join(REPO, "data", "figures")
NAME = "Fig1_boundary_dispersion"

FORECASTS = ["Li2021", "Chen2020", "Chen2022", "Gao2021", "Hou2022", "Koomen2UP", "He2023", "Zhang2023", "Zhao2025"]
LABEL = {"Li2021": "Li 2021", "Chen2020": "Chen 2020", "Chen2022": "Chen 2022", "Gao2021": "Gao 2021",
         "Hou2022": "Hou 2022", "Koomen2UP": "2UP 2023", "He2023": "He 2023", "Zhang2023": "Zhang 2023",
         "Zhao2025": "FU3D 2025"}
COLOR = {"Li2021": "#0077BB", "Chen2020": "#EE7733", "Chen2022": "#DDAA00", "Gao2021": "#009988",
         "Hou2022": "#CC3311", "Koomen2UP": "#33BBEE", "He2023": "#AA3377", "Zhang2023": "#777777",
         "Zhao2025": "#222222"}
DASH = {"He2023": (0, (4, 1.6)), "Zhang2023": (0, (4, 1.6)), "Zhao2025": (0, (4, 1.6))}
YEARS = (2020, 2050, 2100)
PATHWAYS = (1, 3, 5)

plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
    "font.size": 7.4, "axes.labelsize": 7.6, "xtick.labelsize": 7.0, "ytick.labelsize": 7.0,
    "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": 0.7,
    "pdf.fonttype": 42, "savefig.facecolor": "white",
})


def out_path(folder, ext):
    path = os.path.join(folder, f"{NAME}.{ext}")
    k = 1
    while os.path.exists(path):
        path = os.path.join(folder, f"{NAME}_r{k}.{ext}")
        k += 1
    return path


def cv(values):
    return statistics.pstdev(values) / statistics.mean(values)


def panel_letter(fig, x, y, letter):
    fig.text(x, y, letter, fontsize=9.5, fontweight="bold", va="top")


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--table", default=os.path.join(FIG_DATA, "results", "global_boundary_area.csv"),
                    help="global urban area, boundary length and boundary length per unit area by forecast, pathway and year")
    ap.add_argument("--out", default=HERE)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    g = pd.read_csv(a.table)

    fig = plt.figure(figsize=(7.20, 7.30))

    # a-c boundary length trajectories
    for j, s in enumerate(PATHWAYS):
        ax = fig.add_axes([0.105 + j * 0.305, 0.800, 0.258, 0.185])
        for p in FORECASTS:
            q = g[(g["product"] == p) & (g.ssp == s)].set_index("year").reindex(YEARS)
            ax.plot(YEARS, q.perimeter_km / 1e6, color=COLOR[p], lw=1.15, ls=DASH.get(p, "-"),
                    marker="o", ms=2.2, mew=0, zorder=3)
        ax.set_xlim(2012, 2108)
        ax.set_ylim(0, 5.3)
        ax.set_xticks(YEARS)
        ax.grid(axis="y", color="#ECECEC", lw=0.5)
        ax.set_axisbelow(True)
        ax.tick_params(length=2)
        if j == 0:
            ax.set_ylabel("urban boundary length\n(million km)", labelpad=2)
        else:
            ax.set_yticklabels([])
        ax.text(0.035, 0.955, f"SSP{s}", transform=ax.transAxes, fontsize=7.4, va="top")
        panel_letter(fig, 0.105 + j * 0.305 - 0.052, 0.998, "abc"[j])

    handles = [Line2D([], [], color=COLOR[p], lw=1.3, ls=DASH.get(p, "-"), label=LABEL[p]) for p in FORECASTS]
    fig.legend(handles=handles, loc="center", bbox_to_anchor=(0.530, 0.735), ncol=5, frameon=False,
               fontsize=6.9, handlelength=1.9, columnspacing=1.15, handletextpad=0.5, labelspacing=0.35)

    # d boundary length per unit urban area, SSP1
    axd = fig.add_axes([0.118, 0.400, 0.300, 0.270])
    ed = g[g.ssp == 1].pivot_table(index="product", columns="year", values="ed")
    nudge = {"Gao2021": 0.055, "Zhao2025": -0.055, "Li2021": 0.0}
    for p in FORECASTS:
        axd.plot([0, 1, 2], [ed.loc[p, y] for y in YEARS], color=COLOR[p], lw=1.15, ls=DASH.get(p, "-"),
                 marker="o", ms=2.6, mew=0, zorder=3)
        if p in nudge:
            axd.text(2.07, ed.loc[p, 2100] + nudge[p], LABEL[p], fontsize=6.6, color=COLOR[p], va="center")
    axd.set_xlim(-0.08, 2.62)
    axd.set_xticks([0, 1, 2])
    axd.set_xticklabels([str(y) for y in YEARS])
    axd.set_ylim(0.20, 2.15)
    axd.set_ylabel("boundary length per unit urban\narea (km per km$^2$), SSP1", labelpad=2)
    axd.grid(axis="y", color="#ECECEC", lw=0.5)
    axd.set_axisbelow(True)
    axd.tick_params(length=2)
    panel_letter(fig, 0.020, 0.690, "d")

    # e boundary length ratio against urban area ratio
    axe = fig.add_axes([0.596, 0.400, 0.376, 0.270])
    marker = {1: "o", 3: "s", 5: "^"}
    for s in PATHWAYS:
        for p in FORECASTS:
            q = g[(g["product"] == p) & (g.ssp == s)].set_index("year")
            area_ratio = q.loc[2100, "urban_area_km2"] / q.loc[2020, "urban_area_km2"]
            boundary_ratio = q.loc[2100, "perimeter_km"] / q.loc[2020, "perimeter_km"]
            axe.plot(area_ratio, boundary_ratio, marker=marker[s], ms=4.2, mfc=COLOR[p], mec="white", mew=0.5, zorder=4)
    lo, hi = 1.02, 5.25
    axe.plot([lo, hi], [lo, hi], color="#999999", lw=0.8, ls=(0, (3, 2)), zorder=2)
    axe.text(3.55, 3.95, "1:1", fontsize=6.6, color="#777777", rotation=40, ha="center")
    axe.set_xscale("log")
    axe.set_yscale("log")
    axe.set_xlim(lo, hi)
    axe.set_ylim(lo, hi)
    for axis in (axe.xaxis, axe.yaxis):
        axis.set_major_formatter(mticker.FuncFormatter(lambda v, pos: f"{v:g}"))
        axis.set_minor_formatter(mticker.NullFormatter())
    axe.set_xticks([1, 1.5, 2, 3, 4, 5])
    axe.set_yticks([1, 1.5, 2, 3, 4, 5])
    axe.set_xlabel("urban area, 2100 / 2020 (x)", labelpad=2)
    axe.set_ylabel("urban boundary length, 2100 / 2020 (x)", labelpad=2)
    axe.grid(color="#ECECEC", lw=0.5)
    axe.set_axisbelow(True)
    axe.tick_params(length=2, which="both")
    axe.text(0.030, 0.958, "above the line boundary grows faster than area", transform=axe.transAxes,
             fontsize=6.3, color="#555555", va="top")
    axe.legend(handles=[Line2D([], [], marker=marker[s], color="none", mfc="#666666", mec="white", mew=0.5,
                               ms=4.2, label=f"SSP{s}") for s in PATHWAYS],
               loc="lower right", frameon=False, fontsize=6.6, handletextpad=0.2, labelspacing=0.3, borderpad=0.2)
    panel_letter(fig, 0.500, 0.690, "e")

    # f-h dispersion of urban area and boundary length across the nine forecasts
    width, gap, x0 = 0.262, 0.034, 0.118
    axes = []
    for i, s in enumerate(PATHWAYS):
        ax = fig.add_axes([x0 + i * (width + gap), 0.070, width, 0.230], sharey=axes[0] if axes else None)
        axes.append(ax)
        for metric, col, mk, lab in [("urban_area_km2", "#4477AA", "o", "urban area"),
                                     ("perimeter_km", "#CC6677", "s", "boundary length")]:
            base, low, high = [], [], []
            for y in YEARS:
                q = g[(g.ssp == s) & (g.year == y)].set_index("product")[metric]
                base.append(cv([q.loc[p] for p in FORECASTS]))
                loo = [cv([q.loc[p] for p in FORECASTS if p != d]) for d in FORECASTS]
                low.append(min(loo))
                high.append(max(loo))
            ax.fill_between(YEARS, low, high, color=col, alpha=0.18, lw=0)
            ax.plot(YEARS, base, color=col, marker=mk, ms=3.2, lw=1.15, mew=0, label=lab)
        ax.set_xticks(YEARS)
        ax.set_xlim(2010, 2110)
        ax.set_ylim(0, 0.80)
        ax.grid(axis="y", color="#ECECEC", lw=0.5)
        ax.set_axisbelow(True)
        ax.tick_params(length=2)
        ax.set_title(f"SSP{s}", fontsize=7.6, fontweight="bold", pad=3)
        if i == 0:
            ax.set_ylabel("coefficient of variation\nacross the nine products", labelpad=2)
            ax.legend(loc="upper left", bbox_to_anchor=(0.0, 1.0), frameon=False, fontsize=6.9,
                      handlelength=1.6, labelspacing=0.3)
        else:
            plt.setp(ax.get_yticklabels(), visible=False)
        panel_letter(fig, x0 + i * (width + gap) - (0.098 if i == 0 else 0.022), 0.318, "fgh"[i])
    axes[2].text(0.98, 0.03, "shaded band, range over the\nnine leave-one-out ensembles", transform=axes[2].transAxes,
                 ha="right", va="bottom", fontsize=6.2, color="#555555")

    for ext in ("png", "pdf"):
        fig.savefig(out_path(a.out, ext), dpi=300)
    plt.close(fig)


if __name__ == "__main__":
    main()
