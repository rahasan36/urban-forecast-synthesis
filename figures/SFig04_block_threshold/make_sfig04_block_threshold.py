"""Supplementary Fig. 4. Block size, date, and the dispersion of boundary length and urban area.

Blocks of 100 x 100 km.
a     share of blocks where the CV of boundary length across the nine forecasts exceeds the CV of urban area,
      as a function of the minimum urban area every forecast must map in the block; SSP1 bold, SSP3 and SSP5
      thin; 2020, 2050 and 2100.
b     cumulative distributions of the block CV of urban area and of boundary length, SSP1, 2020 (solid) and
      2100 (dashed), blocks with at least 10 km2 in every forecast.
c, d  CV of urban area and of boundary length by quintile of block urban area (mean across forecasts),
      SSP1 2020 (c) and 2100 (d).
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
FIG_DATA = os.path.join(REPO, "data", "figures")
NAME = "SFig04_block_threshold"

YEAR_COLOR = {2020: "#9DBCE3", 2050: "#3B6FB6", 2100: "#1B2A55"}
AREA_COLOR, BOUNDARY_COLOR = "#3B6FB6", "#D6336C"

plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
                     "font.size": 8, "pdf.fonttype": 42, "axes.linewidth": 0.6})


def out_path(folder, ext):
    path = os.path.join(folder, f"{NAME}.{ext}")
    k = 1
    while os.path.exists(path):
        path = os.path.join(folder, f"{NAME}_r{k}.{ext}")
        k += 1
    return path


def load(blocks_dir, ssp, year):
    """Minimum and mean urban area across forecasts, and the CV of urban area and of boundary length, per block."""
    z = np.load(os.path.join(blocks_dir, f"blocks100km_{ssp}_{year}.npz"))
    area = z["area"].astype(float)
    boundary = z["bnd"].astype(float)
    with np.errstate(all="ignore"):
        return area.min(0), area.mean(0), area.std(0) / area.mean(0), boundary.std(0) / boundary.mean(0)


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--blocks", default=os.path.join(FIG_DATA, "blocks100km"),
                    help="folder with blocks100km_<SSP>_<year>.npz")
    ap.add_argument("--out", default=HERE)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    fig = plt.figure(figsize=(7.2, 5.4))

    # a share of blocks with boundary more variable against the area floor
    axa = fig.add_axes([0.085, 0.58, 0.38, 0.37])
    floors = np.unique(np.round(np.logspace(0, np.log10(3000), 40)))
    for ssp, lw in (("SSP1", 1.6), ("SSP3", 0.7), ("SSP5", 0.7)):
        for year in (2020, 2050, 2100):
            mn, _, cv_area, cv_boundary = load(a.blocks, ssp, year)
            share = []
            for f in floors:
                ok = mn >= f
                share.append(100 * (cv_boundary[ok] > cv_area[ok]).mean() if ok.sum() >= 30 else np.nan)
            axa.plot(floors, share, color=YEAR_COLOR[year], lw=lw, ls="-" if ssp != "SSP5" else (0, (3, 1.5)))
    axa.set_xscale("log")
    axa.axhline(50, color="black", lw=0.6, ls=(0, (2, 2)))
    axa.set_ylim(0, 100)
    axa.set_yticks([0, 25, 50, 75, 100])
    axa.set_xlim(1, 3000)
    axa.set_xticks([1, 10, 100, 1000])
    axa.set_xticklabels(["1", "10", "100", "1,000"])
    axa.minorticks_off()
    axa.set_xlabel("Minimum urban area per product in block (km²)")
    axa.set_ylabel("Blocks with boundary more variable (%)")
    axa.spines[["top", "right"]].set_visible(False)
    axa.tick_params(length=2)
    handles = [Line2D([], [], color=YEAR_COLOR[y], lw=1.6, label=str(y)) for y in (2020, 2050, 2100)]
    handles += [Line2D([], [], color="#777777", lw=1.6, label="SSP1"),
                Line2D([], [], color="#777777", lw=0.7, label="SSP3"),
                Line2D([], [], color="#777777", lw=0.7, ls=(0, (3, 1.5)), label="SSP5")]
    axa.legend(handles=handles, loc="lower left", ncol=2, frameon=False, handlelength=1.6, columnspacing=0.8,
               labelspacing=0.25)
    fig.text(0.006, 0.985, "a", fontsize=10, fontweight="bold", va="top")

    # b cumulative distributions of the CV
    axb = fig.add_axes([0.6, 0.58, 0.38, 0.37])
    for year, ls in ((2020, "-"), (2100, (0, (4, 2)))):
        mn, _, cv_area, cv_boundary = load(a.blocks, "SSP1", year)
        ok = mn >= 10
        for v, c in ((cv_area[ok], AREA_COLOR), (cv_boundary[ok], BOUNDARY_COLOR)):
            x = np.sort(v)
            axb.plot(x, np.arange(1, x.size + 1) / x.size, color=c, lw=1.3, ls=ls)
    axb.set_xlim(0, 1.6)
    axb.set_ylim(0, 1)
    axb.set_xlabel("Coefficient of variation among products")
    axb.set_ylabel("Cumulative share of blocks")
    axb.spines[["top", "right"]].set_visible(False)
    axb.tick_params(length=2)
    handles = [Line2D([], [], color=AREA_COLOR, lw=1.3, label="Urban area"),
               Line2D([], [], color=BOUNDARY_COLOR, lw=1.3, label="Boundary length"),
               Line2D([], [], color="#777777", lw=1.3, label="2020"),
               Line2D([], [], color="#777777", lw=1.3, ls=(0, (4, 2)), label="2100")]
    axb.legend(handles=handles, loc="lower right", frameon=False, handlelength=1.8, labelspacing=0.25,
               title="SSP1", title_fontsize=8)
    fig.text(0.515, 0.985, "b", fontsize=10, fontweight="bold", va="top")

    # c, d CV by quintile of block urban area
    for j, year in enumerate((2020, 2100)):
        mn, m, cv_area, cv_boundary = load(a.blocks, "SSP1", year)
        ok = mn >= 10
        q = np.quantile(m[ok], np.linspace(0, 1, 6))
        quintile = np.clip(np.searchsorted(q, m[ok], side="right") - 1, 0, 4)
        ax = fig.add_axes([0.085 + j * 0.515, 0.1, 0.38, 0.33])
        for k in range(5):
            for v, c, off in ((cv_area[ok], AREA_COLOR, -0.17), (cv_boundary[ok], BOUNDARY_COLOR, 0.17)):
                bp = ax.boxplot([v[quintile == k]], positions=[k + off], widths=0.3, patch_artist=True,
                                showfliers=False, whis=(5, 95), medianprops=dict(color="black", lw=0.7),
                                boxprops=dict(lw=0.5), whiskerprops=dict(lw=0.5), capprops=dict(lw=0.5))
                bp["boxes"][0].set_facecolor(c)
                bp["boxes"][0].set_edgecolor("#333333")
        ax.set_xticks(range(5))
        ax.set_xticklabels([f"≤{int(q[k + 1]):,}" for k in range(5)])
        ax.set_xlim(-0.6, 4.6)
        ax.set_ylim(0, 2.0)
        ax.set_xlabel(f"Block urban area, quintile upper bound (km²), {year}")
        ax.set_ylabel("Coefficient of variation\namong products")
        ax.spines[["top", "right"]].set_visible(False)
        ax.tick_params(length=2)
        if j == 0:
            ax.legend(handles=[plt.Rectangle((0, 0), 1, 1, fc=AREA_COLOR, ec="#333333", lw=0.5),
                               plt.Rectangle((0, 0), 1, 1, fc=BOUNDARY_COLOR, ec="#333333", lw=0.5)],
                      labels=["Urban area", "Boundary length"], loc="upper right", frameon=False, handlelength=1.0)
        fig.text(0.006 + j * 0.509, 0.475, "cd"[j], fontsize=10, fontweight="bold", va="top")

    for ext in ("png", "pdf"):
        fig.savefig(out_path(a.out, ext), dpi=300, facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    main()
