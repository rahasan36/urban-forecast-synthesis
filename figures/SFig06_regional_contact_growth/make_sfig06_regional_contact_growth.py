"""Supplementary Fig. 6. Regional change in cropland and forest contact and its spread among forecasts.

Change in urban boundary contact from 2020 to 2100 in the 19 reporting regions under SSP1, SSP3 and SSP5,
with regions ordered by the SSP1 median of panel A.
A, B  ratio of 2100 to 2020 contact length with cropland (A) and forest (B); dots mark the median and bars
      the range across the nine forecasts.
C     ratio of the longest to the shortest forecast cropland contact length in 2100.
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.ticker import FixedFormatter, FixedLocator, NullLocator

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
FIG_DATA = os.path.join(REPO, "data", "figures")
NAME = "SFig06_regional_contact_growth"

COLOR = {1: "#0072B2", 3: "#E69F00", 5: "#CC79A7"}
OFFSET = {1: 0.25, 3: 0.0, 5: -0.25}

plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
                     "font.size": 8, "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
                     "text.color": "black", "axes.labelcolor": "black", "xtick.color": "black",
                     "ytick.color": "black", "pdf.fonttype": 42})


def out_path(folder, ext):
    path = os.path.join(folder, f"{NAME}.{ext}")
    k = 1
    while os.path.exists(path):
        path = os.path.join(folder, f"{NAME}_r{k}.{ext}")
        k += 1
    return path


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--landcover_changes", default=os.path.join(FIG_DATA, "results", "landcover_endpoint_changes.csv"))
    ap.add_argument("--out", default=HERE)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    d = pd.read_csv(a.landcover_changes)
    d = d[(d.year == 2100) & (d.threshold.isna() | (d.threshold == 0.4)) & (d.region != "Global")
          & d.landcover.isin(["Cropland", "Forest"])].copy()
    d["ratio"] = np.where(d.contact_km_2020 > 0, d.contact_km / d.contact_km_2020, np.nan)
    agg = d.groupby(["landcover", "region", "ssp"]).agg(lo=("ratio", "min"), med=("ratio", "median"),
                                                        hi=("ratio", "max"), lmax=("contact_km", "max"),
                                                        lmin=("contact_km", "min")).reset_index()
    agg["spread"] = agg.lmax / agg.lmin
    order = agg[(agg.landcover == "Cropland") & (agg.ssp == 1)].sort_values("med").region.tolist()
    ypos = {r: i for i, r in enumerate(order)}

    fig, axes = plt.subplots(1, 3, figsize=(7.2, 6.6), sharey=True,
                             gridspec_kw={"width_ratios": [1, 1, 0.75], "wspace": 0.08})
    for ax, cover, letter in [(axes[0], "Cropland", "A"), (axes[1], "Forest", "B")]:
        ax.axvline(1, color="#9A9A9A", lw=0.6, ls=(0, (3, 2)), zorder=1)
        for _, r in agg[agg.landcover == cover].iterrows():
            y = ypos[r.region] + OFFSET[r.ssp]
            ax.plot([r.lo, r.hi], [y, y], color=COLOR[r.ssp], lw=1.1, solid_capstyle="butt", zorder=2)
            ax.plot(r.med, y, "o", ms=3.6, color=COLOR[r.ssp], mec="white", mew=0.4, zorder=3)
        ax.set_xscale("log")
        ax.set_xlabel(f"{cover} contact, 2100 / 2020")
        ax.text(-0.02, 1.01, letter, transform=ax.transAxes, fontsize=14, fontweight="bold", ha="right", va="bottom")

    ax = axes[2]
    for _, r in agg[agg.landcover == "Cropland"].iterrows():
        ax.plot(r.spread, ypos[r.region] + OFFSET[r.ssp], "o", ms=3.6, color=COLOR[r.ssp], mec="white", mew=0.4)
    ax.set_xscale("log")
    ax.set_xlabel("Cropland contact in 2100,\nlongest / shortest")
    ax.text(-0.02, 1.01, "C", transform=ax.transAxes, fontsize=14, fontweight="bold", ha="right", va="bottom")

    for ax, ticks, (lo, hi) in [(axes[0], [0.2, 0.5, 1, 2, 5, 10, 20], (0.12, 30)),
                                (axes[1], [0.2, 0.5, 1, 2, 5, 10, 20, 50], (0.12, 80)),
                                (axes[2], [5, 10, 20, 50], (2.5, 65))]:
        ax.set_xlim(lo, hi)
        t = [v for v in ticks if lo <= v <= hi]
        ax.xaxis.set_major_locator(FixedLocator(t))
        ax.xaxis.set_major_formatter(FixedFormatter([f"{v:g}" for v in t]))
        ax.xaxis.set_minor_locator(NullLocator())
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
    axes[0].set_yticks(range(len(order)))
    axes[0].set_yticklabels(order)
    axes[0].set_ylim(-0.6, len(order) - 0.4)
    for ax in axes[1:]:
        ax.tick_params(axis="y", length=0)
    handles = [Line2D([0], [0], color=COLOR[k], lw=1.1, marker="o", ms=3.6, mec="white", mew=0.4, label=f"SSP{k}")
               for k in (1, 3, 5)]
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False, bbox_to_anchor=(0.55, 0.0), fontsize=8,
               handlelength=1.6)
    fig.subplots_adjust(left=0.22, right=0.985, top=0.96, bottom=0.12)

    for ext in ("png", "pdf"):
        fig.savefig(out_path(a.out, ext), dpi=300, facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    main()
