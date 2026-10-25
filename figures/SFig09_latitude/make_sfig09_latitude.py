"""Supplementary Fig. 9: latitudinal distribution of unanimous urban agreement and forecast variation.

(A to C) Area all nine forecasts call urban in each 2-degree latitude band, under SSP1, SSP3
and SSP5. (D to F) Urban area of each forecast as a share of the band's union: the line is
the median across the nine forecasts and the shading the range. Bands with a union below
2,000 km2 are left out.
"""
import argparse
import os

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
FIG_DATA = os.path.join(REPO, "data", "figures")

mpl.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
    "font.size": 12, "axes.spines.top": False, "axes.spines.right": False, "pdf.fonttype": 42,
})

YEARS = [2020, 2050, 2100]
YEAR_COLORS = {2020: "#4477AA", 2050: "#99806A", 2100: "#CC5555"}
MIN_UNION = 2000


def out_name(folder, stem, ext):
    path = os.path.join(folder, f"{stem}.{ext}")
    k = 0
    while os.path.exists(path):
        k += 1
        path = os.path.join(folder, f"{stem}_r{k}.{ext}")
    return path


def label(ax, letter, title=""):
    ax.text(-0.075, 1.045, letter, transform=ax.transAxes, ha="right", weight="bold", fontsize=17)
    if title:
        ax.set_title(title, pad=12)


def main():
    ap = argparse.ArgumentParser(description="Supplementary Fig. 9, latitude profiles of agreement")
    ap.add_argument("--agreement", default=os.path.join(FIG_DATA, "latitude", "latitude_agreement.csv"),
                    help="union and all-nine area by 2-degree band")
    ap.add_argument("--shares", default=os.path.join(FIG_DATA, "latitude", "latitude_product_shares.csv"),
                    help="urban area of each forecast by 2-degree band")
    ap.add_argument("--out", default=HERE)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    lat = pd.read_csv(a.agreement)
    shares = pd.read_csv(a.shares)

    fig, axs = plt.subplots(2, 3, figsize=(12, 10.4), sharey=True)
    fig.subplots_adjust(left=0.09, right=0.975, top=0.925, bottom=0.15, wspace=0.29, hspace=0.35)
    for j, s in enumerate([1, 3, 5]):
        for y in YEARS:
            q = lat[(lat.ssp == s) & (lat.year == y)].sort_values("latitude_mid")
            axs[0, j].step(q.unanimous_km2 / 1000, q.latitude_mid, where="mid", color=YEAR_COLORS[y], lw=1.7)
            z = (shares[(shares.ssp == s) & (shares.year == y)].groupby("latitude_mid")
                 .agg(low=("share_percent", "min"), high=("share_percent", "max"),
                      median=("share_percent", "median"), union=("union_km2", "first")).reset_index())
            z.loc[z.union < MIN_UNION, ["low", "high", "median"]] = np.nan
            axs[1, j].fill_betweenx(z.latitude_mid, z.low, z.high, color=YEAR_COLORS[y], alpha=0.16, linewidth=0)
            axs[1, j].plot(z["median"], z.latitude_mid, color=YEAR_COLORS[y], lw=2.8)
        for i in (0, 1):
            ax = axs[i, j]
            ax.set_ylim(-58, 78)
            ax.grid(color="#DDDDDD", lw=0.5)
            ax.set_axisbelow(True)
            label(ax, chr(65 + i * 3 + j), f"SSP{s}")
            if j == 0:
                ax.set_ylabel("Latitude (°N)")
        axs[0, j].set_xlabel("All-nine area\n(thousand km² per 2° band)")
        axs[1, j].set(xlim=(0, 100), xticks=[0, 25, 50, 75, 100],
                      xlabel="Product urban area as a share\nof band union (%)")
    axs[0, 0].legend([Line2D([], [], color=YEAR_COLORS[y], lw=2.8) for y in YEARS], [str(y) for y in YEARS],
                     loc="upper left", bbox_to_anchor=(15, -10), bbox_transform=axs[0, 0].transData, ncol=1,
                     frameon=False, borderaxespad=0)

    stem = os.path.basename(HERE)
    for ext in ("png", "pdf"):
        path = out_name(a.out, stem, ext)
        fig.savefig(path, dpi=300)
        print("wrote", path)
    plt.close(fig)


if __name__ == "__main__":
    main()
