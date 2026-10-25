"""Supplementary Fig. 8: exact forecast intersections in the nine-forecast ensemble, SSP1.

Each 1-km cell carries a membership code, the set of forecasts that call it urban (bit i
set for forecast i, in the order of the matrix rows). (A to C) Area of the 20 largest
nonempty combinations in 2100, in 2020, 2050 and 2100, with the all-nine combination in
red; the annotation gives the share of that year's union the combinations cover.
(D) Membership matrix of the combinations.
"""
import argparse
import math
import os

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
FIG_DATA = os.path.join(REPO, "data", "figures")

mpl.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
    "font.size": 12, "axes.spines.top": False, "axes.spines.right": False, "pdf.fonttype": 42,
})

YEARS = [2020, 2050, 2100]
FORECASTS = ["Li 2021", "Chen 2020", "Chen 2022", "Gao 2021", "Hou 2022", "2UP 2023", "He 2023",
             "Zhang 2023", "FU3D 2025"]
ALL_NINE = 511


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
    ap = argparse.ArgumentParser(description="Supplementary Fig. 8, exact forecast intersections")
    ap.add_argument("--intersections", default=os.path.join(FIG_DATA, "intersections", "intersections_SSP1.csv"),
                    help="area (km2) of every membership code 1 to 511 by year")
    ap.add_argument("--out", default=HERE)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    d = pd.read_csv(a.intersections)
    hist = {}
    for y in YEARS:
        h = np.zeros(512, np.int64)
        q = d[d.year == y]
        h[q.membership_code.to_numpy()] = q.area_km2.to_numpy()
        hist[y] = h
    codes = sorted(range(1, 512), key=lambda c: (-hist[2100][c], c))[:20]
    if ALL_NINE not in codes:
        codes.append(ALL_NINE)

    fig = plt.figure(figsize=(11.5, 11))
    gs = fig.add_gridspec(4, 1, left=0.16, right=0.965, top=0.945, bottom=0.06, hspace=0.35,
                          height_ratios=[1, 1, 1, 1.35])
    for j, y in enumerate(YEARS):
        ax = fig.add_subplot(gs[j])
        v = np.array([hist[y][c] for c in codes]) / 1000
        upper = math.ceil(float(v.max()) / 100) * 100
        ax.bar(range(len(codes)), v, color=["#A50026" if c == ALL_NINE else "#597C92" for c in codes], width=0.76)
        ax.set(xlim=(-0.7, len(codes) - 0.3), ylim=(0, upper), xticks=[], ylabel="Area\n(thousand km²)")
        ax.set_yticks(np.linspace(0, upper, 5))
        ax.set_axisbelow(True)
        ax.grid(axis="y", color="#DDDDDD", lw=0.6)
        label(ax, chr(65 + j), str(y))
        coverage = sum(hist[y][c] for c in codes) / sum(hist[y][1:]) * 100
        ax.text(0.985, 0.90, f"{coverage:.1f}% of union shown", transform=ax.transAxes, ha="right", fontsize=11)

    ax = fig.add_subplot(gs[3])
    n = len(FORECASTS)
    ax.scatter(np.tile(range(len(codes)), n), np.repeat(range(n), len(codes)), s=17, color="#DDDDDD")
    for k, c in enumerate(codes):
        on = [i for i in range(n) if c & (1 << i)]
        color = "#A50026" if c == ALL_NINE else "#222222"
        ax.plot([k, k], [min(on), max(on)], color=color, lw=1.2)
        ax.scatter([k] * len(on), on, s=29, color=color)
    ax.set(yticks=range(n), yticklabels=FORECASTS, xlim=(-0.7, len(codes) - 0.3), xticks=[])
    ax.invert_yaxis()
    ax.tick_params(axis="y", length=0)
    for sp in ax.spines.values():
        sp.set_visible(False)
    label(ax, "D")

    stem = os.path.basename(HERE)
    for ext in ("png", "pdf"):
        path = out_name(a.out, stem, ext)
        fig.savefig(path, dpi=300)
        print("wrote", path)
    plt.close(fig)


if __name__ == "__main__":
    main()
