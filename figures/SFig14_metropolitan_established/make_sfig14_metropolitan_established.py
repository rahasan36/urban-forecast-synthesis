"""Supplementary Fig. 14: number of forecasts classifying each 1-km cell as urban in 2020 and 2100 under SSP1
in seven large metropolitan regions, in windows of 160 x 160 km.

Each panel gives the share of the urban union that all nine forecasts, and that only one forecast,
classify as urban.
"""
import argparse
import os

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import rasterio
from matplotlib.colors import BoundaryNorm
from matplotlib.patches import Patch

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
FIG_DATA = os.path.join(REPO, "data", "figures")
NAME = "SFig14_metropolitan_established"

FONT = dict(inset=6, title=8, year=9, legend=7, letter=9)
plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
                     "font.size": 7, "pdf.fonttype": 42})

# (title, file name stem in data/figures/metropolitan)
METROS = [
    ("Houston", "houston"),
    ("Dallas–Fort Worth", "dallas_fort_worth"),
    ("London", "london"),
    ("Paris", "paris"),
    ("Tokyo", "tokyo"),
    ("Shanghai", "shanghai"),
    ("Hong Kong", "hong_kong"),
]
WIDTH_KM = 160
LAYOUT = dict(width=7.0, left=0.30, right=0.04, gap=0.06, bottom=0.52, vgap=0.07, title_shift=0.012)
WRAP = {"New York–Philadelphia": "New York–\nPhiladelphia", "Dallas–Fort Worth": "Dallas–\nFort Worth"}
NO_DATA = 65535


def out_name(folder, stem, ext):
    path = os.path.join(folder, f"{stem}.{ext}")
    k = 1
    while os.path.exists(path):
        path = os.path.join(folder, f"{stem}_r{k}.{ext}")
        k += 1
    return path


def product_count(mask):
    """Number of forecasts classifying each cell as urban, from the 9-bit membership mask."""
    lut = np.array([bin(i).count("1") for i in range(512)])
    valid = mask != NO_DATA
    return lut[np.where(valid, mask, 0)], valid


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--windows", default=os.path.join(FIG_DATA, "metropolitan"),
                    help="folder with <city>_<width>km_SSP1_<year>.tif")
    ap.add_argument("--out", default=HERE)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    L = LAYOUT
    n = len(METROS)
    W = L["width"]
    side = (W - L["left"] - L["right"] - (n - 1) * L["gap"]) / n
    top_row = L["bottom"] + side + L["vgap"]
    H = top_row + side + 0.40
    cmap = mpl.colormaps["YlOrRd"].resampled(9)
    norm = BoundaryNorm(np.arange(0.5, 10), 9)

    fig = plt.figure(figsize=(W, H))
    for i, (title, stem) in enumerate(METROS):
        axes = []
        for j, year in enumerate([2020, 2100]):
            y0 = top_row if j == 0 else L["bottom"]
            ax = fig.add_axes([(L["left"] + i * (side + L["gap"])) / W, y0 / H, side / W, side / H])
            axes.append(ax)
            with rasterio.open(os.path.join(a.windows, f"{stem}_{WIDTH_KM}km_SSP1_{year}.tif")) as src:
                count, valid = product_count(src.read(1))
            ax.imshow(np.ma.masked_where((count == 0) | ~valid, count), cmap=cmap, norm=norm, interpolation="nearest")
            ax.set_xticks([])
            ax.set_yticks([])
            for sp in ax.spines.values():
                sp.set_color("#AAAAAA")
                sp.set_linewidth(0.5)
            union = int((count > 0).sum())
            all_nine = 100 * (count == 9).sum() / union
            single = 100 * (count == 1).sum() / union
            ax.text(0.03, 0.03, f"All nine: {all_nine:.1f}%\nOne product: {single:.1f}%", transform=ax.transAxes,
                    ha="left", va="bottom", fontsize=FONT["inset"], linespacing=1.15,
                    bbox=dict(facecolor="white", edgecolor="none", alpha=0.88, pad=1.0))
            print(f"{title}, {year}: union {union:,} km², all nine {all_nine:.1f}%, one forecast {single:.1f}%")
        p0 = axes[0].get_position()
        fig.text(p0.x0, p0.y1 + 0.012, chr(65 + i), weight="bold", fontsize=FONT["letter"], va="bottom")
        fig.text(p0.x0 + p0.width / 2 + L["title_shift"], p0.y1 + 0.012, WRAP.get(title, title), ha="center",
                 va="bottom", fontsize=FONT["title"], linespacing=1.1)
    for j, year in enumerate([2020, 2100]):
        y0 = top_row if j == 0 else L["bottom"]
        fig.text(0.15 / W, (y0 + side / 2) / H, str(year), ha="center", va="center", rotation=90,
                 fontsize=FONT["year"], weight="bold")
    fig.legend([Patch(facecolor=cmap(i)) for i in range(9)], [str(i) for i in range(1, 10)], loc="lower center",
               bbox_to_anchor=((L["left"] + (W - L["left"] - L["right"]) / 2) / W, 0.01), ncol=9, frameon=False,
               title="Number of products classifying a cell as urban", fontsize=FONT["legend"],
               title_fontsize=FONT["legend"], handlelength=1.3, columnspacing=1.6)
    for ext in ("png", "pdf"):
        fig.savefig(out_name(a.out, NAME, ext), dpi=300)
    plt.close(fig)


if __name__ == "__main__":
    main()
