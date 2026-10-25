"""Divergent urban boundaries and high spatial uncertainty in global urban land forecasts

Author: Rakibul Ahasan
"""
import argparse
import json
import os
import string

import matplotlib as mpl
mpl.use("Agg")
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import patheffects as pe
from matplotlib.colors import BoundaryNorm
from shapely.geometry import shape
from shapely.ops import unary_union

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
FIG_DATA = os.path.join(REPO, "data", "figures")

mpl.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
    "font.size": 12, "axes.labelsize": 12,
    "axes.titlesize": 14, "xtick.labelsize": 11, "ytick.labelsize": 12, "axes.spines.top": False,
    "axes.spines.right": False, "pdf.fonttype": 42, "savefig.facecolor": "white",
})

BLOCK = 50000
LEVELS = [1, 5, 25, 100, 500, 1000, 2500]
PANELS = [("unanimous", "YlOrRd", "All-nine urban agreement in 2100"),
          ("inherited_disagreement", "Blues", "Disagreement in 2100 carried from 2020"),
          ("new_disagreement", "Oranges", "Disagreement formed between 2020 and 2100")]
PC = ccrs.PlateCarree()


class CEA6933(ccrs._RectangularProjection):

    def __init__(self):
        super().__init__([("proj", "cea"), ("lat_ts", 30.0), ("lon_0", 0.0)],
                         6378137 * np.pi * np.cos(np.pi / 6), 7290708.89,
                         globe=ccrs.Globe(semimajor_axis=6378137, semiminor_axis=6356752.314245179, ellipse=None))

    @property
    def threshold(self):
        return 1e5


def out_name(folder, stem, ext):
    path = os.path.join(folder, f"{stem}.{ext}")
    k = 0
    while os.path.exists(path):
        k += 1
        path = os.path.join(folder, f"{stem}_r{k}.{ext}")
    return path


def main():
    ap = argparse.ArgumentParser(description="Supplementary Fig. 7, agreement and disagreement in 2100")
    ap.add_argument("--blocks", default=os.path.join(FIG_DATA, "agreement_blocks", "blocks_50km_SSP1_2100.npz"))
    ap.add_argument("--grid", default=os.path.join(FIG_DATA, "agreement_blocks", "grid.json"))
    ap.add_argument("--boundaries", default=os.path.join(FIG_DATA, "regions", "world_boundaries.geojson"))
    ap.add_argument("--out", default=HERE)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    cea = CEA6933()
    geoms = [shape(f["geometry"]).simplify(0.04, preserve_topology=True)
             for f in json.load(open(a.boundaries, encoding="utf-8"))["features"]]
    land = unary_union(geoms)
    tr = json.load(open(a.grid))["transform"]
    blocks = np.load(a.blocks)

    def base(ax):
        ax.set_extent([-179.5, 179.5, -56, 76], PC)
        ax.spines["geo"].set_visible(False)
        coast = ax.add_geometries([land], PC, facecolor="#ECEBE7", edgecolor="#BCBBB7", linewidth=0.3, zorder=1.5)
        coast.set_path_effects([pe.SimplePatchShadow(offset=(8, -8), alpha=0.045, rho=0.15),
                                pe.SimplePatchShadow(offset=(7, -7), alpha=0.065, rho=0.15),
                                pe.SimplePatchShadow(offset=(6, -6), alpha=0.09, rho=0.15),
                                pe.SimplePatchShadow(offset=(5, -5), alpha=0.12, rho=0.15),
                                pe.SimplePatchShadow(offset=(4, -4), alpha=0.15, rho=0.15), pe.Normal()])
        ax.add_geometries(geoms, PC, facecolor="none", edgecolor="#C8C7C2", linewidth=0.18, zorder=2)
        ax.gridlines(draw_labels=False, xlocs=range(-180, 181, 60), ylocs=range(-60, 91, 30),
                     color="#CBCBC7", linewidth=0.35, zorder=0.5)

    def raster(ax, values, cmap_name):
        cmap = mpl.colormaps[cmap_name].resampled(len(LEVELS) - 1)
        cmap.set_bad((0, 0, 0, 0))
        norm = BoundaryNorm(LEVELS, cmap.N)
        v = np.asarray(values, float).copy()
        v[v < LEVELS[0]] = np.nan
        ax.imshow(v, origin="upper",
                  extent=[tr[2], tr[2] + v.shape[1] * BLOCK, tr[5] - v.shape[0] * BLOCK, tr[5]],
                  transform=cea, cmap=cmap, norm=norm, interpolation="nearest", regrid_shape=1300, zorder=3)
        return cmap, norm

    fig = plt.figure(figsize=(10, 12))
    gs = fig.add_gridspec(3, 1, left=0.035, right=0.98, top=0.945, bottom=0.03, hspace=0.12)
    for i, (key, cmap_name, title) in enumerate(PANELS):
        ax = fig.add_subplot(gs[i], projection=ccrs.EqualEarth())
        base(ax)
        cmap, norm = raster(ax, blocks[key], cmap_name)
        ax.text(0, 1.06, string.ascii_uppercase[i] + "  " + title, transform=ax.transAxes, weight="bold", fontsize=14)
        ca = ax.inset_axes([0.017, 0.17, 0.235, 0.035])
        cb = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), cax=ca, orientation="horizontal")
        cb.set_ticks([1, 25, 500, 2500])
        cb.set_ticklabels(["1", "25", "500", "2,500"])
        cb.ax.tick_params(labelsize=10, length=2)
        ax.text(0.017, 0.23, "km² per 50 km block", transform=ax.transAxes, fontsize=11)

    stem = os.path.basename(HERE)
    for ext in ("png", "pdf"):
        path = out_name(a.out, stem, ext)
        fig.savefig(path, dpi=300, bbox_inches="tight", pad_inches=0.04)
        print("wrote", path)
    plt.close(fig)


if __name__ == "__main__":
    main()
