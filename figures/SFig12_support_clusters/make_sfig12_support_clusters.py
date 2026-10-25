"""Supplementary Fig. 12: local clustering of urban support among the nine forecasts in 2100 under SSP1,
under the primary and the stricter specification.

a: Getis-Ord Gi* z score of block disagreement intensity in 50-km blocks with at least 10 km2 of union
urban land (primary specification, circular neighborhood of radius four blocks). b: contested and
consensus clusters under the stricter specification (25 km2 floor, at least nine neighbors, false
discovery rate control). c, d: share of eligible blocks in contested and consensus clusters by region
under the two specifications, with one region order.
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import cartopy.crs as ccrs
from cartopy.io import shapereader
from matplotlib.cm import ScalarMappable
from matplotlib.colors import BoundaryNorm, ListedColormap, TwoSlopeNorm
from matplotlib.patches import Patch

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
FIG_DATA = os.path.join(REPO, "data", "figures")
NAME = "SFig12_support_clusters"

plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Liberation Sans", "Nimbus Sans", "DejaVu Sans"],
    "font.size": 8, "axes.labelsize": 8, "xtick.labelsize": 8, "ytick.labelsize": 8, "axes.linewidth": 0.6,
    "text.color": "black", "axes.labelcolor": "black", "xtick.color": "black", "ytick.color": "black",
    "pdf.fonttype": 42, "ps.fonttype": 42,
})

BLUE, GRAY, RED = "#2B6CA3", "#C9C9C9", "#D7301F"
BLOCK = 50000.0
GRID_X0, GRID_Y0 = -17367529.34, 7289708.89   # upper-left corner of the EPSG:6933 analysis grid (m)
REGIONS = {1: "Asia-Pacific Developed", 2: "Australia & New Zealand", 3: "Caribbean", 4: "Developing Pacific",
           5: "Eastern Africa", 6: "Eastern Asia", 7: "Eurasia", 9: "India & Sri Lanka", 10: "Meso America",
           11: "Middle East", 12: "North Africa", 13: "N. & W. Europe", 14: "Rest of Southern Asia",
           15: "South America", 16: "South-East Asia", 17: "S. & E. Europe", 18: "Southern & Middle Africa",
           19: "USA & Canada", 20: "Western Africa"}


def out_name(folder, stem, ext):
    path = os.path.join(folder, f"{stem}.{ext}")
    k = 1
    while os.path.exists(path):
        path = os.path.join(folder, f"{stem}_r{k}.{ext}")
        k += 1
    return path


def basemap(fig, rect, land):
    ax = fig.add_axes(rect, projection=ccrs.Robinson())
    ax.set_extent([-168, 180, -57, 79], crs=ccrs.PlateCarree())
    ax.add_geometries(land, ccrs.PlateCarree(), facecolor="#FFFFFF", edgecolor="none", zorder=1,
                      path_effects=[pe.withSimplePatchShadow(offset=(1.2, -1.2), shadow_rgbFace="#8E9499",
                                                             alpha=0.55, rho=1.0)])
    ax.add_geometries(land, ccrs.PlateCarree(), facecolor="#F0F0F0", edgecolor="#AFAFAF", linewidth=0.2, zorder=2)
    ax.gridlines(draw_labels=False, linewidth=0.25, color="#E4E4E4", zorder=3)
    ax.spines["geo"].set_visible(False)
    return ax


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--clusters", default=os.path.join(FIG_DATA, "support_clusters"),
                    help="folder with primary_gistar.npz, strict_clusters.npz and strict_region_summary.csv")
    ap.add_argument("--land", default=os.path.join(FIG_DATA, "naturalearth", "ne_110m_land.shp"))
    ap.add_argument("--out", default=HERE)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    g = np.load(os.path.join(a.clusters, "primary_gistar.npz"))
    z, mask, crop, region = g["z"], g["eligible"], g["cropland_boundary_km"], g["region"]
    hot = (z >= 1.96) & mask
    cold = (z <= -1.96) & mask
    none = mask & ~hot & ~cold
    total = crop[mask].sum()
    print("primary: eligible, contested, consensus blocks", int(mask.sum()), int(hot.sum()), int(cold.sum()))
    print("primary: cropland boundary share (%) in consensus, no cluster, contested",
          [round(float(100 * crop[m].sum() / total), 1) for m in (cold, none, hot)])

    s = np.load(os.path.join(a.clusters, "strict_clusters.npz"))
    category, smask, scrop = s["category"], s["eligible"], s["cropland_boundary_km"]
    stotal = scrop[smask].sum()
    print("strict: cropland boundary share (%) in consensus, no cluster, contested",
          [round(float(100 * scrop[smask & (category == c)].sum() / stotal), 1) for c in (1, 0, -1)])
    summary = pd.read_csv(os.path.join(a.clusters, "strict_region_summary.csv")).set_index("region_id")

    rows = []
    for rid, name in REGIONS.items():
        sel = mask & (region == rid)
        if sel.sum() < 15 or rid not in summary.index:
            continue
        rows.append((rid, name, 100 * (hot & sel).sum() / sel.sum(), 100 * (cold & sel).sum() / sel.sum()))
    rows.sort(key=lambda r: r[2])

    land = [geom for geom in shapereader.Reader(a.land).geometries() if geom.centroid.y > -60]
    fig = plt.figure(figsize=(7.2, 5.0))
    xx = GRID_X0 + np.arange(z.shape[1] + 1) * BLOCK
    yy = GRID_Y0 - np.arange(z.shape[0] + 1) * BLOCK

    axa = basemap(fig, [0.005, 0.655, 0.490, 0.335], land)
    norm = TwoSlopeNorm(vmin=-5, vcenter=0, vmax=5)
    cmap = plt.get_cmap("RdYlBu_r")
    axa.pcolormesh(xx, yy, np.ma.array(np.nan_to_num(z), mask=~mask), transform=ccrs.epsg(6933),
                   cmap=cmap, norm=norm, shading="flat", rasterized=True, antialiased=False, zorder=6)
    axb = basemap(fig, [0.505, 0.655, 0.490, 0.335], land)
    axb.pcolormesh(xx, yy, np.ma.array(category, mask=~smask), transform=ccrs.epsg(6933),
                   cmap=ListedColormap([RED, GRAY, BLUE]), norm=BoundaryNorm([-1.5, -0.5, 0.5, 1.5], 3),
                   shading="flat", rasterized=True, antialiased=False, zorder=6)

    cax = fig.add_axes([0.150, 0.598, 0.200, 0.017])
    cb = fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), cax=cax, orientation="horizontal",
                      ticks=[-5, -1.96, 0, 1.96, 5])
    cb.ax.set_xticklabels(["−5", "−2", "0", "2", "5"])
    cb.ax.tick_params(length=2, pad=1.5)
    cb.outline.set_linewidth(0.4)
    fig.text(0.140, 0.6065, "consensus", color=BLUE, ha="right", va="center")
    fig.text(0.360, 0.6065, "contested", color=RED, ha="left", va="center")
    fig.text(0.250, 0.638, "Gi* z score", ha="center", va="center")
    fig.legend(handles=[Patch(fc=BLUE, label="consensus"), Patch(fc=GRAY, label="no cluster"),
                        Patch(fc=RED, label="contested")],
               loc="center", bbox_to_anchor=(0.750, 0.6065), ncol=3, frameon=False,
               handlelength=1.0, handletextpad=0.4, columnspacing=1.2)

    y = np.arange(len(rows))
    names = [r[1] for r in rows]
    primary_hot = np.array([r[2] for r in rows])
    primary_cold = np.array([r[3] for r in rows])
    strict_hot = np.array([summary.loc[r[0], "lower_percent"] for r in rows])
    strict_cold = np.array([summary.loc[r[0], "higher_percent"] for r in rows])
    axc = fig.add_axes([0.235, 0.085, 0.300, 0.450])
    axd = fig.add_axes([0.660, 0.085, 0.300, 0.450], sharey=axc)
    for ax, contested, consensus in ((axc, primary_hot, primary_cold), (axd, strict_hot, strict_cold)):
        ax.barh(y, contested, color=RED, height=0.72, lw=0, zorder=3)
        ax.barh(y, -consensus, color=BLUE, height=0.72, lw=0, zorder=3)
        ax.axvline(0, color="black", lw=0.6, zorder=4)
        ax.set_xlim(-70, 70)
        ax.set_xticks([-60, -30, 0, 30, 60])
        ax.set_xticklabels(["60", "30", "0", "30", "60"])
        ax.set_xlabel("Share of eligible blocks (%)", labelpad=2)
        for side in ("top", "right", "left"):
            ax.spines[side].set_visible(False)
        ax.grid(axis="x", color="#E6E6E6", lw=0.5)
        ax.set_axisbelow(True)
        ax.tick_params(axis="y", length=0, pad=2)
        ax.tick_params(axis="x", length=2)
        ax.text(-66, len(rows) - 0.4, "consensus", color=BLUE, ha="left", va="bottom")
        ax.text(66, len(rows) - 0.4, "contested", color=RED, ha="right", va="bottom")
    axc.set_yticks(y)
    axc.set_yticklabels(names)
    axc.set_ylim(-0.6, len(rows) + 0.6)
    plt.setp(axd.get_yticklabels(), visible=False)
    for letter, (x0, y0) in zip("abcd", [(0.006, 0.985), (0.506, 0.985), (0.006, 0.565), (0.560, 0.565)]):
        fig.text(x0, y0, letter, fontsize=10, fontweight="bold", va="top")

    for ext in ("png", "pdf"):
        fig.savefig(out_name(a.out, NAME, ext), dpi=300, facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    main()
