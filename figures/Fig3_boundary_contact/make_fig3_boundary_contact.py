"""Divergent urban boundaries and high spatial uncertainty in global urban land forecasts

Author: Rakibul Ahasan
"""
import argparse
import json
import os

import matplotlib
matplotlib.use("Agg")
import cartopy.crs as ccrs
import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rasterio
from cartopy.io import shapereader
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle
from rasterio.enums import Resampling
from shapely.geometry import shape

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
FIG_DATA = os.path.join(REPO, "data", "figures")

plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
    "font.size": 8, "axes.labelsize": 8, "xtick.labelsize": 8, "ytick.labelsize": 8,
    "text.color": "black", "axes.labelcolor": "black", "xtick.color": "black", "ytick.color": "black",
    "axes.linewidth": 0.6, "pdf.fonttype": 42,
})

X0, Y0 = -17367529.34, 7289708.89
FACTOR = 20
REGION = {1: "Asia-Pacific Developed", 2: "Australia & NZ", 3: "Caribbean", 4: "Developing Pacific",
          5: "Eastern Africa", 6: "Eastern Asia", 7: "Eurasia", 9: "India & Sri Lanka", 10: "Meso America",
          11: "Middle East", 12: "North Africa", 13: "N. & W. Europe", 14: "Rest of S. Asia",
          15: "South America", 16: "South-East Asia", 17: "S. & E. Europe", 18: "S. & Middle Africa",
          19: "USA & Canada", 20: "Western Africa"}
CONTESTED = (5, 20, 14)
CROP_RGB = np.array([[0.94, 0.94, 0.94], [0.97, 0.74, 0.47], [0.88, 0.43, 0.10]])
FOREST_RGB = np.array([[1, 1, 1], [0.72, 0.88, 0.76], [0.38, 0.70, 0.47]])
BIVARIATE = FOREST_RGB[:, None, :] * CROP_RGB[None, :, :]
GROUPS = {"Cropland": ["Cropland"], "Mosaic cropland": ["Mosaic cropland"], "Forest": ["Forest"],
          "Grass and shrub": ["Grassland", "Shrubland", "Mosaic natural"], "Baseline urban": ["Baseline urban"],
          "Other": ["Bare areas", "Flooded vegetation", "Lichens & mosses", "Permanent snow & ice",
                    "Sparse vegetation", "Water bodies"]}
GROUP_COLORS = {"Cropland": "#E0812B", "Mosaic cropland": "#F3C38C", "Forest": "#2E8B57",
                "Grass and shrub": "#A9CF8F", "Baseline urban": "#7F7F7F", "Other": "#CFCFCF"}


def out_name(folder, stem, ext):
    path = os.path.join(folder, f"{stem}.{ext}")
    k = 0
    while os.path.exists(path):
        k += 1
        path = os.path.join(folder, f"{stem}_r{k}.{ext}")
    return path


def main():
    ap = argparse.ArgumentParser(description="Fig. 3, cropland and forest shares of the urban boundary")
    ap.add_argument("--country_contact", default=os.path.join(FIG_DATA, "countries", "country_contact_SSP1_2100.csv"))
    ap.add_argument("--country_codes", default=os.path.join(FIG_DATA, "countries", "country_codes.json"))
    ap.add_argument("--country_grid", default=os.path.join(FIG_DATA, "countries", "countries_1km.tif"))
    ap.add_argument("--regions", default=os.path.join(FIG_DATA, "regions", "regions_display.geojson"))
    ap.add_argument("--regional_shares", default=os.path.join(FIG_DATA, "boundary_composition", "regional_boundary_shares_SSP1_2100.csv"))
    ap.add_argument("--land", default=os.path.join(FIG_DATA, "naturalearth", "ne_110m_land.shp"))
    ap.add_argument("--out", default=HERE)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    code = json.load(open(a.country_codes))["code"]
    t = pd.read_csv(a.country_contact)
    idx = np.array([code[a3] for a3 in t.a3])
    crop_class = np.digitize(t.crop, [30, 50])
    forest_class = np.digitize(t.forest, [4, 14])
    lut = np.zeros((300, 4))
    lut[idx, :3] = BIVARIATE[forest_class, crop_class]
    lut[idx, 3] = 1

    with rasterio.open(a.country_grid) as src:
        grid = src.read(1, out_shape=(src.height // FACTOR, src.width // FACTOR), resampling=Resampling.nearest)
    grid = np.clip(grid, 0, 299)
    gy, gx = grid.shape
    xe = X0 + np.arange(gx + 1) * FACTOR * 1000
    ye = Y0 - np.arange(gy + 1) * FACTOR * 1000

    land = [g for g in shapereader.Reader(a.land).geometries() if g.centroid.y > -60]
    region_geom = {f["properties"]["region_id"]: shape(f["geometry"])
                   for f in json.load(open(a.regions))["features"]}

    def base(ax, shadow=True):
        ax.set_extent([-168, 180, -57, 79], crs=ccrs.PlateCarree())
        if shadow:
            ax.add_geometries(land, ccrs.PlateCarree(), facecolor="white", edgecolor="none", zorder=1,
                              path_effects=[pe.withSimplePatchShadow(offset=(1.3, -1.3), shadow_rgbFace="#8E9499",
                                                                     alpha=0.5, rho=1.0)])
        ax.add_geometries(land, ccrs.PlateCarree(), facecolor="#C9C9C9", edgecolor="none", zorder=2)
        ax.spines["geo"].set_visible(False)

    fig = plt.figure(figsize=(7.2, 8.2))

    # (a) bivariate map
    ax = fig.add_axes([0.0, 0.56, 1.0, 0.44], projection=ccrs.Robinson())
    base(ax)
    ax.imshow(lut[grid], extent=[xe[0], xe[-1], ye[-1], ye[0]], transform=ccrs.epsg(6933), origin="upper",
              interpolation="nearest", zorder=3, regrid_shape=(1400, 3600))
    for r in CONTESTED:
        ax.add_geometries([region_geom[r]], ccrs.PlateCarree(), facecolor="none", edgecolor="#F5D400",
                          linewidth=1.1, zorder=5)
    fig.text(0.006, 0.99, "a", fontsize=10, fontweight="bold", va="top")
    key = fig.add_axes([0.075, 0.595, 0.085, 0.085 * 7.2 / 8.2])
    for i in range(3):
        for j in range(3):
            key.add_patch(Rectangle((j, i), 1, 1, fc=BIVARIATE[i, j], ec="white", lw=0.6))
    key.set_xlim(0, 3)
    key.set_ylim(0, 3)
    key.set_xticks([0, 1, 2, 3])
    key.set_xticklabels(["0", "30", "50", "100"])
    key.set_yticks([0, 1, 2, 3])
    key.set_yticklabels(["0", "4", "14", "100"])
    key.tick_params(length=2, pad=1)
    for s in key.spines.values():
        s.set_visible(False)
    key.set_xlabel("Cropland share (%)", labelpad=1)
    key.set_ylabel("Forest share (%)", labelpad=1)
    ax.legend(handles=[Line2D([], [], color="#F5D400", lw=1.4, label="Regions with the most\ncontested clusters"),
                       Patch(fc="#C9C9C9", label="No measured boundary")],
              loc="lower left", bbox_to_anchor=(0.33, 0.0), frameon=False, fontsize=8, handlelength=1.4)

    # (b, c) single-class maps
    panels = [("b", "crop", "Oranges", "Cropland share (%)"), ("c", "forest", "Greens", "Forest share (%)")]
    for j, (letter, column, cmap_name, label) in enumerate(panels):
        v = np.full(300, np.nan)
        v[idx] = t[column].to_numpy()
        cmap = ListedColormap([plt.get_cmap(cmap_name)(x) for x in np.linspace(0.18, 0.95, 5)])
        norm = BoundaryNorm([0, 20, 40, 60, 80, 100], 5)
        axs = fig.add_axes([j * 0.505, 0.315, 0.495, 0.235], projection=ccrs.Robinson())
        base(axs, shadow=False)
        axs.pcolormesh(xe, ye, np.ma.masked_invalid(v[grid]), transform=ccrs.epsg(6933), cmap=cmap, norm=norm,
                       shading="flat", rasterized=True, zorder=3)
        fig.text(0.006 + j * 0.505, 0.555, letter, fontsize=10, fontweight="bold", va="top")
        cax = fig.add_axes([0.035 + j * 0.505, 0.345, 0.15, 0.011])
        cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax, orientation="horizontal",
                          ticks=[0, 20, 40, 60, 80, 100])
        cb.outline.set_linewidth(0.4)
        cb.ax.tick_params(length=2)
        cb.set_label(label, labelpad=2)
        cb.ax.xaxis.set_label_position("top")
        cb.set_ticks([0, 50, 100])

    # (d) regional composition
    r = pd.read_csv(a.regional_shares).set_index("region_id")
    p = pd.DataFrame({g: r[cols].sum(1) for g, cols in GROUPS.items()})
    p = 100 * p.div(p.sum(1), axis=0)
    p = p.loc[[i for i in REGION if i in p.index]].sort_values("Cropland").iloc[::-1]
    axd = fig.add_axes([0.085, 0.15, 0.905, 0.11])
    bottom = np.zeros(len(p))
    x = np.arange(len(p))
    for g in GROUPS:
        axd.bar(x, p[g], bottom=bottom, color=GROUP_COLORS[g], width=0.78, lw=0, label=g)
        bottom += p[g].to_numpy()
    axd.set_xticks(x)
    axd.set_xticklabels([REGION[i] for i in p.index], rotation=90)
    axd.set_xlim(-0.6, len(p) - 0.4)
    axd.set_ylim(0, 100)
    axd.set_yticks([0, 50, 100])
    axd.set_ylabel("Share of\nboundary (%)")
    axd.spines[["top", "right"]].set_visible(False)
    axd.tick_params(length=2)
    fig.text(0.006, 0.30, "d", fontsize=10, fontweight="bold", va="top")
    axd.legend(loc="lower center", bbox_to_anchor=(0.5, 1.02), ncol=6, frameon=False, fontsize=8,
               handlelength=0.9, columnspacing=0.9, handletextpad=0.3)

    stem = os.path.basename(HERE)
    for ext in ("png", "pdf"):
        path = out_name(a.out, stem, ext)
        fig.savefig(path, dpi=300, facecolor="white")
        print("wrote", path)
    plt.close(fig)


if __name__ == "__main__":
    main()
