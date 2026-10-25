"""Divergent urban boundaries and high spatial uncertainty in global urban land forecasts

Author: Rakibul Ahasan
"""
import argparse
import itertools
import json
import os

import matplotlib
matplotlib.use("Agg")
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import numpy as np
from cartopy.io import shapereader
from matplotlib.colors import LinearSegmentedColormap, Normalize

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
BLOCK = 20000
EARTH_RADIUS = 6378137.0
COS30 = np.cos(np.radians(30.0))
PAIRS = list(itertools.combinations(range(9), 2))

COLORS = {"low": "#D6336C", "mid": "#E8A317", "high": "#3B6FB6"}
REGION = {1: "Asia-Pacific Developed", 2: "Australia & NZ", 3: "Caribbean", 4: "Developing Pacific",
          5: "Eastern Africa", 6: "Eastern Asia", 7: "Eurasia", 9: "India & Sri Lanka", 10: "Meso America",
          11: "Middle East", 12: "North Africa", 13: "N. & W. Europe", 14: "Rest of S. Asia",
          15: "South America", 16: "South-East Asia", 17: "S. & E. Europe", 18: "S. & Middle Africa",
          19: "USA & Canada", 20: "Western Africa"}
SHORT_COUNTRY = {"United States of America": "United States", "Dem. Rep. Congo": "DR Congo"}


def c2(n):
    return n * (n - 1) / 2.0


def x_to_lon(x):
    return np.degrees(x / (EARTH_RADIUS * COS30))


def y_to_lat(y):
    return np.degrees(np.arcsin(np.clip(y * COS30 / EARTH_RADIUS, -1, 1)))


def pair_jaccard(z, year, unit):
    inter = z[f"inter{year}_{unit}"].astype(float)
    area = z[f"area{year}_{unit}"].astype(float)
    union = np.stack([area[:, i] + area[:, j] for i, j in PAIRS], 1) - inter
    return np.where(union > 0, inter / np.maximum(union, 1), np.nan), inter, union


def out_name(folder, stem, ext):
    path = os.path.join(folder, f"{stem}.{ext}")
    k = 0
    while os.path.exists(path):
        k += 1
        path = os.path.join(folder, f"{stem}_r{k}.{ext}")
    return path


def main():
    ap = argparse.ArgumentParser(description="Fig. 4, pairwise agreement among forecasts")
    ap.add_argument("--agreement", default=os.path.join(FIG_DATA, "pairwise_agreement", "pairwise_agreement_SSP1.npz"))
    ap.add_argument("--country_codes", default=os.path.join(FIG_DATA, "countries", "country_codes.json"))
    ap.add_argument("--land", default=os.path.join(FIG_DATA, "naturalearth", "ne_110m_land.shp"))
    ap.add_argument("--out", default=HERE)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    z = np.load(a.agreement)
    n = np.arange(10)

    # (a) block index and profiles, km2 as means over the 36 pairs
    h = z["hist2100"].astype(float)
    inter = (h * c2(n)).sum(-1) / 36
    union = (h * (36 - c2(9 - n))).sum(-1) / 36
    jb = np.where(union >= 5, inter / np.maximum(union, 1e-9), np.nan)
    gy, gx = jb.shape
    xe = X0 + np.arange(gx + 1) * BLOCK
    ye = Y0 - np.arange(gy + 1) * BLOCK
    lon = x_to_lon(X0 + (np.arange(gx) + 0.5) * BLOCK)
    lat = y_to_lat(Y0 - (np.arange(gy) + 0.5) * BLOCK)
    lon_bins = np.arange(-180, 181, 1.0)
    lat_bins = np.arange(-60, 85, 1.0)
    parts = {"shared": inter, "single": union - inter}
    lon_prof = {k: np.histogram(lon, lon_bins, weights=v.sum(0))[0] / 1e3 for k, v in parts.items()}
    lat_prof = {k: np.histogram(lat, lat_bins, weights=v.sum(1))[0] / 1e3 for k, v in parts.items()}
    part_color = {"shared": COLORS["high"], "single": COLORS["low"]}
    part_label = {"shared": "Urban in both products of a pair", "single": "Urban in one product of a pair"}

    h2020 = z["hist2020_global"].astype(float)
    j2020 = (h2020 * c2(n)).sum() / (h2020 * (36 - c2(9 - n))).sum()
    j2100 = inter.sum() / union.sum()
    print(f"Pooled pairwise Jaccard index, SSP1: 2020 {j2020:.4f}, 2100 {j2100:.4f}")

    cmap = LinearSegmentedColormap.from_list("j", [COLORS["low"], COLORS["mid"], COLORS["high"]])
    norm = Normalize(0, 0.8)

    # (b) regions and (c) countries
    jr21, ir21, ur21 = pair_jaccard(z, 2100, "region")
    _, ir20, ur20 = pair_jaccard(z, 2020, "region")
    pr21 = ir21.sum(1) / np.maximum(ur21.sum(1), 1)
    pr20 = ir20.sum(1) / np.maximum(ur20.sum(1), 1)
    regions = sorted([r for r in REGION if ur21[r].sum() > 0], key=lambda r: -pr21[r])
    _, ic21, uc21 = pair_jaccard(z, 2100, "country")
    _, ic20, uc20 = pair_jaccard(z, 2020, "country")
    pc21 = ic21.sum(1) / np.maximum(uc21.sum(1), 1)
    pc20 = ic20.sum(1) / np.maximum(uc20.sum(1), 1)
    cc = json.load(open(a.country_codes))
    code_to_a3 = {v: k for k, v in cc["code"].items()}
    area = z["area2100_country"].sum(1)
    top = [u for u in np.argsort(-area) if u > 0][:30]
    top = sorted(top, key=lambda u: -pc21[u])

    fig = plt.figure(figsize=(7.2, 7.0))
    left, width, bottom = 0.085, 0.79, 0.48
    height = width * 7.2 / (360 / 137) / 7.0
    ax = fig.add_axes([left, bottom, width, height], projection=ccrs.PlateCarree())
    ax.set_extent([-180, 180, -57, 80], crs=ccrs.PlateCarree())
    land = [g for g in shapereader.Reader(a.land).geometries() if g.centroid.y > -60]
    ax.add_geometries(land, ccrs.PlateCarree(), facecolor="#E9E9E9", edgecolor="none", zorder=1)
    ax.pcolormesh(xe, ye, np.ma.masked_invalid(jb), transform=ccrs.epsg(6933), cmap=cmap, norm=norm,
                  shading="flat", rasterized=True, zorder=3)
    ax.spines["geo"].set_linewidth(0.6)
    ax.set_xticks(range(-150, 181, 50), crs=ccrs.PlateCarree())
    ax.set_yticks([-40, 0, 40, 80], crs=ccrs.PlateCarree())
    ax.tick_params(length=2)
    cax = fig.add_axes([left + 0.02, bottom + 0.05, 0.17, 0.012])
    cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax, orientation="horizontal", ticks=[0, 0.4, 0.8])
    cb.outline.set_linewidth(0.4)
    cb.ax.tick_params(length=2)
    cb.ax.xaxis.set_label_position("top")
    cb.set_label("Pairwise Jaccard index", labelpad=2)
    fig.text(0.006, bottom + height + 0.132, "a", fontsize=10, fontweight="bold", va="top")

    axt = fig.add_axes([left, bottom + height + 0.012, width, 0.12])
    base = np.zeros(len(lon_bins) - 1)
    xm = (lon_bins[:-1] + lon_bins[1:]) / 2
    for k in ("shared", "single"):
        axt.fill_between(xm, base, base + lon_prof[k], color=part_color[k], lw=0, step="mid", label=part_label[k])
        base = base + lon_prof[k]
    axt.set_xlim(-180, 180)
    axt.set_ylim(0, base.max() * 1.05)
    axt.set_xticks([])
    axt.spines[["top", "right", "bottom"]].set_visible(False)
    axt.set_ylabel("10³ km² per\ndegree")
    axt.legend(loc="upper left", frameon=False, handlelength=1.0, borderaxespad=0, labelspacing=0.2)

    axr = fig.add_axes([left + width + 0.012, bottom, 0.09, height])
    base = np.zeros(len(lat_bins) - 1)
    ym = (lat_bins[:-1] + lat_bins[1:]) / 2
    for k in ("shared", "single"):
        axr.fill_betweenx(ym, base, base + lat_prof[k], color=part_color[k], lw=0, step="mid")
        base = base + lat_prof[k]
    axr.set_ylim(-57, 80)
    axr.set_xlim(0, base.max() * 1.05)
    axr.set_yticks([])
    step = 50 if base.max() > 100 else 25
    axr.set_xticks(np.arange(0, base.max(), step)[:2])
    axr.spines[["top", "right", "left"]].set_visible(False)
    axr.set_xlabel("10³ km² per\ndegree")

    axb = fig.add_axes([0.085, 0.215, 0.36, 0.185])
    x = np.arange(len(regions))
    bp = axb.boxplot([jr21[r][~np.isnan(jr21[r])] for r in regions], positions=x, widths=0.6, whis=(0, 100),
                     patch_artist=True, showfliers=False, medianprops=dict(color="black", lw=0.8),
                     boxprops=dict(lw=0.5), whiskerprops=dict(lw=0.5), capprops=dict(lw=0.5))
    for box, r in zip(bp["boxes"], regions):
        box.set_facecolor(cmap(norm(pr21[r])))
        box.set_edgecolor("#444444")
    axb.scatter(x, pr20[regions], marker="D", s=10, facecolor="white", edgecolor="black", lw=0.6, zorder=4,
                label="Pooled index, 2020")
    axb.set_xticks(x)
    axb.set_xticklabels([REGION[r] for r in regions], rotation=90)
    axb.set_xlim(-0.6, len(regions) - 0.4)
    axb.set_ylim(0, 0.9)
    axb.set_ylabel("Pairwise Jaccard index")
    axb.spines[["top", "right"]].set_visible(False)
    axb.tick_params(length=2)
    axb.legend(loc="upper right", frameon=False, handletextpad=0.2, borderaxespad=0)
    fig.text(0.006, 0.425, "b", fontsize=10, fontweight="bold", va="top")

    axc = fig.add_axes([0.545, 0.215, 0.445, 0.185])
    x = np.arange(len(top))
    axc.vlines(x, pc21[top], pc20[top], color="#9A9A9A", lw=0.8)
    axc.scatter(x, pc20[top], s=12, facecolor="white", edgecolor="black", lw=0.6, zorder=3, label="2020")
    axc.scatter(x, pc21[top], s=14, c=[cmap(norm(v)) for v in pc21[top]], edgecolor="black", lw=0.4, zorder=4, label="2100")
    names = [cc["name"][code_to_a3[u]] for u in top]
    axc.set_xticks(x)
    axc.set_xticklabels([SHORT_COUNTRY.get(s, s) for s in names], rotation=90)
    axc.set_xlim(-0.6, len(top) - 0.4)
    axc.set_ylim(0, 0.9)
    axc.set_ylabel("Pooled Jaccard index")
    axc.spines[["top", "right"]].set_visible(False)
    axc.tick_params(length=2)
    axc.legend(loc="upper right", frameon=False, handletextpad=0.2, borderaxespad=0, ncol=2, columnspacing=0.8)
    fig.text(0.47, 0.425, "c", fontsize=10, fontweight="bold", va="top")

    stem = os.path.basename(HERE)
    for ext in ("png", "pdf"):
        path = out_name(a.out, stem, ext)
        fig.savefig(path, dpi=300, facecolor="white")
        print("wrote", path)
    plt.close(fig)


if __name__ == "__main__":
    main()
