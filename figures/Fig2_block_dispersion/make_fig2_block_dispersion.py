"""Figure 2. Block-scale variation in boundary length relative to urban area.

For each 100 x 100 km block in which every forecast maps at least 10 km2 of urban land, R = log2(CV of
boundary length / CV of urban area) across the nine forecasts (population standard deviation over mean).
a  R in 2020 under SSP1, with an inset histogram and the shares on either side of zero.
b  share of blocks with R > 0 by region, 2020 and 2100, SSP1 (regions with at least 20 eligible blocks in both years).
c  R by decile of block urban area (mean across forecasts), 2020 and 2100, SSP1.
d  share of blocks with R > 0 by pathway and year, bars at the 10 km2 floor and dots at a 50 km2 floor.
"""
import argparse
import os

import cartopy.crs as ccrs
import matplotlib
matplotlib.use("Agg")
import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import numpy as np
from cartopy.io import shapereader
from matplotlib.colors import BoundaryNorm, ListedColormap

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
FIG_DATA = os.path.join(REPO, "data", "figures")
NAME = "Fig2_block_dispersion"

GRID_X0, GRID_Y0 = -17367529.34, 7289708.89     # upper-left corner of the 1-km EPSG:6933 analysis grid (m)
BLOCK = 100000                                  # block size (m)
REGION = {1: "Asia-Pacific Developed", 2: "Australia & NZ", 3: "Caribbean", 4: "Developing Pacific",
          5: "Eastern Africa", 6: "Eastern Asia", 7: "Eurasia", 9: "India & Sri Lanka", 10: "Meso America",
          11: "Middle East", 12: "North Africa", 13: "N. & W. Europe", 14: "Rest of S. Asia", 15: "South America",
          16: "South-East Asia", 17: "S. & E. Europe", 18: "S. & Middle Africa", 19: "USA & Canada",
          20: "Western Africa"}
BREAKS = [-9, -1, -0.5, -0.15, 0.15, 0.5, 1, 9]
COLORS = ["#2B5C9E", "#5F8FC9", "#A9C4E4", "#E9E9E9", "#F2B3C6", "#E2678F", "#B0244F"]

plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
                     "font.size": 8, "pdf.fonttype": 42, "axes.linewidth": 0.6,
                     "mathtext.default": "regular"})


def out_path(folder, ext):
    path = os.path.join(folder, f"{NAME}.{ext}")
    k = 1
    while os.path.exists(path):
        path = os.path.join(folder, f"{NAME}_r{k}.{ext}")
        k += 1
    return path


def block_ratio(blocks_dir, ssp, year, floor=10):
    """log2 ratio of the CV of boundary length to the CV of urban area, and the mean urban area, per block."""
    z = np.load(os.path.join(blocks_dir, f"blocks100km_{ssp}_{year}.npz"))
    area = z["area"].astype(float)
    boundary = z["bnd"].astype(float)
    eligible = area.min(0) >= floor
    with np.errstate(all="ignore"):
        r = np.log2((boundary.std(0) / boundary.mean(0)) / (area.std(0) / area.mean(0)))
    r[~eligible] = np.nan
    return r, area.mean(0)


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--blocks", default=os.path.join(FIG_DATA, "blocks100km"),
                    help="folder with blocks100km_<SSP>_<year>.npz and region_ids_100km.npy")
    ap.add_argument("--land", default=os.path.join(FIG_DATA, "naturalearth", "ne_110m_land.shp"))
    ap.add_argument("--out", default=HERE)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    r20, m20 = block_ratio(a.blocks, "SSP1", 2020)
    r21, m21 = block_ratio(a.blocks, "SSP1", 2100)
    nby, nbx = r20.shape
    regions = np.zeros((nby, nbx), int)
    g = np.load(os.path.join(a.blocks, "region_ids_100km.npy"))
    regions[:g.shape[0], :g.shape[1]] = g[:nby, :nbx]
    cmap = ListedColormap(COLORS)
    norm = BoundaryNorm(BREAKS, 7)

    fig = plt.figure(figsize=(7.2, 6.95))

    # a map of R in 2020
    ax = fig.add_axes([0.0, 0.50, 1.0, 0.50], projection=ccrs.Robinson())
    ax.set_extent([-168, 180, -57, 79], crs=ccrs.PlateCarree())
    land = [geom for geom in shapereader.Reader(a.land).geometries() if geom.centroid.y > -60]
    ax.add_geometries(land, ccrs.PlateCarree(), facecolor="white", edgecolor="none", zorder=1,
                      path_effects=[pe.withSimplePatchShadow(offset=(1.3, -1.3), shadow_rgbFace="#8E9499",
                                                             alpha=0.5, rho=1.0)])
    ax.add_geometries(land, ccrs.PlateCarree(), facecolor="#F4F4F4", edgecolor="none", zorder=2)
    ax.spines["geo"].set_visible(False)
    xe = GRID_X0 + np.arange(nbx + 1) * BLOCK
    ye = GRID_Y0 - np.arange(nby + 1) * BLOCK
    ax.pcolormesh(xe, ye, np.ma.masked_invalid(r20), transform=ccrs.epsg(6933), cmap=cmap, norm=norm,
                  shading="flat", rasterized=True, zorder=3)
    fig.text(0.006, 0.995, "a", fontsize=10, fontweight="bold", va="top")

    v = r20[np.isfinite(r20)]
    counts = np.histogram(np.clip(v, -8, 8), BREAKS)[0]
    share = 100 * counts / counts.sum()
    hx = fig.add_axes([0.075, 0.575, 0.19, 0.115])
    hx.bar(np.arange(7), share, color=COLORS, edgecolor="#9A9A9A", lw=0.3, width=0.8)
    hx.axvline(3, color="black", lw=0.6, ls=(0, (2, 2)))
    hx.set_xticks(np.arange(-0.5, 7, 1))
    hx.set_xticklabels(["", "−1", "−0.5", "−0.15", "0.15", "0.5", "1", ""], rotation=90)
    hx.set_ylabel("Blocks (%)")
    hx.spines[["top", "right"]].set_visible(False)
    hx.tick_params(length=2)
    hx.set_xlim(-0.6, 6.6)
    neg = 100 * (v < 0).mean()
    pos = 100 * (v > 0).mean()
    top = hx.get_ylim()[1] * 1.04
    hx.text(1.0, top, f"{neg:.0f}%\narea more\nvariable", ha="center", va="bottom", color=COLORS[0])
    hx.text(5.2, top, f"{pos:.0f}%\nboundary\nmore variable", ha="center", va="bottom", color=COLORS[6])
    hx.set_xlabel("log$_2$(CV boundary / CV area)", labelpad=1)

    # b share of blocks with R > 0 by region
    def share_positive(r, rid):
        return 100 * (r[np.isfinite(r) & (regions == rid)] > 0).mean()

    rids = [rid for rid in REGION
            if (np.isfinite(r20) & (regions == rid)).sum() >= 20 and (np.isfinite(r21) & (regions == rid)).sum() >= 20]
    s20 = {rid: share_positive(r20, rid) for rid in rids}
    s21 = {rid: share_positive(r21, rid) for rid in rids}
    rids = sorted(rids, key=lambda rid: s20[rid])
    axb = fig.add_axes([0.205, 0.07, 0.19, 0.37])
    y = np.arange(len(rids))
    axb.barh(y + 0.2, [s20[rid] for rid in rids], height=0.38, color=COLORS[6], lw=0, label="2020")
    axb.barh(y - 0.2, [s21[rid] for rid in rids], height=0.38, color=COLORS[4], lw=0, label="2100")
    axb.axvline(50, color="black", lw=0.6, ls=(0, (2, 2)))
    axb.set_yticks(y)
    axb.set_yticklabels([REGION[rid] for rid in rids])
    axb.set_xlim(0, 100)
    axb.set_xticks([0, 50, 100])
    axb.set_xlabel("Blocks with boundary\nmore variable (%)")
    axb.spines[["top", "right"]].set_visible(False)
    axb.tick_params(length=2)
    axb.set_ylim(-0.6, len(rids) - 0.4)
    axb.legend(loc="lower right", frameon=False, handlelength=1.0, borderaxespad=0)
    fig.text(0.006, 0.465, "b", fontsize=10, fontweight="bold", va="top")

    # c R by decile of block urban area
    axc = fig.add_axes([0.475, 0.11, 0.25, 0.33])
    for r, m, col, off in ((r20, m20, COLORS[6], -0.2), (r21, m21, COLORS[4], 0.2)):
        ok = np.isfinite(r)
        q = np.quantile(m[ok], np.linspace(0, 1, 11))
        dec = np.clip(np.searchsorted(q, m[ok], side="right") - 1, 0, 9)
        data = [r[ok][dec == i] for i in range(10)]
        bp = axc.boxplot(data, positions=np.arange(10) + off, widths=0.34, patch_artist=True, showfliers=False,
                         medianprops=dict(color="black", lw=0.7), boxprops=dict(lw=0.4),
                         whiskerprops=dict(lw=0.4), capprops=dict(lw=0.4))
        for box in bp["boxes"]:
            box.set_facecolor(col)
            box.set_edgecolor("#555555")
    axc.axhline(0, color="black", lw=0.6, ls=(0, (2, 2)))
    axc.set_xticks(np.arange(10))
    axc.set_xticklabels([str(i) for i in range(1, 11)])
    axc.set_xlim(-0.6, 9.6)
    axc.set_xlabel("Decile of block urban area")
    axc.set_ylabel("log$_2$(CV boundary / CV area)")
    axc.spines[["top", "right"]].set_visible(False)
    axc.tick_params(length=2)
    axc.set_ylim(-2, 2.5)
    fig.text(0.415, 0.465, "c", fontsize=10, fontweight="bold", va="top")

    # d pathway and year
    axd = fig.add_axes([0.81, 0.11, 0.17, 0.33])
    shade = {"SSP1": "#4D4D4D", "SSP3": "#8C8C8C", "SSP5": "#C4C4C4"}
    for i, year in enumerate((2020, 2050, 2100)):
        for j, ssp in enumerate(("SSP1", "SSP3", "SSP5")):
            r10, _ = block_ratio(a.blocks, ssp, year)
            r50, _ = block_ratio(a.blocks, ssp, year, 50)
            x = i + (j - 1) * 0.27
            axd.bar(x, 100 * (r10[np.isfinite(r10)] > 0).mean(), width=0.25, color=shade[ssp], lw=0,
                    label=ssp if i == 0 else None)
            axd.scatter(x, 100 * (r50[np.isfinite(r50)] > 0).mean(), s=9, color="white", edgecolor="black",
                        lw=0.6, zorder=3)
    axd.axhline(50, color="black", lw=0.6, ls=(0, (2, 2)))
    axd.set_xticks([0, 1, 2])
    axd.set_xticklabels(["2020", "2050", "2100"])
    axd.set_ylim(0, 100)
    axd.set_yticks([0, 50, 100])
    axd.set_ylabel("Blocks with boundary more variable (%)")
    axd.spines[["top", "right"]].set_visible(False)
    axd.tick_params(length=2)
    axd.legend(loc="upper right", frameon=False, handlelength=0.9, borderaxespad=0, labelspacing=0.2)
    fig.text(0.745, 0.465, "d", fontsize=10, fontweight="bold", va="top")

    for ext in ("png", "pdf"):
        fig.savefig(out_path(a.out, ext), dpi=300, facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    main()
