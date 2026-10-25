"""Divergent urban boundaries and high spatial uncertainty in global urban land forecasts

Author: Rakibul Ahasan
"""
import argparse
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import numpy as np
import cartopy.crs as ccrs
from cartopy.io import shapereader
from matplotlib.patches import Polygon as MplPolygon
from pyproj import Transformer
from shapely.geometry import box, shape
from shapely.ops import transform as shapely_transform

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
FIG_DATA = os.path.join(REPO, "data", "figures")
NAME = "SFig11_disagreement_timing"

plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
    "font.size": 8, "text.color": "black", "pdf.fonttype": 42,
})

SSPS = ["SSP1", "SSP3", "SSP5"]
GRID_X0, GRID_Y0 = -17367529.34, 7289708.89
BLOCK = 5000.0
HALF = 110000.0
KEYS = ["d2020", "d2050", "d2100"]
LABELS = ["Disputed\nby 2020", "First disputed\n2020 to 2050", "First disputed\n2050 to 2100"]
COLORS = np.array([[0.23, 0.43, 0.74], [0.91, 0.66, 0.12], [0.84, 0.20, 0.42]])
WINDOWS = [("b", "New York–Philadelphia", -74.6, 40.25), ("c", "Lagos–Ibadan", 3.6, 6.95),
           ("d", "Rhine–Ruhr", 7.0, 51.4), ("e", "Dhaka", 90.4, 23.78)]

PANEL_W, PANEL_H = 7.2, 5.35
GAP = 0.2
FIG_H = 3 * PANEL_H + 2 * GAP

TO_6933 = Transformer.from_crs(4326, 6933, always_xy=True)
TO_LONLAT = Transformer.from_crs(6933, 4326, always_xy=True)


def out_name(folder, stem, ext):
    path = os.path.join(folder, f"{stem}.{ext}")
    k = 1
    while os.path.exists(path):
        path = os.path.join(folder, f"{stem}_r{k}.{ext}")
        k += 1
    return path


def page_rect(k, x, y, w, h):
    bottom = FIG_H - k * (PANEL_H + GAP) - PANEL_H
    return [x, (bottom + y * PANEL_H) / FIG_H, w, h * PANEL_H / FIG_H]


def aggregate(a, f):
    h, w = a.shape
    a = np.pad(a, ((0, (-h) % f), (0, (-w) % f)))
    return a.reshape(a.shape[0] // f, f, a.shape[1] // f, f).sum(axis=(1, 3))


def mix(shares):
    return np.clip(shares @ COLORS, 0, 1)


def window_slice(lon, lat):
    cx, cy = TO_6933.transform(lon, lat)
    x0, x1, y0, y1 = cx - HALF, cx + HALF, cy - HALF, cy + HALF
    c0 = int((x0 - GRID_X0) // BLOCK)
    c1 = int((x1 - GRID_X0) // BLOCK)
    r0 = int((GRID_Y0 - y1) // BLOCK)
    r1 = int((GRID_Y0 - y0) // BLOCK)
    return (x0, x1, y0, y1), (r0, r1, c0, c1)


def nice(x):
    p = 10 ** np.floor(np.log10(x))
    return float(np.round(x / p) * p)


def dot_caps(layers):
    pooled = {"global": [], **{name: [] for _, name, _, _ in WINDOWS}}
    for p in SSPS:
        comps = [layers[p][k] for k in KEYS]
        total = sum(aggregate(c, 60) for c in comps)
        pooled["global"].append(total[total > 0])
        for _, name, lon, lat in WINDOWS:
            _, (r0, r1, c0, c1) = window_slice(lon, lat)
            t = sum(c[r0:r1, c0:c1] for c in comps)
            pooled[name].append(t[t > 0])
    return {k: float(np.percentile(np.concatenate(v), 98)) for k, v in pooled.items()}


def draw_page(fig, k, comps, caps, land, borders):
    totals = np.array([c.sum() for c in comps])
    print(SSPS[k], "disputed by date (km²):", totals.round().tolist(), "shares:", (totals / totals.sum()).round(3).tolist())

    # a: global map, 300-km lattice
    grids = [aggregate(c, 60) for c in comps]
    total = sum(grids)
    iy, ix = np.nonzero(total > 0)
    lon, lat = TO_LONLAT.transform(GRID_X0 + (ix + 0.5) * BLOCK * 60, GRID_Y0 - (iy + 0.5) * BLOCK * 60)
    lon, lat = np.asarray(lon), np.asarray(lat)
    keep = lat > -57
    shares = np.stack([g[iy, ix] for g in grids], -1) / total[iy, ix][:, None]
    cap = caps["global"]
    dmax = 7.0
    size = (dmax * np.sqrt(np.minimum(total[iy, ix], cap) / cap)) ** 2

    ax = fig.add_axes(page_rect(k, 0.0, 0.40, 1.0, 0.60), projection=ccrs.Robinson())
    ax.set_extent([-168, 180, -57, 79], crs=ccrs.PlateCarree())
    ax.add_geometries(land, ccrs.PlateCarree(), facecolor="white", edgecolor="none", zorder=1,
                      path_effects=[pe.withSimplePatchShadow(offset=(1.5, -1.5), shadow_rgbFace="#8E9499",
                                                             alpha=0.55, rho=1.0)])
    ax.add_geometries(land, ccrs.PlateCarree(), facecolor="#F4F4F4", edgecolor="#BDBDBD", linewidth=0.25, zorder=2)
    ax.gridlines(draw_labels=False, linewidth=0.25, color="#E2E2E2", zorder=3)
    order = np.argsort(-size[keep])
    ax.scatter(lon[keep][order], lat[keep][order], s=size[keep][order], c=mix(shares[keep][order]),
               linewidths=0, transform=ccrs.PlateCarree(), zorder=5)
    ax.spines["geo"].set_edgecolor("#BDBDBD")
    ax.spines["geo"].set_linewidth(0.5)
    for letter, _, wlon, wlat in WINDOWS:
        (x0, x1, y0, y1), _ = window_slice(wlon, wlat)
        bx, by = TO_LONLAT.transform([x0, x1, x1, x0, x0], [y0, y0, y1, y1, y0])
        ax.plot(bx, by, color="black", lw=0.8, transform=ccrs.PlateCarree(), zorder=7)
        ax.text(wlon + 3.5, wlat + 3.0, letter, fontsize=9, fontweight="bold", transform=ccrs.PlateCarree(), zorder=8,
                path_effects=[pe.withStroke(linewidth=2, foreground="white")])
    x, y, _, _ = page_rect(k, 0.006, 0.985, 0, 0)
    fig.text(x, y, "a", fontsize=10, fontweight="bold", va="top")

    key = fig.add_axes(page_rect(k, 0.035, 0.455, 0.155, 0.17))
    key.set_axis_off()
    n = 200
    u, v = np.meshgrid(np.linspace(0, 1, n), np.linspace(0, 1, n))
    s3 = v / (np.sqrt(3) / 2)
    s2 = u - 0.5 * s3
    s1 = 1 - s2 - s3
    img = mix(np.clip(np.stack([s1, s2, s3], -1), 0, 1))
    alpha = ((s1 >= 0) & (s2 >= 0) & (s3 >= 0) & (v <= np.sqrt(3) / 2)).astype(float)
    key.imshow(np.dstack([img, alpha]), origin="lower", extent=[0, 1, 0, 1], interpolation="bilinear")
    key.set_xlim(-0.05, 1.05)
    key.set_ylim(-0.05, 1.0)
    key.set_aspect("equal")
    gsh = totals / totals.sum()
    key.plot(gsh[1] + 0.5 * gsh[2], gsh[2] * np.sqrt(3) / 2, "o", ms=4.5, mfc="white", mec="black", mew=0.8)
    key.text(0.0, -0.06, LABELS[0], ha="center", va="top", fontsize=8)
    key.text(1.0, -0.06, LABELS[1], ha="center", va="top", fontsize=8)
    key.text(0.5, 0.92, LABELS[2], ha="center", va="bottom", fontsize=8)

    sk = fig.add_axes(page_rect(k, 0.62, 0.415, 0.20, 0.07))
    sk.set_axis_off()
    sk.set_xlim(0, 3)
    sk.set_ylim(0, 1)
    for i, value in enumerate([nice(cap * 0.05), nice(cap * 0.3), nice(cap)]):
        sk.scatter(i + 0.5, 0.62, s=(dmax * np.sqrt(min(value, cap) / cap)) ** 2, color="#9A9A9A", linewidths=0)
        sk.text(i + 0.5, 0.05, (f"≥{value:,.0f}" if i == 2 else f"{value:,.0f}"), ha="center", va="bottom", fontsize=8)
    sk.text(1.5, 1.12, "Area per dot (km²)", ha="center", va="bottom", fontsize=8)

    # b to e: 5-km dots in 220 x 220 km windows
    width = 0.236
    gap = (1 - 4 * width) / 5
    height = width * PANEL_W / PANEL_H
    for j, (letter, name, wlon, wlat) in enumerate(WINDOWS):
        (x0, x1, y0, y1), (r0, r1, c0, c1) = window_slice(wlon, wlat)
        sub = [c[r0:r1, c0:c1] for c in comps]
        t = sum(sub)
        ry, rx = np.nonzero(t > 0)
        px = GRID_X0 + (c0 + rx + 0.5) * BLOCK
        py = GRID_Y0 - (r0 + ry + 0.5) * BLOCK
        sh = np.stack([q[ry, rx] for q in sub], -1) / t[ry, rx][:, None]
        capw = caps[name]
        ds = 2.9 * np.sqrt(np.minimum(t[ry, rx], capw) / capw)
        left = gap + j * (width + gap)
        axw = fig.add_axes(page_rect(k, left, 0.035, width, height))
        axw.set_facecolor("#EEF2F5")
        frame = box(x0, y0, x1, y1)
        for g in borders:
            if not g.intersects(frame):
                continue
            part = g.intersection(frame)
            polys = [part] if part.geom_type == "Polygon" else list(getattr(part, "geoms", []))
            for pg in polys:
                if pg.geom_type == "Polygon":
                    axw.add_patch(MplPolygon(np.asarray(pg.exterior.coords), closed=True, fc="white",
                                             ec="#B8B8B8", lw=0.4, zorder=1))
        axw.scatter(px, py, s=ds ** 2, c=mix(sh), linewidths=0, zorder=3)
        axw.set_xlim(x0, x1)
        axw.set_ylim(y0, y1)
        axw.set_aspect("equal")
        axw.set_xticks([])
        axw.set_yticks([])
        for sp in axw.spines.values():
            sp.set_edgecolor("#9A9A9A")
            sp.set_linewidth(0.6)
        axw.plot([x1 - 112000, x1 - 12000], [y0 + 12000] * 2, color="black", lw=1.2, zorder=5)
        axw.text(x1 - 62000, y0 + 16000, "100 km", ha="center", va="bottom", fontsize=8, zorder=5,
                 bbox=dict(fc="white", ec="none", pad=1.0, alpha=0.85))
        tx, ty, _, _ = page_rect(k, left, 0.035 + height + 0.008, 0, 0)
        fig.text(tx, ty, letter, fontsize=10, fontweight="bold", va="bottom")
        fig.text(tx + 0.025, ty, name, fontsize=9, fontweight="bold", va="bottom")


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--blocks", default=os.path.join(FIG_DATA, "blocks5km"),
                    help="folder with disagreement_onset_<SSP>.npz")
    ap.add_argument("--borders", default=os.path.join(FIG_DATA, "regions", "world_boundaries.geojson"))
    ap.add_argument("--land", default=os.path.join(FIG_DATA, "naturalearth", "ne_110m_land.shp"))
    ap.add_argument("--out", default=HERE)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    layers = {}
    for p in SSPS:
        z = np.load(os.path.join(a.blocks, f"disagreement_onset_{p}.npz"))
        layers[p] = {k: z[k].astype(float) for k in KEYS}
    caps = dot_caps(layers)
    land = [g for g in shapereader.Reader(a.land).geometries() if g.centroid.y > -60]
    with open(a.borders) as f:
        features = json.load(f)["features"]
    borders = [shapely_transform(lambda x, y, z=None: TO_6933.transform(x, y), shape(ft["geometry"])) for ft in features]

    fig = plt.figure(figsize=(PANEL_W, FIG_H))
    for k, p in enumerate(SSPS):
        draw_page(fig, k, [layers[p][key] for key in KEYS], caps, land, borders)
    print("dot-size caps (km²):", {key: round(v, 1) for key, v in caps.items()})
    for ext in ("png", "pdf"):
        fig.savefig(out_name(a.out, NAME, ext), dpi=300, facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    main()
