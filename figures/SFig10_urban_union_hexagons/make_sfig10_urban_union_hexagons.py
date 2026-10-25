"""Divergent urban boundaries and high spatial uncertainty in global urban land forecasts

Author: Rakibul Ahasan
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import cartopy.crs as ccrs
from matplotlib.collections import PatchCollection
from matplotlib.colors import LinearSegmentedColormap, LogNorm
from matplotlib.patches import Circle, RegularPolygon
from pyproj import Transformer

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
FIG_DATA = os.path.join(REPO, "data", "figures")
NAME = "SFig10_urban_union_hexagons"

plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
    "font.size": 8, "text.color": "black", "pdf.fonttype": 42,
})

SSPS = ["SSP1", "SSP3", "SSP5"]
LETTERS = {"SSP1": "a", "SSP3": "b", "SSP5": "c"}
GRID_X0, GRID_Y0 = -17367529.34, 7289708.89
BLOCK = 5000.0
SPACING = 350000.0
SIDE = SPACING / np.sqrt(3)
APOTHEM = SPACING / 2
RMAX = 0.92 * APOTHEM
RING = 0.6
ROB = ccrs.Robinson()

PANEL_W, PANEL_H = 7.2, 5.35
KEEP = 0.61
GAP = 0.2
FIG_H = 3 * KEEP * PANEL_H + 2 * GAP

NORM = LogNorm(10, 30000)
CMAP = LinearSegmentedColormap.from_list("union2020", ["#F2C14E", "#E07A3F", "#D6336C", "#7B2A8C", "#2B1A5E"])


def out_name(folder, stem, ext):
    path = os.path.join(folder, f"{stem}.{ext}")
    k = 1
    while os.path.exists(path):
        path = os.path.join(folder, f"{stem}_r{k}.{ext}")
        k += 1
    return path


def panel_rect(k, x, y, w, h):
    top = FIG_H - k * (KEEP * PANEL_H + GAP)
    y_in = top - (1 - y) * PANEL_H
    return [x, y_in / FIG_H, w, h * PANEL_H / FIG_H]


def hexagon_index(valid):
    iy, ix = np.nonzero(valid > 0)
    to_lonlat = Transformer.from_crs(6933, 4326, always_xy=True)
    lon, lat = to_lonlat.transform(GRID_X0 + (ix + 0.5) * BLOCK, GRID_Y0 - (iy + 0.5) * BLOCK)
    xy = ROB.transform_points(ccrs.PlateCarree(), np.asarray(lon), np.asarray(lat))
    x, y = xy[:, 0], xy[:, 1]
    qf = (np.sqrt(3) / 3 * x - y / 3) / SIDE
    rf = (2 / 3 * y) / SIDE
    sf = -qf - rf
    q, r, s = np.round(qf), np.round(rf), np.round(sf)
    dq, dr, ds = abs(q - qf), abs(r - rf), abs(s - sf)
    fix_q = (dq > dr) & (dq > ds)
    q = np.where(fix_q, -r - s, q)
    r = np.where(~fix_q & (dr > ds), -q - s, r)
    key = (q.astype(np.int64) + 2000) * 4000 + (r.astype(np.int64) + 2000)
    keys, inverse = np.unique(key, return_inverse=True)
    hq = keys // 4000 - 2000
    hr = keys % 4000 - 2000
    hx = SIDE * np.sqrt(3) * (hq + hr / 2)
    hy = SIDE * 1.5 * hr
    return iy, ix, inverse, len(keys), hx, hy


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--blocks", default=os.path.join(FIG_DATA, "blocks5km"),
                    help="folder with valid_cells.npz and urban_union_<SSP>.npz")
    ap.add_argument("--out", default=HERE)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    valid = np.load(os.path.join(a.blocks, "valid_cells.npz"))["valid"].astype(float)
    iy, ix, inverse, n_hex, hx, hy = hexagon_index(valid)

    def hex_sum(arr):
        return np.bincount(inverse, weights=arr[iy, ix], minlength=n_hex)

    south_limit = ROB.transform_point(0, -57, ccrs.PlateCarree())[1]
    land = (hex_sum(valid) >= 2000) & (hy > south_limit)

    data = {}
    for p in SSPS:
        z = np.load(os.path.join(a.blocks, f"urban_union_{p}.npz"))
        data[p] = {year: hex_sum(z[f"u{year}"].astype(float)) for year in (2020, 2050, 2100)}
    cap = float(np.percentile(np.concatenate([data[p][2100][land & (data[p][2100] > 0)] for p in SSPS]), 98))

    def radius(u):
        return RMAX * np.sqrt(np.minimum(u, cap) / cap)

    fig = plt.figure(figsize=(PANEL_W, FIG_H))
    for k, p in enumerate(SSPS):
        d = data[p]
        ax = fig.add_axes(panel_rect(k, 0.0, 0.40, 1.0, 0.60), projection=ROB)
        ax.set_extent([-168, 180, -57, 79], crs=ccrs.PlateCarree())
        ax.spines["geo"].set_edgecolor("#BDBDBD")
        ax.spines["geo"].set_linewidth(0.5)
        hexes = [RegularPolygon((x, y), 6, radius=SIDE, orientation=0) for x, y in zip(hx[land], hy[land])]
        ax.add_collection(PatchCollection(hexes, facecolor="#EFEFEF", edgecolor="white", linewidth=0.4, zorder=1))

        sel = np.where(land & (d[2100] > 0))[0]
        colors = CMAP(NORM(np.clip(d[2020][sel], 10, 30000)))
        rings = [Circle((hx[i], hy[i]), radius(d[2100][i])) for i in sel]
        ax.add_collection(PatchCollection(rings, facecolor="none", edgecolor=colors, linewidth=RING, zorder=3))
        sel50 = sel[d[2050][sel] > 0]
        colors50 = CMAP(NORM(np.clip(d[2020][sel50], 10, 30000)))
        discs = [Circle((hx[i], hy[i]), radius(d[2050][i])) for i in sel50]
        ax.add_collection(PatchCollection(discs, facecolor=colors50, edgecolor="none", zorder=4))

        x0, y0, _, h0 = panel_rect(k, 0.006, 0.985, 0, 0)
        fig.text(x0, y0, LETTERS[p], fontsize=10, fontweight="bold", va="top")

        cax = fig.add_axes(panel_rect(k, 0.04, 0.47, 0.17, 0.016))
        cb = fig.colorbar(plt.cm.ScalarMappable(norm=NORM, cmap=CMAP), cax=cax, orientation="horizontal",
                          ticks=[10, 100, 1000, 10000])
        cb.ax.set_xticklabels(["10", "100", "1,000", "10,000"])
        cb.outline.set_linewidth(0.4)
        cb.ax.tick_params(length=2, labelsize=8)
        cb.ax.xaxis.set_label_position("top")
        cb.set_label("Urban land in 2020 (km²)", labelpad=2)

        fig.canvas.draw()
        bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
        xlim = ax.get_xlim()
        m_per_in = (xlim[1] - xlim[0]) / bbox.width

        def fig_xy(lon, lat):
            x, y = ROB.transform_point(lon, lat, ccrs.PlateCarree())
            return fig.transFigure.inverted().transform(ax.transData.transform((x, y)))

        (fx0, fy1), (fx1, fy0) = fig_xy(-31, 3), fig_xy(9, -47)
        key_w, key_h = (fx1 - fx0) * PANEL_W, (fy1 - fy0) * FIG_H
        kax = fig.add_axes([fx0, fy0, fx1 - fx0, fy1 - fy0])
        kax.set_axis_off()
        kax.set_xlim(0, key_w * m_per_in)
        kax.set_ylim(0, key_h * m_per_in)
        kax.set_aspect("equal")
        kax.text(0.02 * m_per_in, key_h * m_per_in, "Urban land\nper hexagon\n(km²)", ha="left", va="top", fontsize=8)
        vals = [1000, 8000, round(cap, -3)]
        cx = 0.02 * m_per_in + RMAX
        tx = cx + RMAX + 0.07 * m_per_in
        rows = [(f"{vals[0]:,.0f}", vals[0], None), (f"{vals[1]:,.0f}", vals[1], None),
                (f"≥{vals[2]:,.0f}", vals[2], None), ("2100 / 2050", vals[2], vals[2] * 0.3)]
        for j, (label, u, u_inner) in enumerate(rows):
            yc = (key_h - 0.52 - j * 0.165) * m_per_in
            kax.add_patch(Circle((cx, yc), radius(u), facecolor="none", edgecolor="#6E6E6E", linewidth=RING))
            if u_inner:
                kax.add_patch(Circle((cx, yc), radius(u_inner), facecolor="#6E6E6E", edgecolor="none"))
            last = kax.text(tx, yc, label, ha="left", va="center", fontsize=8)

        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        to_fig = fig.transFigure.inverted()
        y_cb = min(to_fig.transform((0, t.get_window_extent(renderer).y0))[1]
                   for t in cb.ax.get_xticklabels() if t.get_text())
        y_key = to_fig.transform((0, last.get_window_extent(renderer).y0))[1]
        pos = kax.get_position()
        kax.set_position([pos.x0, pos.y0 + (y_cb - y_key), pos.width, pos.height])

        print(p, {year: int(d[year].sum()) for year in (2020, 2050, 2100)}, "km² in the union")
    print("land hexagons", int(land.sum()), "| largest circle at", round(cap), "km²")

    for ext in ("png", "pdf"):
        fig.savefig(out_name(a.out, NAME, ext), dpi=300, facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    main()
