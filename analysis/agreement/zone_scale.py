"""Amount and placement sums for zones of 5, 10, 20, 50 and 100 km (SSP1, 2100).

For each zone size, the sums over zones of O (pair intersections), M (the largest pair intersections
the zone amounts allow), E (pair intersections expected under independent placement in the zone) and
S1 (summed forecast area), by region (the region of the zone's center cell) and, for 10, 20 and 50 km
zones, by distance band from the cells all nine forecasts call urban in 2020.

Usage: python analysis/agreement/zone_scale.py [--part 0|1] [--combine]
Output: analysis/output/agreement/zone_scale_SSP1_2100.npz (or its part files)."""
import argparse
import sys
import time
from pathlib import Path

import numpy as np
import rasterio
from rasterio.windows import Window
from scipy.ndimage import distance_transform_edt

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from agreement_stats import strata
from common import (COMMON_VALID, COS30, EARTH_RADIUS, HEIGHT, PRODUCTS, REGION_RASTER, WIDTH, Y0, forecast_path,
                    output_dir, urban_rank)

ap = argparse.ArgumentParser()
ap.add_argument("--part", type=int, choices=[0, 1])
ap.add_argument("--combine", action="store_true")
a = ap.parse_args()
OUT = output_dir("agreement")
if a.combine:
    p = [np.load(OUT / f"zone_scale_SSP1_2100_part{i}.npz") for i in (0, 1)]
    np.savez(OUT / "zone_scale_SSP1_2100.npz", **{k: p[0][k] + p[1][k] for k in p[0].files})
    print("combined")
    raise SystemExit(0)

BAND, HALO, NR = 500, 60, 21
ZS, ZB = (5, 10, 20, 50, 100), (10, 20, 50)
cv = rasterio.open(COMMON_VALID)
reg = rasterio.open(REGION_RASTER)
fh = {y: [rasterio.open(forecast_path(p, 1, y)) for p in PRODUCTS] for y in (2020, 2100)}
acc = {f"{q}_r{z}": np.zeros(NR) for z in ZS for q in "OMES"}
acc.update({f"{q}_b{z}": np.zeros(12) for z in ZB for q in "OMES"})
t0 = time.time()


def masks(y, r0, h):
    win = Window(0, r0, WIDTH, h)
    return [f.read(1, window=win) >= urban_rank(p) for p, f in zip(PRODUCTS, fh[y])]


def add(s, grp, N, tag):
    E = np.where(s["A"] > 0, (s["S1"] ** 2 - s["S2"]) / (2 * np.maximum(s["A"], 1)), 0)
    for q, val in (("O", s["O"]), ("M", s["M"]), ("E", E), ("S", s["S1"])):
        acc[f"{q}_{tag}"] += np.bincount(grp, weights=val, minlength=N)


rows = {None: range(0, HEIGHT, BAND), 0: range(0, 7000, BAND), 1: range(7000, HEIGHT, BAND)}[a.part]
for r0 in rows:
    h = min(BAND, HEIGHT - r0)
    e0, e1 = max(0, r0 - HALO), min(HEIGHT, r0 + h + HALO)
    vE = cv.read(1, window=Window(0, e0, WIDTH, e1 - e0)).astype(bool)
    mE = masks(2020, e0, e1 - e0)
    A9 = np.logical_and.reduce(mE) & vE
    del mE
    yc = Y0 - ((r0 + h / 2) + 0.5) * 1000.0
    cl = np.cos(np.arcsin(np.clip(yc * COS30 / EARTH_RADIUS, -1, 1)))
    d = distance_transform_edt(~A9, sampling=(COS30 / cl, cl / COS30))[r0 - e0:r0 - e0 + h] if A9.any() else np.full((h, WIDTH), 1e9)
    del A9, vE
    b = np.where(d == 0, 0, np.minimum(11, np.ceil(d / 5.0))).astype(np.int64)
    del d
    win = Window(0, r0, WIDTH, h)
    v = cv.read(1, window=win).astype(bool)
    ur = reg.read(1, window=win).astype(np.int64)
    ur[ur >= NR] = 0
    m = masks(2100, r0, h)
    for z in ZS:
        nbx, nby = (WIDTH + z - 1) // z, (h + z - 1) // z
        blk = (np.arange(h)[:, None] // z) * nbx + (np.arange(WIDTH)[None, :] // z)
        cy = np.minimum(np.arange(nby) * z + z // 2, h - 1)
        cx = np.minimum(np.arange(nbx) * z + z // 2, WIDTH - 1)
        zr = ur[cy][:, cx].ravel()
        add(strata(blk, v, m, nby * nbx), zr, NR, f"r{z}")
        if z in ZB:
            add(strata(blk * 12 + b, v, m, nby * nbx * 12), np.arange(nby * nbx * 12) % 12, 12, f"b{z}")
        del blk
    print("rows", r0, round(time.time() - t0, 1), flush=True)
name = "zone_scale_SSP1_2100.npz" if a.part is None else f"zone_scale_SSP1_2100_part{a.part}.npz"
np.savez(OUT / name, **acc)
print("saved", name, flush=True)
