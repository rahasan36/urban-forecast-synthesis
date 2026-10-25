"""Pairwise agreement tallies under SSP1 for 2020 and 2100 (Fig. 4, Table S12 and the amount and
placement split).

hist<Y>[row, col, n]      valid 1-km cells by the number n of forecasts calling them urban, in
                          20-km blocks, for Y = 2020 and 2100
inter<Y>_<unit>[u, pair]  intersection of each of the 36 forecast pairs (itertools order), by country
                          or region
area<Y>_<unit>[u, i]      urban cells of each forecast, by country or region
rel2020[region, n, obs]   valid cells by 2020 forecast count and observed 2020 urban land (class 190)
st_*                      2100 strata = 20-km block x distance band from the cells all nine forecasts
                          call urban in 2020 (bands 0 = inside, 1 to 10 = 5-km steps to 50 km, 11 =
                          beyond), with A, S1, S2, O, M and Un from agreement_stats.strata, only for
                          strata with S1 > 0
Distances are Euclidean on ground spacing, with the row and column spacing of each band of rows.

Usage: python analysis/agreement/aggregate.py            (whole grid)
       python analysis/agreement/aggregate.py --part 0   (rows 0 to 6999; --part 1 for the rest)
       python analysis/agreement/aggregate.py --combine  (joins the two part files)
Output: analysis/output/agreement/aggregate_SSP1.npz (or aggregate_SSP1_part<P>.npz)."""
import argparse
import itertools
import sys
import time
from pathlib import Path

import numpy as np
import rasterio
from rasterio.windows import Window
from scipy.ndimage import distance_transform_edt

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from agreement_stats import strata
from common import (COMMON_VALID, COS30, COUNTRY_RASTER, EARTH_RADIUS, HEIGHT, LANDCOVER_2020, PRODUCTS, REGION_RASTER,
                    URBAN_LANDCOVER_CODE, WIDTH, Y0, forecast_path, output_dir, urban_rank)

ap = argparse.ArgumentParser()
ap.add_argument("--part", type=int, choices=[0, 1])
ap.add_argument("--combine", action="store_true")
a = ap.parse_args()
OUT = output_dir("agreement")
if a.combine:
    p = [np.load(OUT / f"aggregate_SSP1_part{i}.npz") for i in (0, 1)]
    joined = {k: (np.concatenate([p[0][k], p[1][k]]) if k.startswith("st_") else p[0][k] + p[1][k]) for k in p[0].files}
    np.savez_compressed(OUT / "aggregate_SSP1.npz", **joined)
    print("combined")
    raise SystemExit(0)

BAND, HALO, B = 500, 60, 20
NBY, NBX = (HEIGHT + B - 1) // B, (WIDTH + B - 1) // B
PAIRS = list(itertools.combinations(range(9), 2))
NC, NR = 260, 21
cv = rasterio.open(COMMON_VALID)
cty = rasterio.open(COUNTRY_RASTER)
reg = rasterio.open(REGION_RASTER)
lc = rasterio.open(LANDCOVER_2020)
fh = {y: [rasterio.open(forecast_path(p, 1, y)) for p in PRODUCTS] for y in (2020, 2100)}
acc = {f"hist{y}": np.zeros((NBY, NBX, 10), np.int32) for y in (2020, 2100)}
for y in (2020, 2100):
    for u, N in (("country", NC), ("region", NR)):
        acc[f"inter{y}_{u}"] = np.zeros((N, 36), np.int64)
        acc[f"area{y}_{u}"] = np.zeros((N, 9), np.int64)
acc["rel2020"] = np.zeros((NR, 10, 2), np.int64)
ST = {k: [] for k in ("by", "bx", "b", "A", "S1", "S2", "O", "M", "Un")}
t0 = time.time()


def masks(y, r0, h):
    win = Window(0, r0, WIDTH, h)
    return [f.read(1, window=win) >= urban_rank(p) for p, f in zip(PRODUCTS, fh[y])]


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
    uc = cty.read(1, window=win).astype(np.int64)
    uc[uc >= NC] = 0
    ur = reg.read(1, window=win).astype(np.int64)
    ur[ur >= NR] = 0
    nby = (h + B - 1) // B
    blk = (np.arange(h)[:, None] // B) * NBX + (np.arange(WIDTH)[None, :] // B)
    for y in (2020, 2100):
        m = masks(y, r0, h)
        n = np.zeros((h, WIDTH), np.int64)
        for x in m:
            n += x
        hb = np.bincount((blk * 10 + n)[v], minlength=nby * NBX * 10).reshape(nby, NBX, 10)
        acc[f"hist{y}"][r0 // B:r0 // B + nby] += hb.astype(np.int32)
        for u, uu, N in (("country", uc, NC), ("region", ur, NR)):
            for i, x in enumerate(m):
                acc[f"area{y}_{u}"][:, i] += np.bincount(uu[x & v], minlength=N)
            for k, (i, j) in enumerate(PAIRS):
                acc[f"inter{y}_{u}"][:, k] += np.bincount(uu[m[i] & m[j] & v], minlength=N)
        if y == 2020:
            o = (lc.read(1, window=win) == URBAN_LANDCOVER_CODE).astype(np.int64)
            acc["rel2020"] += np.bincount(((ur * 10 + n) * 2 + o)[v], minlength=NR * 20).reshape(NR, 10, 2)
        else:
            K = nby * NBX * 12
            s = strata(blk * 12 + b, v, m, K)
            nz = np.where(s["S1"] > 0)[0]
            ST["by"].append(r0 // B + (nz // 12) // NBX)
            ST["bx"].append((nz // 12) % NBX)
            ST["b"].append(nz % 12)
            for k in ("A", "S1", "S2", "O", "M", "Un"):
                ST[k].append(s[k][nz])
        del m, n
    print("rows", r0, round(time.time() - t0, 1), flush=True)
for k in ST:
    acc["st_" + k] = np.concatenate(ST[k])
name = "aggregate_SSP1.npz" if a.part is None else f"aggregate_SSP1_part{a.part}.npz"
np.savez(OUT / name, **acc)
print("saved", name, round(time.time() - t0, 1), flush=True)
