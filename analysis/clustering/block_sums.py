"""Divergent urban boundaries and high spatial uncertainty in global urban land forecasts

Author: Rakibul Ahasan
"""
import argparse
import sys
import time
from pathlib import Path

import numpy as np
import rasterio
from rasterio.windows import Window

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import (CLASS_CODES, COMMON_VALID, COS30, EARTH_RADIUS, HEIGHT, LANDCOVER_2016, PRODUCTS, WIDTH, Y0,
                    forecast_path, output_dir, urban_rank)

B = 50
NBY, NBX = (HEIGHT + B - 1) // B, (WIDTH + B - 1) // B
KEYS = ["union", "unan", "dis", "bt", "bc", "bf", "nsum"]
look = np.zeros(256, np.uint8)
for i, (k, codes) in enumerate(CLASS_CODES.items()):
    look[codes] = i + 1
CROP = list(CLASS_CODES).index("Cropland") + 1
FOR = list(CLASS_CODES).index("Forest") + 1


def blk(a):
    pad = (-a.shape[0]) % B
    if pad:
        a = np.vstack([a, np.zeros((pad, a.shape[1]), a.dtype)])
    padx = (-a.shape[1]) % B
    if padx:
        a = np.hstack([a, np.zeros((a.shape[0], padx), a.dtype)])
    return a.reshape(a.shape[0] // B, B, a.shape[1] // B, B).sum(axis=(1, 3))


def accumulate(ssp, year, row_starts):
    acc = {k: np.zeros((NBY, NBX), np.float64) for k in KEYS}
    cv = rasterio.open(COMMON_VALID)
    lcr = rasterio.open(LANDCOVER_2016)
    fh = [rasterio.open(forecast_path(p, ssp, year)) for p in PRODUCTS]
    t0 = time.time()
    for r0 in row_starts:
        h = min(B * 10, HEIGHT - r0)
        win = Window(0, r0, WIDTH, h)
        v = cv.read(1, window=win).astype(bool)
        lc = look[lcr.read(1, window=win)]
        n = np.zeros((h, WIDTH), np.uint8)
        for p_, f in zip(PRODUCTS, fh):
            n += ((f.read(1, window=win) >= urban_rank(p_)) & v).astype(np.uint8)
        U = (n > 0) & v
        UN = (n == 9) & v
        rows = np.arange(r0, r0 + h)
        y = Y0 - (rows + 0.5) * 1000.0
        cl = np.cos(np.arcsin(np.clip(y * COS30 / EARTH_RADIUS, -1, 1)))
        wh = (1000.0 * COS30 / cl) / 1000.0
        wv = (1000.0 * cl / COS30) / 1000.0
        bt = np.zeros((h, WIDTH))
        bc = np.zeros((h, WIDTH))
        bf = np.zeros((h, WIDTH))
        a, b = U[:, :-1], U[:, 1:]
        e = (a ^ b) & v[:, :-1] & v[:, 1:]
        nb = np.where(a, lc[:, 1:], lc[:, :-1])
        wgt = np.repeat(wh[:, None], WIDTH - 1, axis=1)
        np.add.at(bt, (slice(None), slice(0, WIDTH - 1)), np.where(e, wgt, 0))
        bc[:, :WIDTH - 1] += np.where(e & (nb == CROP), wgt, 0)
        bf[:, :WIDTH - 1] += np.where(e & (nb == FOR), wgt, 0)
        a, b = U[:-1, :], U[1:, :]
        e = (a ^ b) & v[:-1, :] & v[1:, :]
        nb = np.where(a, lc[1:, :], lc[:-1, :])
        wgt = np.repeat(wv[:-1, None], WIDTH, axis=1)
        bt[:h - 1, :] += np.where(e, wgt, 0)
        bc[:h - 1, :] += np.where(e & (nb == CROP), wgt, 0)
        bf[:h - 1, :] += np.where(e & (nb == FOR), wgt, 0)
        b0 = r0 // B
        for k, arr in (("nsum", np.where(v, n, 0).astype(np.float64)), ("union", U.astype(np.float64)),
                       ("unan", UN.astype(np.float64)), ("dis", (U & ~UN).astype(np.float64)),
                       ("bt", bt), ("bc", bc), ("bf", bf)):
            s = blk(arr)
            acc[k][b0:b0 + s.shape[0], :s.shape[1]] += s
        print("rows", r0, round(time.time() - t0, 1), flush=True)
    return acc


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("ssp", choices=["SSP1", "SSP3", "SSP5"])
    ap.add_argument("year", type=int, choices=[2020, 2050, 2100])
    a = ap.parse_args()
    acc = accumulate(int(a.ssp[3]), a.year, range(0, HEIGHT, B * 10))
    np.savez_compressed(output_dir("clustering") / f"block_sums_{a.ssp}_{a.year}.npz", **acc)
    print("saved")
