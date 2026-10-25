"""Urban area and boundary length of each forecast in 100 x 100 km blocks, for one SSP and year.

Area counts the valid 1-km cells a forecast classifies as urban. Boundary length counts the grid
sides joining an urban and a nonurban cell, both valid, weighted by their ground length on the
equal-area grid (cos30/cos(lat) km for sides between columns, cos(lat)/cos30 km for sides between
rows) and assigned to the block of the upper or left cell. Used for the block dispersion of
Fig. 2 and Figs. S3 and S4.

Usage: python analysis/measure/blocks.py SSP1 2100
Output: analysis/output/measure/blocks100_<SSP>_<year>.npz with area[9, rows, cols] and bnd[9, rows, cols]."""
import argparse
import sys
import time
from pathlib import Path

import numpy as np
import rasterio
from rasterio.windows import Window

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import COMMON_VALID, COS30, EARTH_RADIUS, HEIGHT, PRODUCTS, WIDTH, Y0, forecast_path, output_dir, urban_rank

ap = argparse.ArgumentParser()
ap.add_argument("ssp", choices=["SSP1", "SSP3", "SSP5"])
ap.add_argument("year", type=int, choices=[2020, 2050, 2100])
a = ap.parse_args()
s = int(a.ssp[3])
B = 100
NBY, NBX = (HEIGHT + B - 1) // B, (WIDTH + B - 1) // B
cv = rasterio.open(COMMON_VALID)
fh = [rasterio.open(forecast_path(p, s, a.year)) for p in PRODUCTS]
area = np.zeros((9, NBY, NBX))
bnd = np.zeros((9, NBY, NBX))
t0 = time.time()
BAND = 500
colblk = np.arange(WIDTH) // B
for r0 in range(0, HEIGHT, BAND):
    h = min(BAND, HEIGHT - r0)
    hh = min(h + 1, HEIGHT - r0)
    win = Window(0, r0, WIDTH, hh)
    v = cv.read(1, window=win).astype(bool)
    rows = r0 + np.arange(hh)
    yc = Y0 - (rows + 0.5) * 1000.0
    cl = np.cos(np.arcsin(np.clip(yc * COS30 / EARTH_RADIUS, -1, 1)))
    rh = (COS30 / cl)[:, None]
    cw = (cl / COS30)[:, None]
    rb = rows // B
    for k, (p, f) in enumerate(zip(PRODUCTS, fh)):
        u = (f.read(1, window=win) >= urban_rank(p)) & v
        uh = u[:h]
        vh = v[:h]
        key = rb[:h, None] * NBX + colblk[None, :]
        area[k].ravel()[:] += np.bincount(key[uh], minlength=NBY * NBX)
        sx = (uh[:, 1:] ^ uh[:, :-1]) & vh[:, 1:] & vh[:, :-1]
        bnd[k].ravel()[:] += np.bincount(key[:, :-1][sx], weights=np.broadcast_to(rh[:h], sx.shape)[sx], minlength=NBY * NBX)
        if hh > 1:
            n = hh - 1
            sy = (u[1:hh] ^ u[:n]) & v[1:hh] & v[:n]
            sy = sy[:h] if n >= h else sy
            m = sy.shape[0]
            bnd[k].ravel()[:] += np.bincount(key[:m][sy], weights=np.broadcast_to(cw[1:m + 1], sy.shape)[sy], minlength=NBY * NBX)
    print("rows", r0, round(time.time() - t0, 1), flush=True)
out = output_dir("measure") / f"blocks100_{a.ssp}_{a.year}.npz"
np.savez_compressed(out, area=area.astype(np.float32), bnd=bnd.astype(np.float32))
print("saved", out, "area", area.sum(axis=(1, 2)).round().tolist(), flush=True)
