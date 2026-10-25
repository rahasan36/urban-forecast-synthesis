"""Boundary length and its cropland and forest contact for each forecast, by zone (country or region).

The zone of a boundary side is the zone of its nonurban cell. Side lengths use the spherical scale
factors of the equal-area grid. For countries, the script also writes the median across forecasts of
boundary length and of the cropland and forest shares of the boundary for the countries in which all
nine forecasts have a boundary (Fig. 3).

Usage: python analysis/measure/zone_contact.py country SSP1 2100
       python analysis/measure/zone_contact.py region SSP1 2100
Outputs (analysis/output/measure): zone_contact_<zone>_<SSP>_<year>.npz with acc[forecast, zone,
(all, cropland, forest)] in km, and for countries country_contact_<SSP>_<year>.csv."""
import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio
from rasterio.windows import Window

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import (CLASS_CODES, COMMON_VALID, COS30, COUNTRY_CODES, COUNTRY_RASTER, EARTH_RADIUS, HEIGHT, LANDCOVER_2016,
                    PRODUCTS, REGION_RASTER, WIDTH, Y0, forecast_path, output_dir, urban_rank)

ap = argparse.ArgumentParser()
ap.add_argument("zone", choices=["country", "region"])
ap.add_argument("ssp", choices=["SSP1", "SSP3", "SSP5"])
ap.add_argument("year", type=int, choices=[2020, 2050, 2100])
a = ap.parse_args()
s = int(a.ssp[3])
look = np.zeros(256, np.uint8)
for i, (k, codes) in enumerate(CLASS_CODES.items()):
    look[codes] = i + 1
CROP = list(CLASS_CODES).index("Cropland") + 1
FOR = list(CLASS_CODES).index("Forest") + 1
zsrc = rasterio.open(COUNTRY_RASTER if a.zone == "country" else REGION_RASTER)
NZ = 65536 if zsrc.dtypes[0] == "uint16" else 256
acc = np.zeros((len(PRODUCTS), NZ, 3), np.float64)
cv = rasterio.open(COMMON_VALID)
lcr = rasterio.open(LANDCOVER_2016)
fh = [rasterio.open(forecast_path(p, s, a.year)) for p in PRODUCTS]
STEP = 500
t0 = time.time()
for r0 in range(0, HEIGHT, STEP):
    h = min(STEP, HEIGHT - r0)
    he = min(h + 1, HEIGHT - r0)
    win = Window(0, r0, WIDTH, he)
    v = cv.read(1, window=win).astype(bool)
    lc = look[lcr.read(1, window=win)]
    z = zsrc.read(1, window=win).astype(np.int64)
    rows = np.arange(r0, r0 + he)
    y = Y0 - (rows + 0.5) * 1000.0
    cl = np.cos(np.arcsin(np.clip(y * COS30 / EARTH_RADIUS, -1, 1)))
    wh = (1000.0 * COS30 / cl) / 1000.0
    wv = (1000.0 * cl / COS30) / 1000.0
    for pi, p_ in enumerate(PRODUCTS):
        U = (fh[pi].read(1, window=win) >= urban_rank(p_)) & v
        for orient in ("between_columns", "between_rows"):
            if orient == "between_columns":
                a_, b_ = U[:h, :-1], U[:h, 1:]
                vv = v[:h, :-1] & v[:h, 1:]
                zn = np.where(a_, z[:h, 1:], z[:h, :-1])
                cn = np.where(a_, lc[:h, 1:], lc[:h, :-1])
                wgt = np.repeat(wh[:h, None], WIDTH - 1, axis=1)
            else:
                a_, b_ = U[:-1, :], U[1:, :]
                vv = v[:-1, :] & v[1:, :]
                zn = np.where(a_, z[1:, :], z[:-1, :])
                cn = np.where(a_, lc[1:, :], lc[:-1, :])
                wgt = np.repeat(wv[:he - 1, None], WIDTH, axis=1)
            e = (a_ ^ b_) & vv
            idx, wq, cq = zn[e], wgt[e], cn[e]
            acc[pi, :, 0] += np.bincount(idx, weights=wq, minlength=NZ)
            acc[pi, :, 1] += np.bincount(idx[cq == CROP], weights=wq[cq == CROP], minlength=NZ)
            acc[pi, :, 2] += np.bincount(idx[cq == FOR], weights=wq[cq == FOR], minlength=NZ)
    print("rows", r0, round(time.time() - t0, 1), flush=True)
out = output_dir("measure")
np.savez_compressed(out / f"zone_contact_{a.zone}_{a.ssp}_{a.year}.npz", acc=acc)

if a.zone == "country":
    cc = json.loads(COUNTRY_CODES.read_text(encoding="utf-8"))
    inv = {v: k for k, v in cc["code"].items()}
    rows = []
    for c, a3 in inv.items():
        t = acc[:, c, 0]
        if (t > 0).sum() < 9:
            continue
        crop = 100 * acc[:, c, 1] / t
        forest = 100 * acc[:, c, 2] / t
        rows.append(dict(a3=a3, name=cc["name"][a3], n_prod=int((t > 0).sum()), median_boundary_km=np.median(t),
                         crop=np.median(crop), forest=np.median(forest), diff=np.median(crop) - np.median(forest),
                         crop_min=crop.min(), crop_max=crop.max()))
    d = pd.DataFrame(rows).sort_values("diff")
    d.to_csv(out / f"country_contact_{a.ssp}_{a.year}.csv", index=False)
    print(len(d), "countries with a boundary in all nine forecasts")
