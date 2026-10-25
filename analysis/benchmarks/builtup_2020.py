"""Divergent urban boundaries and high spatial uncertainty in global urban land forecasts

Author: Rakibul Ahasan
"""
import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import rasterio
from rasterio.crs import CRS
from rasterio.enums import Resampling
from rasterio.transform import Affine
from rasterio.vrt import WarpedVRT
from rasterio.windows import Window

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import (CELL, COMMON_VALID, HEIGHT, PRODUCTS, SOURCE, SSPS, WIDTH, X0, Y0, forecast_path, output_dir,
                    urban_rank)

ap = argparse.ArgumentParser()
ap.add_argument("--source", default=str(SOURCE / "GHSL" / "GHS_BUILT_S_E2020_GLOBE_R2023A_54009_100_V1_0.tif"))
a = ap.parse_args()
BAND, NB = 200, 100
T = Affine(CELL, 0, X0, 0, -CELL, Y0)
hist = np.zeros(NB + 1, dtype=np.int64)
ph = {f"{p}|SSP{x}": np.zeros(NB + 1, dtype=np.int64) for p in PRODUCTS for x in SSPS}
t0 = time.time()
cv = rasterio.open(COMMON_VALID)
fh = {f"{p}|SSP{x}": rasterio.open(forecast_path(p, x, 2020)) for p in PRODUCTS for x in SSPS}
with rasterio.open(a.source, overview_level=0) as src:
    with WarpedVRT(src, crs=CRS.from_epsg(6933), transform=T, width=WIDTH, height=HEIGHT,
                   resampling=Resampling.average, src_nodata=65535, nodata=65535) as vrt:
        for r0 in range(0, HEIGHT, BAND):
            h = min(BAND, HEIGHT - r0)
            win = Window(0, r0, WIDTH, h)
            x = vrt.read(1, window=win).astype("float32")
            x[x == 65535] = 0.0
            f = np.clip(x / 10000.0, 0, 1)
            v = cv.read(1, window=win).astype(bool)
            b = np.minimum((f * NB).astype(np.int32), NB)
            hist += np.bincount(b[v], minlength=NB + 1)
            for k, fp in fh.items():
                m = (fp.read(1, window=win) >= urban_rank(k.split("|")[0])) & v
                ph[k] += np.bincount(b[m], minlength=NB + 1)
            print("rows", r0, round(time.time() - t0, 1), flush=True)
(output_dir("benchmarks") / "builtup_2020_histograms.json").write_text(
    json.dumps({"hist": hist.tolist(), "ph": {k: v.tolist() for k, v in ph.items()}}))
print("cells at or above 0.13:", int(hist[13:].sum()))
