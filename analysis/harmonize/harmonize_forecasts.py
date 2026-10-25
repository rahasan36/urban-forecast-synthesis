"""Reproject the 81 source forecast maps (nine forecasts, SSP1, SSP3 and SSP5, 2020, 2050 and 2100)
onto the 1-km equal-area grid (EPSG:6933) and build the common valid domain.

Categorical forecasts are written as 0 (nonurban), 1 (urban) and 255 (missing), with nearest-neighbor
resampling. The two fractional forecasts (Gao 2021 and He 2023, He 2023 percentages divided by 100)
are area-averaged and stored as threshold ranks 0 to 3, the number of the thresholds 0.30, 0.40 and
0.50 the urban fraction reaches, with 255 for missing. The common valid domain holds the cells that
are valid in all 81 maps and in the 2016 land-cover baseline, between 55.908°S and 75.004°N.

Inputs: the source files listed in source_index.csv under data/rasters/source, the region raster
(its grid defines the analysis grid) and the 2016 land-cover baseline.
Outputs: data/rasters/harmonized/<forecast>_SSP<s>_<year>.tif with a JSON histogram beside each,
data/rasters/common_valid.tif and data/rasters/grid.json."""
import concurrent.futures
import csv
import json
import sys
import time
from pathlib import Path

import numpy as np
import rasterio
from pyproj import Transformer
from rasterio.enums import Resampling
from rasterio.vrt import WarpedVRT
from rasterio.windows import Window

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import (COMMON_VALID, FRACTIONAL, GRID, HARMONIZED, LANDCOVER_2016, PRODUCTS,
                    REGION_RASTER, SOURCE, URBAN_CODE, forecast_path, output_dir)

HARMONIZED.mkdir(parents=True, exist_ok=True)
with rasterio.open(REGION_RASTER) as r:
    PROFILE = r.profile
    H, W = r.shape
    TR = r.transform
PROFILE.update(driver="GTiff", count=1, dtype="uint8", nodata=255, compress="lzw", tiled=True,
               blockxsize=512, blockysize=256)
with open(Path(__file__).with_name("source_index.csv"), encoding="utf-8") as fh:
    records = [dict(r, ssp=int(r["ssp"]), year=int(r["year"])) for r in csv.DictReader(fh)]
assert len(records) == 81


def harmonize(d):
    out = forecast_path(d["product"], d["ssp"], d["year"])
    meta = out.with_suffix(".json")
    if meta.exists():
        return json.loads(meta.read_text())
    t = time.time()
    hist = np.zeros(256, dtype="int64")
    p = d["product"]
    src_file = SOURCE / d["source_file"]
    with rasterio.Env(GDAL_CACHEMAX=128 * 1024 ** 2):
        with rasterio.open(src_file) as src:
            frac = p in FRACTIONAL
            nd = src.nodata
            if p == "Zhang2023":
                nd = 0
            opts = dict(crs="EPSG:6933", transform=TR, width=W, height=H, src_crs=src.crs or "EPSG:4326",
                        resampling=Resampling.average if frac else Resampling.nearest, dtype="float32",
                        nodata=-9999, warp_mem_limit=128)
            if nd is not None:
                opts["src_nodata"] = nd
            with WarpedVRT(src, **opts) as v, rasterio.open(out, "w", **PROFILE) as dst:
                for row in range(0, H, 256):
                    win = Window(0, row, W, min(256, H - row))
                    a = v.read(1, window=win)
                    valid = (a != -9999) & np.isfinite(a)
                    if frac:
                        if p == "He2023":
                            a /= 100
                        assert not np.any((a[valid] < -1e-5) | (a[valid] > 1.0001)), (p, float(a[valid].min()), float(a[valid].max()))
                        b = ((a >= .3).astype("uint8") + (a >= .4) + (a >= .5)).astype("uint8")
                    else:
                        b = (a == URBAN_CODE[p]).astype("uint8")
                    b[~valid] = 255
                    dst.write(b, 1, window=win)
                    hist += np.bincount(b.ravel(), minlength=256)
    result = dict(product=p, ssp=d["ssp"], year=d["year"], source=d["source_file"], histogram=hist.tolist(),
                  seconds=round(time.time() - t, 2))
    meta.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print("Harmonized", out.name, result["seconds"], "s", flush=True)
    return result


with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
    results = list(pool.map(harmonize, records))
(output_dir("harmonize") / "harmonization.json").write_text(json.dumps(results, indent=2), encoding="utf-8")

# The common domain spans all forecasts, scenarios and years; missing data never counts as nonurban.
srcs = [rasterio.open(forecast_path(d["product"], d["ssp"], d["year"])) for d in records]
inv = Transformer.from_crs(6933, 4326, always_xy=True)
lats = np.array(inv.transform(np.zeros(H), TR.f + (np.arange(H) + .5) * TR.e)[1])
counts = np.zeros(256, dtype="int64")
with rasterio.open(LANDCOVER_2016) as lc, rasterio.open(COMMON_VALID, "w", **dict(PROFILE, nodata=None)) as dst:
    for row in range(0, H, 256):
        win = Window(0, row, W, min(256, H - row))
        a = lc.read(1, window=win)
        band = lats[row:row + int(win.height)]
        v = (a > 0) & ((band >= -55.908211257) & (band <= 75.00416667))[:, None]
        for r in srcs:
            v &= r.read(1, window=win) != 255
        counts += np.bincount(a[v], minlength=256)
        dst.write(v.astype("uint8"), 1, window=win)
for r in srcs:
    r.close()
grid = dict(width=W, height=H, transform=list(TR)[:6], crs="EPSG:6933", nominal_cell_km2=1,
            valid_cells=int(counts.sum()), baseline_code_counts=counts.tolist(),
            max_source_latitude=75.00416667, min_source_latitude=-55.908211257, products=PRODUCTS,
            fraction_thresholds=[.3, .4, .5], primary_threshold=.4)
GRID.write_text(json.dumps(grid, indent=2), encoding="utf-8")
print("COMMON VALID", int(counts.sum()), "km2", flush=True)
