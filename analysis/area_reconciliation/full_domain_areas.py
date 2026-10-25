"""Urban area of every source map on its full source domain and on the common analysis domain.

Each of the 81 analyzed maps and the base-year map distributed with eight of the forecasts is
reprojected onto a global 1-km EPSG:6933 grid that covers every source extent (origin -17,368,529.34 m,
7,342,708.89 m; 34,738 x 14,686 cells), with nearest-neighbor resampling for categorical maps and
area averaging for the fractional maps (urban at a fraction of 0.40, He 2023 divided by 100). The
script counts urban cells on the full source domain and inside the common valid domain, and for the
fractional maps also sums the unthresholded fractions. The difference between the two counts is the
urban area removed by the restriction to the common domain (Methods).

Inputs: harmonize/source_index.csv and baseline_sources.csv (files under data/rasters/source) and
data/rasters/common_valid.tif. Output: analysis/output/area_reconciliation/maps/*.json"""
import concurrent.futures
import csv
import json
import sys
import time
from pathlib import Path

import numpy as np
import rasterio
from rasterio.enums import Resampling
from rasterio.transform import Affine
from rasterio.vrt import WarpedVRT
from rasterio.windows import Window

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import COMMON_VALID, FRACTIONAL, SOURCE, URBAN_CODE, output_dir

OUT = output_dir("area_reconciliation") / "maps"
OUT.mkdir(exist_ok=True)
TR = Affine(1000, 0, -17368529.34, 0, -1000, 7342708.89)
W, H = 34738, 14686
HERE = Path(__file__).resolve().parent
jobs = []
with open(HERE / "baseline_sources.csv", encoding="utf-8") as fh:
    jobs += [dict(r, kind="baseline", ssp=None, year=int(r["year"])) for r in csv.DictReader(fh)]
with open(HERE.parent / "harmonize" / "source_index.csv", encoding="utf-8") as fh:
    jobs += [dict(r, kind="forecast", ssp=int(r["ssp"]), year=int(r["year"])) for r in csv.DictReader(fh)]


def one(j):
    name = f"{j['product']}_{j['kind']}_{j['ssp']}_{j['year']}"
    save = OUT / (name + ".json")
    if save.exists():
        return json.loads(save.read_text())
    t = time.time()
    p = j["product"]
    frac = p in FRACTIONAL
    full = common = src_valid = 0
    fraction_full = fraction_common = 0.
    with rasterio.Env(GDAL_CACHEMAX=1536 * 1024 ** 2), rasterio.open(SOURCE / j["source_file"]) as src, \
            rasterio.open(COMMON_VALID) as dom:
        nd = 0 if p == "Zhang2023" else src.nodata
        opts = dict(crs="EPSG:6933", transform=TR, width=W, height=H, src_crs=src.crs or "EPSG:4326",
                    resampling=Resampling.average if frac else Resampling.nearest, dtype="float32", nodata=-9999,
                    warp_mem_limit=96)
        if nd is not None:
            opts["src_nodata"] = nd
        with WarpedVRT(src, **opts) as v, WarpedVRT(dom, crs="EPSG:6933", transform=TR, width=W, height=H,
                                                     resampling=Resampling.nearest, src_nodata=0, nodata=0) as dv:
            for row in range(0, H, 128):
                win = Window(0, row, W, min(128, H - row))
                a = v.read(1, window=win)
                valid = (a != -9999) & np.isfinite(a)
                mask = dv.read(1, window=win) == 1
                if frac:
                    if p == "He2023":
                        a /= 100.
                    assert not np.any((a[valid] < -1e-5) | (a[valid] > 1.0001)), p
                    b = (a >= .4) & valid
                    fraction_full += np.sum(a, where=valid, dtype=np.float64)
                    fraction_common += np.sum(a, where=valid & mask, dtype=np.float64)
                else:
                    b = (a == URBAN_CODE[p]) & valid
                full += int(np.count_nonzero(b))
                common += int(np.count_nonzero(b & mask))
                src_valid += int(np.count_nonzero(valid))
    result = dict(product=p, kind=j["kind"], ssp=j["ssp"], year=j["year"], source_file=j["source_file"], crs="EPSG:6933",
                  cell_km2=1, resampling="average" if frac else "nearest",
                  urban_rule="fraction >= 0.4" if frac else f"value = {URBAN_CODE[p]}", raster_full_km2=full,
                  raster_common_km2=common, fraction_full_km2=fraction_full if frac else None,
                  fraction_common_km2=fraction_common if frac else None, source_valid_grid_km2=src_valid,
                  seconds=round(time.time() - t, 1))
    save.write_text(json.dumps(result, indent=2))
    print(name, full, common, result["seconds"], "s", flush=True)
    return result


with concurrent.futures.ThreadPoolExecutor(4) as pool:
    results = list(pool.map(one, jobs))
print("COMPLETE", len(results), flush=True)
