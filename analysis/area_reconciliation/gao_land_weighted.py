"""Source-native urban area of Gao 2021: the urban fraction of each source cell multiplied by the
authors' land-area layer (land_area_km2.tif) on the native grid, for the 2000 base year and the nine
analyzed maps. A diagnostic beside the thresholded footprint, not an author-tabulated total.

Input: the land-area layer under data/rasters/source/Gao2021 (or --land) and the Gao 2021 maps listed
in baseline_sources.csv and harmonize/source_index.csv.
Output: analysis/output/area_reconciliation/gao_land_weighted.json"""
import argparse
import concurrent.futures
import csv
import json
import sys
from pathlib import Path

import numpy as np
import rasterio
from rasterio.windows import Window

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import SOURCE, output_dir

ap = argparse.ArgumentParser()
ap.add_argument("--land", default=str(SOURCE / "Gao2021" / "ssp-1-km-downscaled-urban-land-extent-projection-base-year-ssp-2000-2100-ancillary-layers-geotiff" / "ancillary-layers-geotiff" / "land_area_km2.tif"))
a = ap.parse_args()
HERE = Path(__file__).resolve().parent
jobs = []
with open(HERE / "baseline_sources.csv", encoding="utf-8") as fh:
    jobs += [dict(r, kind="baseline", ssp=None) for r in csv.DictReader(fh) if r["product"] == "Gao2021"]
with open(HERE.parent / "harmonize" / "source_index.csv", encoding="utf-8") as fh:
    jobs += [dict(r, kind="forecast") for r in csv.DictReader(fh) if r["product"] == "Gao2021"]


def run(j):
    total = 0.
    with rasterio.open(SOURCE / j["source_file"]) as src, rasterio.open(a.land) as land:
        assert src.crs == land.crs
        offset = ~land.transform * src.transform
        coloff, rowoff = round(offset.c), round(offset.f)
        assert abs(offset.a - 1) < 1e-8 and abs(offset.e - 1) < 1e-8 and abs(offset.c - coloff) < 1e-5 and abs(offset.f - rowoff) < 1e-5, offset
        for row in range(0, src.height, 128):
            win = Window(0, row, src.width, min(128, src.height - row))
            x = src.read(1, window=win, masked=True)
            lnd = land.read(1, window=Window(coloff, row + rowoff, win.width, win.height), masked=True)
            valid = (~np.ma.getmaskarray(x)) & (~np.ma.getmaskarray(lnd)) & (x.data >= 0) & (x.data <= 1) & (lnd.data > 0)
            total += np.sum(x.data.astype("float64") * lnd.data, where=valid, dtype="float64")
    result = dict(product="Gao2021", kind=j["kind"], ssp=j["ssp"], year=int(j["year"]), derived_native_land_weighted_km2=total)
    print(result, flush=True)
    return result


with concurrent.futures.ThreadPoolExecutor(2) as pool:
    results = list(pool.map(run, jobs))
(output_dir("area_reconciliation") / "gao_land_weighted.json").write_text(json.dumps(results, indent=2))
