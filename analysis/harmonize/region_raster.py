"""Rasterize the 20 reporting regions of the IPCC AR6 country grouping onto the analysis grid.

Region codes follow the alphabetical order of the region names (1 Asia-Pacific Developed to
20 Western Africa). Input: Countries_AR6_20Regions.shp under data/rasters/source/AR6_regions
(or --shapefile). Outputs: data/rasters/baseline/AR6_20regions_cea1km.tif and its code table."""
import argparse
import json
import sys
from pathlib import Path

import rasterio
import shapely
import shapely.ops
from pyogrio.raw import read
from pyproj import Transformer
from rasterio.crs import CRS
from rasterio.features import rasterize
from rasterio.transform import Affine

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import BASELINE, CELL, HEIGHT, REGION_CODES, REGION_RASTER, SOURCE, WIDTH, X0, Y0

ap = argparse.ArgumentParser()
ap.add_argument("--shapefile", default=str(SOURCE / "AR6_regions" / "Countries_AR6_20Regions.shp"))
a = ap.parse_args()
BASELINE.mkdir(parents=True, exist_ok=True)

meta, idx, geom, fields = read(a.shapefile)
names = list(fields[0])
code = {n: i + 1 for i, n in enumerate(sorted(names))}
tr = Transformer.from_crs("EPSG:4326", "EPSG:6933", always_xy=True)
shapes = [(shapely.ops.transform(lambda x, y, z=None: tr.transform(x, y), shapely.from_wkb(wkb)), code[nm])
          for nm, wkb in zip(names, geom)]
T = Affine(CELL, 0, X0, 0, -CELL, Y0)
prof = dict(driver="GTiff", width=WIDTH, height=HEIGHT, count=1, dtype="uint8", crs=CRS.from_epsg(6933),
            transform=T, compress="lzw", tiled=True, blockxsize=512, blockysize=512, BIGTIFF="YES")
with rasterio.open(REGION_RASTER, "w", **prof) as d:
    B = 1000
    for r0 in range(0, HEIGHT, B):
        h = min(B, HEIGHT - r0)
        wt = Affine(T[0], 0, T[2], 0, T[4], T[5] + r0 * T[4])
        arr = rasterize(shapes, out_shape=(h, WIDTH), transform=wt, fill=0, dtype="uint8", all_touched=False)
        d.write(arr, 1, window=rasterio.windows.Window(0, r0, WIDTH, h))
REGION_CODES.write_text(json.dumps({"code": code, "order": sorted(names)}, indent=1), encoding="utf-8")
print("wrote", REGION_RASTER)
