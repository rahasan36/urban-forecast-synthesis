"""Divergent urban boundaries and high spatial uncertainty in global urban land forecasts

Author: Rakibul Ahasan
"""
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
from common import BASELINE, CELL, COUNTRY_CODES, COUNTRY_RASTER, HEIGHT, REPO, WIDTH, X0, Y0

ap = argparse.ArgumentParser()
ap.add_argument("--shapefile", default=str(REPO / "data" / "figures" / "naturalearth" / "ne_50m_admin_0_countries.shp"))
a = ap.parse_args()
BASELINE.mkdir(parents=True, exist_ok=True)

meta, idx, geom, fields = read(a.shapefile)
f = list(meta["fields"])
NAME = fields[f.index("NAME")]
A3 = fields[f.index("ADM0_A3")]
tr = Transformer.from_crs("EPSG:4326", "EPSG:6933", always_xy=True)
code = {}
shapes = []
for nm, a3, wkb in zip(NAME, A3, geom):
    if a3 == "ATA":
        continue
    if a3 not in code:
        code[a3] = len(code) + 1
    g = shapely.ops.transform(lambda x, y, z=None: tr.transform(x, y), shapely.from_wkb(wkb))
    shapes.append((g, code[a3]))
name = {a3: str(n) for n, a3 in zip(NAME, A3)}
T = Affine(CELL, 0, X0, 0, -CELL, Y0)
prof = dict(driver="GTiff", width=WIDTH, height=HEIGHT, count=1, dtype="uint16", crs=CRS.from_epsg(6933),
            transform=T, compress="lzw", tiled=True, blockxsize=512, blockysize=512, BIGTIFF="YES")
with rasterio.open(COUNTRY_RASTER, "w", **prof) as d:
    B = 1000
    for r0 in range(0, HEIGHT, B):
        h = min(B, HEIGHT - r0)
        wt = Affine(T[0], 0, T[2], 0, T[4], T[5] + r0 * T[4])
        arr = rasterize(shapes, out_shape=(h, WIDTH), transform=wt, fill=0, dtype="uint16", all_touched=False)
        d.write(arr, 1, window=rasterio.windows.Window(0, r0, WIDTH, h))
COUNTRY_CODES.write_text(json.dumps({"code": code, "name": name}, indent=1), encoding="utf-8")
print("wrote", COUNTRY_RASTER, len(code), "countries")
