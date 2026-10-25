"""Project the ESA CCI / C3S land-cover map of 2016 or 2020 (300 m, geographic) onto the 1-km
equal-area analysis grid by nearest-neighbor resampling.

Usage: python analysis/harmonize/landcover_baseline.py 2016
Input: the land-cover netCDF file under data/rasters/source/ESA_CCI (or --source).
Output: data/rasters/baseline/ESA_CCI_<year>_lccs_cea1km.tif"""
import argparse
import sys
from pathlib import Path

import rasterio
from rasterio.crs import CRS
from rasterio.enums import Resampling
from rasterio.transform import Affine
from rasterio.vrt import WarpedVRT
from rasterio.windows import Window

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import BASELINE, CELL, HEIGHT, SOURCE, WIDTH, X0, Y0

FILES = {"2016": "C3S-LC-L4-LCCS-Map-300m-P1Y-2016-v2.1.1.nc",
         "2020": "C3S-LC-L4-LCCS-Map-300m-P1Y-2020-v2.1.1.nc"}

ap = argparse.ArgumentParser()
ap.add_argument("year", choices=sorted(FILES))
ap.add_argument("--source", help="land-cover netCDF file")
a = ap.parse_args()
src_file = Path(a.source) if a.source else SOURCE / "ESA_CCI" / FILES[a.year]
out = BASELINE / f"ESA_CCI_{a.year}_lccs_cea1km.tif"
BASELINE.mkdir(parents=True, exist_ok=True)

T = Affine(CELL, 0, X0, 0, -CELL, Y0)
prof = dict(driver="GTiff", width=WIDTH, height=HEIGHT, count=1, dtype="uint8", crs=CRS.from_epsg(6933),
            transform=T, compress="lzw", tiled=True, blockxsize=512, blockysize=512, BIGTIFF="YES")
BAND = 1000
with rasterio.open(f"netcdf:{src_file}:lccs_class") as s:
    # The netCDF variable reports an identity transform; the true one is the global 300 m grid.
    st = Affine(360.0 / s.width, 0, -180.0, 0, -180.0 / s.height, 90.0)
    with WarpedVRT(s, src_crs=CRS.from_epsg(4326), src_transform=st, crs=CRS.from_epsg(6933), transform=T,
                   width=WIDTH, height=HEIGHT, resampling=Resampling.nearest, src_nodata=0, nodata=0) as vrt, \
            rasterio.open(out, "w", **prof) as d:
        for r0 in range(0, HEIGHT, BAND):
            h = min(BAND, HEIGHT - r0)
            d.write(vrt.read(1, window=Window(0, r0, WIDTH, h)), 1, window=Window(0, r0, WIDTH, h))
            print("rows", r0, flush=True)
print("wrote", out)
