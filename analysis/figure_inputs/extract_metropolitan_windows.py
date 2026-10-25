"""Cut the metropolitan windows drawn in Supplementary Figs. 13 to 15 out of the SSP1 ensemble rasters.

The ensemble rasters store, for each 1-km cell of the EPSG:6933 analysis grid, a 9-bit mask of the
forecasts that classify the cell as urban (65535 marks cells outside the common valid domain). Each
window is a square centered on the given point, 220 km wide for Supplementary Fig. 13 and 160 km wide
for Supplementary Figs. 14 and 15, written as a compressed GeoTIFF to data/figures/metropolitan.
"""
import argparse
import os

import rasterio
from pyproj import Transformer
from rasterio.windows import from_bounds

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))

WINDOWS = {
    220: [("new_york_philadelphia", -74.5, 40.25), ("lagos_ibadan", 3.65, 6.97), ("dhaka", 90.4, 23.85),
          ("sao_paulo", -46.63, -23.55), ("beijing_tianjin", 116.8, 39.4), ("jakarta", 106.85, -6.2)],
    160: [("houston", -95.37, 29.76), ("dallas_fort_worth", -97.0, 32.8), ("london", -0.13, 51.51),
          ("paris", 2.35, 48.86), ("tokyo", 139.76, 35.68), ("shanghai", 121.47, 31.23),
          ("hong_kong", 114.17, 22.32), ("dhaka", 90.4, 23.85), ("lagos", 3.4, 6.5), ("kinshasa", 15.3, -4.33),
          ("dar_es_salaam", 39.27, -6.8), ("addis_ababa", 38.75, 9.03), ("karachi", 67.01, 24.86)],
}
YEARS = [2020, 2100]


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--ensemble_dir", default=os.path.join(REPO, "data", "rasters", "ensemble"),
                    help="folder with ensemble_SSP1_<year>.tif")
    ap.add_argument("--out", default=os.path.join(REPO, "data", "figures", "metropolitan"))
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    to_grid = Transformer.from_crs(4326, 6933, always_xy=True)
    for year in YEARS:
        with rasterio.open(os.path.join(a.ensemble_dir, f"ensemble_SSP1_{year}.tif")) as src:
            for width_km, places in WINDOWS.items():
                half = width_km * 500
                for slug, lon, lat in places:
                    x, y = to_grid.transform(lon, lat)
                    win = from_bounds(x - half, y - half, x + half, y + half, src.transform)
                    win = win.round_offsets().round_lengths()
                    data = src.read(1, window=win)
                    profile = src.profile.copy()
                    profile.update(height=data.shape[0], width=data.shape[1], transform=src.window_transform(win),
                                   compress="deflate", predictor=2, tiled=False)
                    profile.pop("blockxsize", None)
                    profile.pop("blockysize", None)
                    path = os.path.join(a.out, f"{slug}_{width_km}km_SSP1_{year}.tif")
                    with rasterio.open(path, "w", **profile) as dst:
                        dst.write(data, 1)
                    print(path, data.shape)


if __name__ == "__main__":
    main()
