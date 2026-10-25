"""Divergent urban boundaries and high spatial uncertainty in global urban land forecasts

Author: Rakibul Ahasan
"""
import json
import sys
import time
from pathlib import Path

import rasterio
from rasterio.windows import Window

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import (COMMON_VALID, HEIGHT, LANDCOVER_2020, PRODUCTS, SSPS, URBAN_LANDCOVER_CODE, WIDTH, forecast_path,
                    output_dir, urban_rank)

BAND = 1000
s = {"obs": 0, "valid": 0, "c": {f"{p}|SSP{x}": {"area": 0, "inter": 0} for p in PRODUCTS for x in SSPS}}
t0 = time.time()
lc = rasterio.open(LANDCOVER_2020)
cv = rasterio.open(COMMON_VALID)
fh = {f"{p}|SSP{x}": rasterio.open(forecast_path(p, x, 2020)) for p in PRODUCTS for x in SSPS}
for r0 in range(0, HEIGHT, BAND):
    h = min(BAND, HEIGHT - r0)
    win = Window(0, r0, WIDTH, h)
    v = cv.read(1, window=win).astype(bool)
    o = (lc.read(1, window=win) == URBAN_LANDCOVER_CODE) & v
    s["obs"] += int(o.sum())
    s["valid"] += int(v.sum())
    for k, fp in fh.items():
        a = (fp.read(1, window=win) >= urban_rank(k.split("|")[0])) & v
        s["c"][k]["area"] += int(a.sum())
        s["c"][k]["inter"] += int((a & o).sum())
    print("rows", r0, round(time.time() - t0, 1), flush=True)
(output_dir("benchmarks") / "observed_2020_counts.json").write_text(json.dumps(s))
print("observed urban land", s["obs"], "km2")
