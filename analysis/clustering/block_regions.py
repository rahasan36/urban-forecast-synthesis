"""Divergent urban boundaries and high spatial uncertainty in global urban land forecasts

Author: Rakibul Ahasan
"""
import sys
from pathlib import Path

import numpy as np
import rasterio
from rasterio.windows import Window

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import COMMON_VALID, HEIGHT, REGION_RASTER, WIDTH, output_dir

B = 50
HH, WW = (HEIGHT + B - 1) // B, (WIDTH + B - 1) // B


def aggregate(a):
    return np.pad(a, ((0, 0), (0, WW * B - WIDTH))).reshape(a.shape[0], WW, B).sum((0, 2))


region_all = np.zeros((HH, WW), dtype="uint8")
region_valid = np.zeros((HH, WW), dtype="uint8")
with rasterio.open(REGION_RASTER) as rs, rasterio.open(COMMON_VALID) as vs:
    for bi, row in enumerate(range(0, HEIGHT, B)):
        w = Window(0, row, WIDTH, min(B, HEIGHT - row))
        a = rs.read(1, window=w)
        v = vs.read(1, window=w) > 0
        for out, sel in ((region_all, np.ones_like(v)), (region_valid, v)):
            counts = np.array([aggregate((a == i) & sel) for i in range(1, 21)])
            out[bi] = np.where(counts.max(axis=0) > 0, counts.argmax(axis=0) + 1, 0)
np.savez_compressed(output_dir("clustering") / "block_regions.npz", region_all=region_all, region_valid=region_valid)
print("saved block regions")
