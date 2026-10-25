"""Divergent urban boundaries and high spatial uncertainty in global urban land forecasts

Author: Rakibul Ahasan
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio
from pyproj import Transformer
from rasterio.windows import Window

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import PRODUCTS, SSPS, YEARS, ensemble_path, output_dir

OUT = output_dir("agreement")
rows, persist = [], []
for s in SSPS:
    for y in YEARS:
        b = np.load(OUT / f"blocks_SSP{s}_{y}.npz")
        h = b["histogram"]
        for code in range(1, 512):
            rows.append(dict(ssp=s, year=y, code=code, n_products=code.bit_count(),
                             products=" + ".join(p for i, p in enumerate(PRODUCTS) if code & (1 << i)), area_km2=int(h[code])))
        if y > 2020:
            a = b["unanimous"]
            p = b["persistent_agreement"]
            v = a >= 20
            persist.append(dict(ssp=s, year=y, future_agreement_km2=int(a.sum()), persistent_agreement_km2=int(p.sum()),
                                global_percent=float(100 * p.sum() / a.sum()), map_minimum_agreement_km2=20,
                                map_block_n=int(v.sum()), map_block_median_percent=float(np.median(100 * p[v] / a[v]))))
pd.DataFrame(rows).to_csv(OUT / "membership_exact.csv", index=False)
pd.DataFrame(persist).to_csv(OUT / "persistent_agreement.csv", index=False)

result = []
tf = Transformer.from_crs(6933, 4326, always_xy=True)
for s in SSPS:
    for y in YEARS:
        with rasterio.open(ensemble_path(s, y)) as ds:
            sums = np.zeros((90, 2), dtype="int64")
            for row in range(0, ds.height, 256):
                n = min(256, ds.height - row)
                a = ds.read(1, window=Window(0, row, ds.width, n))
                counts = np.column_stack((((a > 0) & (a != 65535)).sum(1), (a == 511).sum(1)))
                yy = ds.transform.f + (np.arange(row, row + n) + .5) * ds.transform.e
                _, lat = tf.transform(np.zeros(n), yy)
                band = np.clip(np.floor((lat + 90) / 2).astype(int), 0, 89)
                np.add.at(sums, band, counts)
            for k, (u, i) in enumerate(sums):
                result.append(dict(ssp=s, year=y, latitude_lower=-90 + 2 * k, latitude_upper=-88 + 2 * k,
                                   latitude_mid=-89 + 2 * k, union_km2=int(u), unanimous_km2=int(i),
                                   share_percent=float(100 * i / u) if u else None))
        print("Latitude profile", s, y, flush=True)
lat = pd.DataFrame(result)
lat.to_csv(OUT / "latitude_agreement.csv", index=False)
ens = pd.read_csv(OUT / "ensemble_summary.csv").query("n_products==9")
membership = pd.read_csv(OUT / "membership_exact.csv")
for r in ens.itertuples():
    q = lat[(lat.ssp == r.ssp) & (lat.year == r.year)]
    assert q.union_km2.sum() == r.union_km2
    assert q.unanimous_km2.sum() == r.unanimous_km2
    assert membership[(membership.ssp == r.ssp) & (membership.year == r.year)].area_km2.sum() == r.union_km2
print("Intersections and latitude profiles match the ensemble totals")
