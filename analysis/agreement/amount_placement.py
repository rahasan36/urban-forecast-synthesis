"""Split of the shortfall of pooled pairwise agreement into amount and placement (SSP1, 2100).

Within each stratum the largest intersection a forecast pair can reach is min(u_i, u_j), which gives
the ceiling J_max of the pooled Jaccard index, and independent placement gives J_ind. The shortfall
1 - J splits into an amount part (1 - J_max: the forecasts assign different amounts of urban land to
the stratum) and a placement part (J_max - J: the same amount in different cells), Eqs. S8 to S10.
Pooled indices are X / (8 S1 - X) for X = O, M or E summed over the strata.

Results: by distance band from the 2020 unanimous cells (20-km zone x band strata), by zone size
(5 to 100 km) and by region (20-km zone x band strata, each zone assigned to the region of its
center). Run after aggregate.py and zone_scale.py.
Output: analysis/output/agreement/amount_placement_SSP1_2100.json"""
import json
import sys
from pathlib import Path

import numpy as np
import rasterio
from rasterio.enums import Resampling

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import REGION_RASTER, REGIONS, REGIONS_MAIN, output_dir

OUT = output_dir("agreement")
z = np.load(OUT / "aggregate_SSP1.npz")
s = {k[3:]: z[k].astype(float) for k in z.files if k.startswith("st_")}
s["E"] = np.where(s["A"] > 0, (s["S1"] ** 2 - s["S2"]) / (2 * np.maximum(s["A"], 1)), 0)


def split(sel):
    S = s["S1"][sel].sum()

    def f(X):
        return X / (8 * S - X)
    return f(s["O"][sel].sum()), f(s["M"][sel].sum()), f(s["E"][sel].sum())


LOG = {}
b = s["b"].astype(int)
X = np.arange(12)
SB = np.array([split(b == k) for k in X])
LOG["band_J_Jmax_Jind"] = np.round(SB, 3).tolist()
un = np.array([s["Un"][b == k].sum() for k in X])
LOG["union_by_band_km2"] = un.astype(int).tolist()
LOG["share_beyond_10km"] = round(float(un[3:].sum() / un.sum()), 3)

a = np.load(OUT / "zone_scale_SSP1_2100.npz")
ZS = (5, 10, 20, 50, 100)
shares = {}
for zz in ZS:
    S = a[f"S_r{zz}"].sum()

    def f(X):
        return X / (8 * S - X)
    o, m = f(a[f"O_r{zz}"].sum()), f(a[f"M_r{zz}"].sum())
    shares[str(zz)] = [round(float((1 - m) / (1 - o)), 3), round(float((m - o) / (1 - o)), 3)]
LOG["zone_amount_placement_share"] = shares

with rasterio.open(REGION_RASTER) as d:
    RG = d.read(1, out_shape=(d.height // 20, d.width // 20), resampling=Resampling.nearest)
rg = RG[s["by"].astype(int), s["bx"].astype(int)]
SR = {r: split(rg == r) for r in REGIONS_MAIN if (rg == r).any()}
LOG["regions_J_Jmax_Jind"] = {REGIONS[r]: np.round(SR[r], 3).tolist() for r in sorted(SR, key=lambda r: -SR[r][0])}
(OUT / "amount_placement_SSP1_2100.json").write_text(json.dumps(LOG, indent=1))
print(json.dumps({k: LOG[k] for k in ("share_beyond_10km", "zone_amount_placement_share")}))
