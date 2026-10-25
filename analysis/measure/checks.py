"""Accounting checks across the measurements, and gross and net changes of overlap.

Checks that class contact sums to the boundary length and class overlap to the urban area of each
map, that the boundary segments tallied by the number of forecasts sharing them (agreement/ensemble.py)
add up to the separately measured boundary lengths, and that urban area falls as the threshold of the
fractional forecasts rises. Gross gains and losses of overlap from 2020 to 2100 are tallied from the
ensemble rasters and checked against the net changes.

Run after measure/summaries.py and agreement/ensemble.py. Outputs (analysis/output/measure):
gross_and_net_changes.csv and numerical_checks.json."""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio
from rasterio.windows import Window

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import (CLASS_CODES, CLASSES, COMMON_VALID, FRACTIONAL, LANDCOVER_2016, OUTPUT, PRODUCTS, REGION_RASTER,
                    REGIONS, REGIONS_MAIN, SSPS, YEARS, ensemble_path, output_dir)

OUT = output_dir("measure")
g = pd.read_csv(OUT / "global_metrics.csv")
lc = pd.read_csv(OUT / "landcover_metrics.csv")
cl = pd.read_csv(OUT / "landcover_endpoint_changes.csv")
assert len(g[g.threshold == .4]) == 81
assert len(g) == 117
gg = lc[lc.region == "Global"].groupby(["product", "ssp", "year", "threshold"])
idx = ["product", "ssp", "year", "threshold"]
assert np.allclose(gg.contact_km.sum().sort_index(), g.set_index(idx).perimeter_km.sort_index())
assert np.allclose(gg.contact_share_percent.sum(), 100)
assert np.allclose(gg.overlap_km2.sum().sort_index(), g.set_index(idx).urban_area_km2.sort_index())
bd = pd.read_csv(OUTPUT / "agreement" / "direct_boundary_agreement.csv")
for s in SSPS:
    for year in YEARS:
        for cls in ["All boundaries", "Cropland", "Forest"]:
            a = bd[(bd.ssp == s) & (bd.year == year) & (bd.landcover == cls)]
            direct = (a.n_products * a.segment_length_km).sum()
            if cls == "All boundaries":
                independent = g[(g.ssp == s) & (g.year == year) & (g.threshold == .4)].perimeter_km.sum()
            else:
                independent = lc[(lc.ssp == s) & (lc.year == year) & (lc.threshold == .4) & (lc.region == "Global")
                                 & (lc.landcover == cls)].contact_km.sum()
            assert abs(direct - independent) < .01, (s, year, cls, direct, independent)
for prod in FRACTIONAL:
    for s in SSPS:
        for y in YEARS:
            a = g[(g["product"] == prod) & (g.ssp == s) & (g.year == y)].sort_values("threshold").urban_area_km2.to_numpy()
            assert np.all(np.diff(a) <= 0)
with rasterio.open(COMMON_VALID) as r:
    v = r.read(1).astype(bool)
    H, W = r.shape
with rasterio.open(LANDCOVER_2016) as r:
    lcmap = r.read(1)
with rasterio.open(REGION_RASTER) as r:
    reg = r.read(1)
lookup = np.zeros(256, dtype="uint8")
for c, (_, codes) in enumerate(CLASS_CODES.items()):
    lookup[codes] = c
lookup[190] = 12
lcmap = lookup[lcmap]
rows = []
for s in SSPS:
    accum = np.zeros((9, 2, 21, 13), dtype="int64")
    with rasterio.open(ensemble_path(s, 2020)) as a, rasterio.open(ensemble_path(s, 2100)) as b:
        for row in range(0, H, 250):
            win = Window(0, row, W, min(250, H - row))
            x = a.read(1, window=win)
            y = b.read(1, window=win)
            sl = slice(row, row + int(win.height))
            vv, rr, ll = v[sl], reg[sl], lcmap[sl]
            for k in range(9):
                old = (x & (1 << k)) > 0
                new = (y & (1 << k)) > 0
                for t, sel in enumerate([(~old) & new & vv, old & (~new) & vv]):
                    accum[k, t] += np.bincount(rr[sel].astype("int32") * 13 + ll[sel], minlength=21 * 13).reshape(21, 13)
    for k, prod in enumerate(PRODUCTS):
        for ri in [-1] + REGIONS_MAIN:
            a = accum[k].sum(1) if ri == -1 else accum[k, :, ri, :]
            region = "Global" if ri == -1 else REGIONS[ri]
            for c, name in enumerate(CLASSES):
                net = cl[(cl["product"] == prod) & (cl.ssp == s) & (cl.region == region) & (cl.landcover == name)
                         & (cl.threshold == .4)].overlap_km2_change.iloc[0]
                assert a[0, c] - a[1, c] == net, (prod, s, region, name)
                rows.append(dict(product=prod, ssp=s, region=region, landcover=name, gross_gain_km2=int(a[0, c]),
                                 gross_loss_km2=int(a[1, c]), net_change_km2=int(a[0, c] - a[1, c])))
pd.DataFrame(rows).to_csv(OUT / "gross_and_net_changes.csv", index=False)
result = dict(primary_maps=81, threshold_map_summaries=117, landcover_contact_sums_pass=True,
              area_class_sums_pass=True, threshold_area_monotonicity_pass=True,
              direct_boundary_membership_matches_independent_product_lengths=True,
              gross_minus_loss_equals_net_all_regions_classes=True)
(OUT / "numerical_checks.json").write_text(json.dumps(result, indent=2))
print("CHECKED", result, flush=True)
