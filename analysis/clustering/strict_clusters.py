"""Stricter specification of local clustering of urban support (Fig. S12B) and block support.

Block support is the summed urban votes divided by nine times the union area of a 50-km block,
computed from the membership rasters for SSP1, SSP3 and SSP5 in 2020 and 2100 and checked against
the block unions of agreement/ensemble.py and the measured urban areas. The stricter specification
uses support itself (positive z means higher support), blocks with at least 25 km2 of union urban
land, a 9 x 9 block neighborhood including the focal block, at least nine eligible blocks in the
neighborhood, two-sided normal p-values and Benjamini-Yekutieli false discovery rate control at
0.05 (Eq. S13). The script also measures the outer boundary of the SSP1 2100 union footprint
adjoining 2016 cropland, with WGS84 ellipsoidal side lengths, by block.

Run after agreement/ensemble.py, measure/measure_maps.py and block_regions.py.
Outputs (analysis/output/clustering): support_SSP<s>_<year>.npz, union_cropland_boundary_2100.npz,
support_clusters.npz and strict_clusters_settings.json"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio
from pyproj import Geod, Transformer
from rasterio.windows import Window
from scipy.ndimage import convolve
from scipy.stats import false_discovery_control, norm

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import CLASS_CODES, GRID, LANDCOVER_2016, OUTPUT, ensemble_path, output_dir

OUT = output_dir("clustering")
g = json.loads(GRID.read_text())
H, W = g["height"], g["width"]
B = 50
HH, WW = (H + B - 1) // B, (W + B - 1) // B
vote = np.zeros(65536, dtype="uint8")
vote[:512] = [i.bit_count() for i in range(512)]
checks = []


def aggregate(a):
    return np.pad(a, ((0, 0), (0, WW * B - W))).reshape(a.shape[0], WW, B).sum((0, 2))


glob = pd.read_csv(OUTPUT / "measure" / "global_metrics.csv")
for ssp in (1, 3, 5):
    for year in (2020, 2100):
        out = OUT / f"support_SSP{ssp}_{year}.npz"
        if out.exists():
            continue
        u = np.zeros((HH, WW), dtype="int64")
        ns = u.copy()
        with rasterio.open(ensemble_path(ssp, year)) as ds:
            for bi, row in enumerate(range(0, H, B)):
                a = ds.read(1, window=Window(0, row, W, min(B, H - row)))
                u[bi] = aggregate((a > 0) & (a != 65535))
                ns[bi] = aggregate(vote[a])
        reference = np.load(OUTPUT / "agreement" / f"blocks_SSP{ssp}_{year}.npz")
        assert np.array_equal(u, reference["union"])
        total = glob.loc[(glob.ssp == ssp) & (glob.year == year) & (glob.threshold == .4), "urban_area_km2"].sum()
        assert ns.sum() == total
        a = np.divide(ns, 9 * u, out=np.full_like(u, np.nan, dtype=float), where=u > 0)
        assert np.nanmin(a) >= 1 / 9 - 1e-12 and np.nanmax(a) <= 1
        np.savez_compressed(out, union=u, vote_sum=ns, support=a)
        checks.append(dict(ssp=ssp, year=year, union_km2=int(u.sum()), summed_product_area_km2=int(ns.sum())))
        print("support verified", ssp, year, flush=True)

rid = np.load(OUT / "block_regions.npz")["region_valid"]

out = OUT / "union_cropland_boundary_2100.npz"
if not out.exists():
    inv = Transformer.from_crs(6933, 4326, always_xy=True)
    geod = Geod(ellps="WGS84")
    ys = g["transform"][5] + np.arange(H + 1) * g["transform"][4]
    lo, la = inv.transform(np.zeros(H + 1), ys)
    lo1, la1 = inv.transform(np.full(H + 1, 1000.), ys)
    hz = np.asarray(geod.inv(lo, la, lo1, la1)[2]) / 1000
    vt = np.asarray(geod.inv(lo[:-1], la[:-1], lo[1:], la[1:])[2]) / 1000
    codes = CLASS_CODES["Cropland"]
    bc = np.zeros((HH, WW))
    side_count = 0
    with rasterio.open(ensemble_path(1, 2100)) as es, rasterio.open(LANDCOVER_2016) as ls:
        for row in range(0, H, 250):
            end = min(row + 250, H)
            readend = min(end + 1, H)
            win = Window(0, row, W, readend - row)
            a = es.read(1, window=win)
            lc = ls.read(1, window=win)
            v = a != 65535
            u = (a > 0) & v
            c = np.isin(lc, codes)
            n = end - row
            edge = v[:n, :-1] & v[:n, 1:] & ((u[:n, :-1] & ~u[:n, 1:] & c[:n, 1:]) | (~u[:n, :-1] & u[:n, 1:] & c[:n, :-1]))
            rr, cc = np.nonzero(edge)
            urbancc = cc + (~u[rr, cc]).astype(int)
            bc += np.bincount(((rr + row) // B) * WW + urbancc // B, weights=vt[rr + row], minlength=HH * WW).reshape(HH, WW)
            side_count += len(rr)
            n = min(end, H - 1) - row
            edge = v[:n] & v[1:n + 1] & ((u[:n] & ~u[1:n + 1] & c[1:n + 1]) | (~u[:n] & u[1:n + 1] & c[:n]))
            rr, cc = np.nonzero(edge)
            urbanrr = rr + (~u[rr, cc]).astype(int) + row
            bc += np.bincount((urbanrr // B) * WW + cc // B, weights=hz[rr + row + 1], minlength=HH * WW).reshape(HH, WW)
            side_count += len(rr)
    np.savez_compressed(out, cropland_boundary_km=bc)
    print("union cropland boundary", bc.sum(), side_count, flush=True)
else:
    bc = np.load(out)["cropland_boundary_km"]

z = np.load(OUT / "support_SSP1_2100.npz")
a = z["support"]
u = z["union"]
mask = u >= 25
n = int(mask.sum())
x = a[mask]
sd = x.std(ddof=0)
kernel = np.ones((9, 9))
k = convolve(mask.astype(float), kernel, mode="constant", cval=0)
local = convolve(np.where(mask, a, 0.), kernel, mode="constant", cval=0)
den = sd * np.sqrt(np.maximum((n * k - k * k) / (n - 1), 0))
zs = np.divide(local - x.mean() * k, den, out=np.full_like(a, np.nan), where=den > 0)
tested = mask & (k >= 9) & np.isfinite(zs)
p = 2 * norm.sf(np.abs(zs[tested]))
q = false_discovery_control(p, method="by")
sig = np.zeros_like(mask)
sig[tested] = q <= .05
cat = np.zeros_like(u, dtype="int8")
cat[sig & (zs > 0)] = 1
cat[sig & (zs < 0)] = -1
np.savez_compressed(OUT / "support_clusters.npz", support=a, union=u, mask=mask, tested=tested, z=zs, category=cat,
                    neighborhood_n=k, cropland_boundary_km=bc, region=rid)
settings = dict(block_grid_km=50, block_area_km2=2500, crs="EPSG:6933", transform=g["transform"], min_union_km2=25,
                neighborhood="9 x 9 projected blocks, including self", min_neighbors_including_self=9,
                test="two-sided normal approximation", multiple_testing="Benjamini-Yekutieli FDR q <= 0.05",
                support_definition="sum of urban votes / (9 * urban union area)", eligible_blocks=n,
                tested_blocks=int(tested.sum()), higher_support_blocks=int((cat == 1).sum()),
                lower_support_blocks=int((cat == -1).sum()), union_crop_boundary_all_km=float(bc.sum()),
                union_crop_boundary_eligible_km=float(bc[mask].sum()),
                boundary_coverage_percent=float(100 * bc[mask].sum() / bc.sum()))
(OUT / "strict_clusters_settings.json").write_text(json.dumps(settings, indent=2))
if checks:
    (OUT / "support_checks.json").write_text(json.dumps(checks, indent=2))
print(json.dumps(settings, indent=2), flush=True)
