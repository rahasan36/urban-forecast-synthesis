"""Nine-forecast membership, pairwise agreement, temporal transitions and shared boundary segments.

For each scenario and year the script writes the membership raster, in which bit k (Li 2021 first,
FU3D 2025 last) marks the cells forecast k classifies as urban, 65535 marks invalid cells and 511
unanimous cells. From its histogram it derives the union, unanimous, single-forecast and majority
areas, the pair intersections and unions of the 36 forecast pairs with the quantity and allocation
components of their disagreement (Pontius and Millones, 2011), and the leave-one-out ensembles.
Cells are classed as no urban votes, disputed (one to eight votes) or unanimous; transitions from
2020 separate inherited disagreement (disputed in 2020 and later) from new disagreement. Boundary
segments are tallied by the number of forecasts whose boundary includes them, for all boundaries and
for sides adjoining 2016 cropland and forest, with WGS84 ellipsoidal side lengths.

Outputs: data/rasters/harmonized/ensemble_SSP<s>_<year>.tif and, in analysis/output/agreement,
ensemble_summary.csv, pairwise_agreement.csv, leave_one_out_agreement.csv, temporal_transitions.csv,
direct_boundary_agreement.csv and blocks_SSP<s>_<year>.npz (50-km block sums of the union,
unanimous, majority, single-forecast, valid, inherited and new disagreement and persistent agreement
cells, with the membership histogram)."""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio
from pyproj import Geod, Transformer
from rasterio.windows import Window

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import (CLASS_CODES, COMMON_VALID, GRID, LANDCOVER_2016, PRODUCTS, SSPS, YEARS, ensemble_path,
                    forecast_path, output_dir, urban_mask)

OUT = output_dir("agreement")
grid = json.loads(GRID.read_text())
H, W = grid["height"], grid["width"]
COUNT = np.array([i.bit_count() for i in range(512)], dtype="uint8")
with rasterio.open(COMMON_VALID) as r:
    VALID = r.read(1).astype(bool)
    profile = r.profile
    TR = r.transform
with rasterio.open(LANDCOVER_2016) as r:
    LC = r.read(1)
profile.update(dtype="uint16", nodata=65535, compress="lzw")


def blocks(a, size=50):
    hh = (a.shape[0] + size - 1) // size
    ww = (a.shape[1] + size - 1) // size
    b = np.pad(a, ((0, hh * size - a.shape[0]), (0, ww * size - a.shape[1])))
    return b.reshape(hh, size, ww, size).sum((1, 3))


summary, pairs, loo, transition, edge_records = [], [], [], [], []
inv = Transformer.from_crs(6933, 4326, always_xy=True)
geod = Geod(ellps="WGS84")
y = TR.f + np.arange(H + 1) * TR.e
lo, la = inv.transform(np.zeros(H + 1), y)
lo1, la1 = inv.transform(np.ones(H + 1) * 1000, y)
hz = np.asarray(geod.inv(lo, la, lo1, la1)[2]) / 1000
vt = np.asarray(geod.inv(lo[:-1], la[:-1], lo[1:], la[1:])[2]) / 1000
for s in SSPS:
    previous = None
    for year in YEARS:
        p = ensemble_path(s, year)
        if not p.exists():
            srcs = [rasterio.open(forecast_path(prod, s, year)) for prod in PRODUCTS]
            with rasterio.open(p, "w", **profile) as dst:
                for row in range(0, H, 250):
                    win = Window(0, row, W, min(250, H - row))
                    a = np.zeros((int(win.height), W), dtype="uint16")
                    for k, (prod, r) in enumerate(zip(PRODUCTS, srcs)):
                        a |= urban_mask(r.read(1, window=win), prod).astype("uint16") << k
                    a[~VALID[row:row + int(win.height)]] = 65535
                    dst.write(a, 1, window=win)
            for r in srcs:
                r.close()
        with rasterio.open(p) as r:
            a = r.read(1)
        hist = np.bincount(a[VALID], minlength=512)
        a[~VALID] = 0
        cnt = COUNT[a]
        union = int(hist[1:].sum())
        inter = int(hist[511])
        for subset, title in [(511, "Expanded nine"), (15, "Original four on common domain")]:
            code = np.arange(512) & subset
            pop = np.array([int(k).bit_count() for k in code])
            n = subset.bit_count()
            u = int(hist[pop > 0].sum())
            it = int(hist[pop == n].sum())
            summary.append(dict(ssp=s, year=year, ensemble=title, n_products=n, union_km2=u, unanimous_km2=it,
                                unanimous_share_percent=100 * it / u, single_km2=int(hist[pop == 1].sum()),
                                majority_km2=int(hist[pop >= n // 2 + 1].sum())))
        for k, prod in enumerate(PRODUCTS):
            code = np.arange(512) & (511 ^ (1 << k))
            n = np.array([int(i).bit_count() for i in code])
            u = int(hist[n > 0].sum())
            it = int(hist[n == 8].sum())
            loo.append(dict(ssp=s, year=year, omitted=prod, union_km2=u, unanimous_km2=it, unanimous_share_percent=100 * it / u))
            for j in range(k + 1, 9):
                x = (np.arange(512) & (1 << k)) > 0
                z = (np.arange(512) & (1 << j)) > 0
                A = int(hist[x].sum())
                B = int(hist[z].sum())
                I = int(hist[x & z].sum())
                U = A + B - I
                D = A + B - 2 * I
                pairs.append(dict(ssp=s, year=year, product_a=prod, product_b=PRODUCTS[j], area_a_km2=A, area_b_km2=B,
                                  intersection_km2=I, union_km2=U, jaccard=I / U, disagreement_km2=D,
                                  quantity_disagreement_km2=abs(A - B), allocation_disagreement_km2=D - abs(A - B),
                                  quantity_share_percent=100 * abs(A - B) / D if D else 0))
        state = np.where(cnt == 0, 0, np.where(cnt == 9, 2, 1)).astype("uint8")
        maps = dict(union=blocks(cnt > 0), unanimous=blocks(cnt == 9), majority=blocks(cnt >= 5), single=blocks(cnt == 1),
                    valid=blocks(VALID), histogram=hist)
        if year == 2020:
            previous = state.copy()
        else:
            trans = np.bincount(previous[VALID] * 3 + state[VALID], minlength=9).reshape(3, 3)
            for i in range(3):
                for j in range(3):
                    transition.append(dict(ssp=s, year=year, state_2020=i, state_later=j, area_km2=int(trans[i, j])))
            maps.update(inherited_disagreement=blocks((previous == 1) & (state == 1)),
                        new_disagreement=blocks((previous != 1) & (state == 1)),
                        persistent_agreement=blocks((previous == 2) & (state == 2)))
        # Boundary segments are tallied directly, not inferred from agreement on urban cells.
        eh = {c: np.zeros(10) for c in ["All boundaries", "Cropland", "Forest"]}
        for row in range(0, H, 250):
            end = min(row + 250, H)
            for orient in ["v", "h"]:
                if orient == "v":
                    left, right = a[row:end, :-1], a[row:end, 1:]
                    valid = VALID[row:end, :-1] & VALID[row:end, 1:]
                    lc1, lc2 = LC[row:end, :-1], LC[row:end, 1:]
                    weights = vt[row:end]
                else:
                    en = min(end, H - 1)
                    left, right = a[row:en], a[row + 1:en + 1]
                    valid = VALID[row:en] & VALID[row + 1:en + 1]
                    lc1, lc2 = LC[row:en], LC[row + 1:en + 1]
                    weights = hz[row + 1:en + 1]
                edge = left ^ right
                for c in eh:
                    if c == "All boundaries":
                        codes = edge
                    else:
                        v1 = np.isin(lc1, CLASS_CODES[c])
                        v2 = np.isin(lc2, CLASS_CODES[c])
                        codes = ((left & ~right) * v2) | ((right & ~left) * v1)
                    rr, cc = np.nonzero((codes > 0) & valid)
                    eh[c] += np.bincount(COUNT[codes[rr, cc]], weights=weights[rr], minlength=10)
        for c, h in eh.items():
            for n in range(1, 10):
                edge_records.append(dict(ssp=s, year=year, landcover=c, n_products=n, segment_length_km=float(h[n])))
        np.savez_compressed(OUT / f"blocks_SSP{s}_{year}.npz", **maps)
        print("Ensemble", s, year, "union", union, "all nine", inter, flush=True)
pd.DataFrame(summary).to_csv(OUT / "ensemble_summary.csv", index=False)
pd.DataFrame(pairs).to_csv(OUT / "pairwise_agreement.csv", index=False)
pd.DataFrame(loo).to_csv(OUT / "leave_one_out_agreement.csv", index=False)
pd.DataFrame(transition).to_csv(OUT / "temporal_transitions.csv", index=False)
pd.DataFrame(edge_records).to_csv(OUT / "direct_boundary_agreement.csv", index=False)
print("COMPLETE ensemble", flush=True)
