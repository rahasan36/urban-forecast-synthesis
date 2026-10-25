"""Urban area, overlap with the 2016 land-cover baseline and boundary contact of every forecast map.

For each forecast, scenario, year and threshold (0.30, 0.40 and 0.50 for the fractional forecasts,
0.40 elsewhere), the script counts urban cells by region and baseline class and measures the
boundary: every horizontal or vertical grid side that joins an urban and a nonurban cell, both in
the common valid domain, contributes its WGS84 ellipsoidal ground length and is assigned to the
region and baseline class of the nonurban cell. Sides that touch an invalid cell are excluded and
their length is recorded separately.

Outputs (analysis/output/measure): maps/metrics_<forecast>_SSP<s>_<year>_<threshold>.json,
global_metrics.csv, regional_metrics.csv, landcover_metrics.csv and baseline_availability.json."""
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio
from pyproj import Geod, Transformer

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import (CLASS_CODES, CLASSES, COMMON_VALID, FRACTIONAL, GRID, LANDCOVER_2016, PRODUCTS, REGION_RASTER,
                    REGIONS, REGIONS_MAIN, SSPS, YEARS, forecast_path, output_dir, urban_mask)

OUT = output_dir("measure")
MAPS = OUT / "maps"
MAPS.mkdir(exist_ok=True)
grid = json.loads(GRID.read_text())
H, W = grid["height"], grid["width"]
with rasterio.open(LANDCOVER_2016) as r:
    LC = r.read(1)
    TR = r.transform
with rasterio.open(REGION_RASTER) as r:
    REG = r.read(1)
with rasterio.open(COMMON_VALID) as r:
    VALID = r.read(1).astype(bool)
lookup = np.full(256, 13, dtype="uint8")
for k, (name, codes) in enumerate(CLASS_CODES.items()):
    lookup[codes] = k
lookup[190] = 12
LC = lookup[LC]
assert np.all(LC[VALID] < 13)
regarea = np.bincount(REG[VALID], minlength=21)

inv = Transformer.from_crs(6933, 4326, always_xy=True)
geod = Geod(ellps="WGS84")
y = TR.f + np.arange(H + 1) * TR.e
lo, la = inv.transform(np.zeros(H + 1), y)
lo1, la1 = inv.transform(np.ones(H + 1) * 1000, y)
horizontal = np.asarray(geod.inv(lo, la, lo1, la1)[2]) / 1000
vertical = np.asarray(geod.inv(lo[:-1], la[:-1], lo[1:], la[1:])[2]) / 1000


def edges(u):
    """Boundary length by region and class of the nonurban cell, excluded length and side count."""
    out = np.zeros((21, 13), dtype="float64")
    excluded = 0.
    total_count = 0
    for i in range(0, H, 256):
        j = min(i + 256, H)
        a = u[i:j]
        v = VALID[i:j]
        cl = LC[i:j]
        reg = REG[i:j]
        dif = a[:, :-1] != a[:, 1:]
        ok = v[:, :-1] & v[:, 1:]
        rr, cc = np.nonzero(dif & ok)
        nn = cc + a[rr, cc].astype("int32")
        keys = reg[rr, nn].astype("int32") * 13 + cl[rr, nn]
        out += np.bincount(keys, weights=vertical[i + rr], minlength=21 * 13).reshape(21, 13)
        total_count += len(rr)
        rr, cc = np.nonzero(dif & (v[:, :-1] ^ v[:, 1:]))
        excluded += float(vertical[i + rr].sum())
        end = min(j, H - 1)
        a = u[i:end]
        b = u[i + 1:end + 1]
        dif = a != b
        ok = VALID[i:end] & VALID[i + 1:end + 1]
        rr, cc = np.nonzero(dif & ok)
        nn = i + rr + a[rr, cc].astype("int32")
        keys = REG[nn, cc].astype("int32") * 13 + LC[nn, cc]
        out += np.bincount(keys, weights=horizontal[i + rr + 1], minlength=21 * 13).reshape(21, 13)
        total_count += len(rr)
        rr, cc = np.nonzero(dif & (VALID[i:end] ^ VALID[i + 1:end + 1]))
        excluded += float(horizontal[i + rr + 1].sum())
    return out, excluded, total_count


def measure(prod, s, yr, threshold):
    dest = MAPS / f"metrics_{prod}_SSP{s}_{yr}_{threshold:.1f}.json"
    if dest.exists():
        return
    t = time.time()
    with rasterio.open(forecast_path(prod, s, yr)) as r:
        u = urban_mask(r.read(1), prod, threshold)
    uv = u & VALID
    area = np.bincount(REG[uv], minlength=21)
    overlap = np.bincount(REG[uv].astype("int32") * 13 + LC[uv], minlength=21 * 13).reshape(21, 13)
    e, excluded, count = edges(u)
    assert int(overlap.sum()) == int(area.sum())
    obj = dict(product=prod, ssp=s, year=yr, threshold=threshold, area_by_region=area.tolist(),
               overlap_by_region_class=overlap.tolist(), edge_by_region_class=e.tolist(),
               excluded_boundary_km=excluded, shared_edge_count=count)
    dest.write_text(json.dumps(obj), encoding="utf-8")
    print("Measured", prod, s, yr, threshold, round(time.time() - t, 1), "s", flush=True)


for prod in PRODUCTS:
    for s in SSPS:
        for yr in YEARS:
            for th in ([.3, .4, .5] if prod in FRACTIONAL else [.4]):
                measure(prod, s, yr, th)

global_rows, region_rows, land_rows = [], [], []
for f in sorted(MAPS.glob("metrics_*.json")):
    d = json.loads(f.read_text())
    key = {k: d[k] for k in ["product", "ssp", "year", "threshold"]}
    area = np.asarray(d["area_by_region"])
    overlap = np.asarray(d["overlap_by_region_class"])
    edge = np.asarray(d["edge_by_region_class"])
    for ri in [-1] + REGIONS_MAIN:
        ar = area.sum() if ri == -1 else area[ri]
        ov = overlap.sum(0) if ri == -1 else overlap[ri]
        ed = edge.sum(0) if ri == -1 else edge[ri]
        record = dict(**key, region="Global" if ri == -1 else REGIONS[ri], region_id=ri, urban_area_km2=int(ar),
                      perimeter_km=float(ed.sum()), valid_area_km2=int(regarea.sum() if ri == -1 else regarea[ri]))
        if ri == -1:
            record.update(excluded_boundary_km=d["excluded_boundary_km"], shared_edge_count=d["shared_edge_count"])
            global_rows.append(record)
        else:
            region_rows.append(record)
        for c, name in enumerate(CLASSES):
            land_rows.append(dict(**record, landcover=name, overlap_km2=int(ov[c]), contact_km=float(ed[c]),
                                  contact_share_percent=float(ed[c] / ed.sum() * 100) if ed.sum() else None))
pd.DataFrame(global_rows).to_csv(OUT / "global_metrics.csv", index=False)
pd.DataFrame(region_rows).to_csv(OUT / "regional_metrics.csv", index=False)
pd.DataFrame(land_rows).to_csv(OUT / "landcover_metrics.csv", index=False)
(OUT / "baseline_availability.json").write_text(json.dumps(dict(
    region_area=regarea.tolist(), classes=CLASSES, class_area=np.bincount(LC[VALID], minlength=13).tolist(),
    regions=REGIONS), indent=2))
print("COMPLETE metrics", len(global_rows), len(region_rows), len(land_rows), flush=True)
