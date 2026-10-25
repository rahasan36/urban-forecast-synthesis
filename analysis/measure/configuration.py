"""Urban landscape configuration: patch count, mean patch area and edge density of every forecast map.

Patches are eight-neighbor connected urban components, counted globally and within each region
(clipped at region and common-domain boundaries). The global grid is labeled in strips of rows and
the strip seams are stitched across all three eight-neighbor contacts, so the strip counts equal
whole-map counts; verify() checks this on random maps. Mean patch area divides urban area by the
patch count, and edge density divides boundary length by the valid landscape area.

Run after measure_maps.py. Outputs (analysis/output/measure): patches/patches_*.json and
configuration_metrics.csv."""
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio
from scipy import ndimage
from skimage.measure import label

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import COMMON_VALID, PRODUCTS, REGION_RASTER, SSPS, YEARS, forecast_path, output_dir, urban_mask

OUT = output_dir("measure")
PATCHES = OUT / "patches"


def component_counts(u, regions=None, chunk=256):
    counts = np.zeros(21, dtype="int64")
    previous = None
    previous_region = None
    offset = 0
    parent = {}

    def find(x):
        parent.setdefault(x, x)
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for row in range(0, len(u), chunk):
        b = u[row:row + chunk]
        if regions is None:
            labs, n = ndimage.label(b, structure=np.ones((3, 3), dtype="uint8"))
            counts[0] += n
            reg = np.zeros(b.shape, dtype="uint8")
        else:
            reg = regions[row:row + chunk]
            labs, n = label(b.astype("uint8") * reg, connectivity=2, return_num=True)
            rr, cc = np.nonzero(labs)
            ids, ix = np.unique(labs[rr, cc], return_index=True)
            counts += np.bincount(reg[rr[ix], cc[ix]], minlength=21)
        first = labs[0].astype("int64")
        first[first > 0] += offset
        if previous is not None:
            for shift in [-1, 0, 1]:
                if shift == -1:
                    p, q, r, pr = previous[:-1], first[1:], reg[0, 1:], previous_region[:-1]
                elif shift == 1:
                    p, q, r, pr = previous[1:], first[:-1], reg[0, :-1], previous_region[1:]
                else:
                    p, q, r, pr = previous, first, reg[0], previous_region
                yes = (p > 0) & (q > 0) & (r == pr)
                for x, y, k in np.unique(np.column_stack((p[yes], q[yes], r[yes])), axis=0):
                    x = find(int(x))
                    y = find(int(y))
                    if x != y:
                        parent[y] = x
                        counts[int(k)] -= 1
        previous = labs[-1].astype("int64")
        previous[previous > 0] += offset
        previous_region = reg[-1].copy()
        offset += n
    return counts


def verify():
    rng = np.random.default_rng(20260917)
    for n in range(10):
        u = rng.random((31, 43)) > .6
        reg = np.tile(np.repeat(np.arange(1, 5), 11)[:43], (31, 1)).astype("uint8")
        assert component_counts(u, chunk=7)[0] == ndimage.label(u, structure=np.ones((3, 3)))[1]
        a = component_counts(u, reg, chunk=7)
        for k in range(1, 5):
            assert a[k] == ndimage.label(u & (reg == k), structure=np.ones((3, 3)))[1], (a, k)
    print("Strip-wise counts match whole-map eight-neighbor counts on 10 random maps")


if __name__ == "__main__":
    verify()
    PATCHES.mkdir(exist_ok=True)
    with rasterio.open(COMMON_VALID) as r:
        v = r.read(1).astype(bool)
    with rasterio.open(REGION_RASTER) as r:
        reg = r.read(1)
    for prod in PRODUCTS:
        for s in SSPS:
            for y in YEARS:
                p = PATCHES / f"patches_{prod}_SSP{s}_{y}.json"
                if p.exists():
                    continue
                t = time.time()
                with rasterio.open(forecast_path(prod, s, y)) as r:
                    u = urban_mask(r.read(1), prod) & v
                global_n = int(component_counts(u)[0])
                regional = component_counts(u, reg)
                p.write_text(json.dumps(dict(product=prod, ssp=s, year=y, global_patches=global_n,
                                             regional_patches=regional.tolist())))
                print("Configuration", prod, s, y, round(time.time() - t, 1), "s", flush=True)
    g = pd.read_csv(OUT / "global_metrics.csv")
    rr = pd.read_csv(OUT / "regional_metrics.csv")
    m = pd.concat([g, rr])
    m = m[m.threshold == .4]
    rows = []
    for p in PATCHES.glob("patches_*.json"):
        d = json.loads(p.read_text())
        for _, r in m[(m["product"] == d["product"]) & (m.ssp == d["ssp"]) & (m.year == d["year"])].iterrows():
            n = d["global_patches"] if r.region == "Global" else d["regional_patches"][int(r.region_id)]
            rows.append(dict(product=d["product"], ssp=d["ssp"], year=d["year"], region=r.region, patch_count=n,
                             mean_patch_area_ha=r.urban_area_km2 * 100 / n if n else None,
                             edge_density_m_per_ha=r.perimeter_km * 10 / r.valid_area_km2,
                             urban_area_km2=r.urban_area_km2, perimeter_km=r.perimeter_km))
    pd.DataFrame(rows).to_csv(OUT / "configuration_metrics.csv", index=False)
    print("COMPLETE configuration", len(rows), flush=True)
