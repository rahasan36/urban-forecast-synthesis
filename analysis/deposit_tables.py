"""Divergent urban boundaries and high spatial uncertainty in global urban land forecasts

Author: Rakibul Ahasan
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import LABEL, OUTPUT, REPO

THRESHOLD = 0.4
CLASSES = {"Cropland": "cropland", "Forest": "forest"}
EXCLUDED_REGIONS = {"Global", "Greenland, Bermuda + others"}


def area_tables(area_dir, out):
    forecasts = pd.read_csv(area_dir / "area_comparison_forecasts.csv")
    forecasts = forecasts.drop(columns=["product"], errors="ignore")
    forecasts["excluded_percent"] = forecasts["excluded_percent"].round(3)
    forecasts.to_csv(out / "urban_area_by_domain.csv", index=False)
    baselines = pd.read_csv(area_dir / "area_comparison_baselines.csv")
    baselines = baselines.drop(columns=["product"], errors="ignore")
    baselines["excluded_percent"] = baselines["excluded_percent"].round(3)
    baselines.to_csv(out / "urban_area_base_year_maps.csv", index=False)
    return forecasts


def regional_contact(endpoint_file, out):
    d = pd.read_csv(endpoint_file)
    d = d[(d.threshold == THRESHOLD) & (d.year == 2100) & d.landcover.isin(CLASSES)
          & ~d.region.isin(EXCLUDED_REGIONS)]
    rows = []
    for (region, ssp), g in d.groupby(["region", "ssp"], sort=True):
        row = {"region": region, "ssp": int(ssp)}
        for landcover, key in CLASSES.items():
            c = g[g.landcover == landcover]
            start = c.contact_km_2020.to_numpy(float)
            end = c.contact_km.to_numpy(float)
            ratio = end[start > 0] / start[start > 0]
            row[f"{key}_ratio_median"] = np.median(ratio)
            row[f"{key}_ratio_min"] = ratio.min()
            row[f"{key}_ratio_max"] = ratio.max()
            row[f"{key}_forecasts_increasing"] = int((ratio > 1).sum())
            row[f"{key}_spread_2100"] = end.max() / end[end > 0].min()
        rows.append(row)
    table = pd.DataFrame(rows).round(3)
    table.to_csv(out / "regional_contact_change.csv", index=False)
    return table


def intersections(membership_file, out):
    m = pd.read_csv(membership_file)
    m["forecasts"] = m["products"].apply(lambda s: " + ".join(LABEL[p] for p in s.split(" + ")))
    m = m.rename(columns={"code": "membership_code", "n_products": "n_forecasts"})
    m = m[["ssp", "year", "membership_code", "n_forecasts", "forecasts", "area_km2"]]
    m.to_csv(out / "exact_intersections.csv", index=False)
    return m


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--area_dir", default=OUTPUT / "area_reconciliation", type=Path)
    ap.add_argument("--endpoint", default=OUTPUT / "measure" / "landcover_endpoint_changes.csv", type=Path)
    ap.add_argument("--membership", default=OUTPUT / "agreement" / "membership_exact.csv", type=Path)
    ap.add_argument("--out", default=REPO / "data" / "tables", type=Path)
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    f = area_tables(a.area_dir, a.out)
    print(f"excluded share across {len(f)} maps: median {f.excluded_percent.median():.2f}%, "
          f"maximum {f.excluded_percent.max():.2f}%")
    t = regional_contact(a.endpoint, a.out)
    both = ((t.cropland_ratio_min > 1) & (t.forest_ratio_min > 1)).sum()
    print(f"{len(t)} region-pathway rows; cropland median ratio {t.cropland_ratio_median.min():.2f} to "
          f"{t.cropland_ratio_median.max():.2f}, forest {t.forest_ratio_median.min():.2f} to "
          f"{t.forest_ratio_median.max():.2f}, cropland spread {t.cropland_spread_2100.min():.1f} to "
          f"{t.cropland_spread_2100.max():.1f}")
    m = intersections(a.membership, a.out)
    print(f"{len(m)} intersection rows")


if __name__ == "__main__":
    main()
