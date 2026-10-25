"""Table S9: Kendall's coefficient of concordance of regional ranks across the nine forecasts.

For urban area (from measure/measure_maps.py) and for cropland and forest contact (from
measure/zone_contact.py by region), the 19 reporting regions are ranked within each forecast and
W = 12 S / (m^2 (n^3 - n)), with m = 9 forecasts, n = 19 regions and S the sum of squared deviations
of the regions' rank sums from their mean. W ranges from 0 (no common ordering) to 1 (identical
orderings).

Run after measure_maps.py and zone_contact.py region <SSP> <year> for the nine scenario-year
combinations. Output: analysis/output/agreement/table_concordance.csv"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import OUTPUT, PRODUCTS, REGIONS, REGIONS_MAIN, SSPS, YEARS, output_dir


def kendall_w(M):
    """M: regions x forecasts."""
    r = np.apply_along_axis(rankdata, 0, M)
    m, n = M.shape[1], M.shape[0]
    Ri = r.sum(1)
    return 12 * ((Ri - Ri.mean()) ** 2).sum() / (m ** 2 * (n ** 3 - n))


rg = pd.read_csv(OUTPUT / "measure" / "regional_metrics.csv")
rg = rg[rg.threshold == .4]
rows = []
for s in SSPS:
    row = dict(metric="Urban area", ssp=s)
    for y in YEARS:
        piv = rg[(rg.ssp == s) & (rg.year == y)].pivot(index="region", columns="product", values="urban_area_km2")
        row[str(y)] = round(kendall_w(piv.loc[[REGIONS[r] for r in REGIONS_MAIN], PRODUCTS].values), 3)
    rows.append(row)
for k, metric in ((1, "Cropland contact"), (2, "Forest contact")):
    for s in SSPS:
        row = dict(metric=metric, ssp=s)
        for y in YEARS:
            acc = np.load(OUTPUT / "measure" / f"zone_contact_region_SSP{s}_{y}.npz")["acc"]
            row[str(y)] = round(kendall_w(acc[:, REGIONS_MAIN, k].T), 3)
        rows.append(row)
table = pd.DataFrame(rows)
table.to_csv(output_dir("agreement") / "table_concordance.csv", index=False)
print(table.to_string(index=False))
