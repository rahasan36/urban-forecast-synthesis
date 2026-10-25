"""Comparison of urban area on each forecast's full source domain and on the common analysis domain.

Assembles the outputs of full_domain_areas.py (and gao_land_weighted.py when present) into one row
per analyzed map and one per distributed base-year map, with the area outside the common domain and
its share of the full-source area, and checks the common-domain areas against measure_maps.py.
Reported values: the Li 2021 country demand sums from the authors' Cntry_SSP_UrbanTS.xlsx (146
countries; model demand, not allocated raster area) when the workbook is under
data/rasters/source/Li2021, and the Gao 2021 base-year total of 610,742 km2 (Gao and O'Neill, 2020,
Supplementary Table 1). Other totals were not located and stay empty.
Outputs (analysis/output/area_reconciliation): area_comparison_forecasts.csv and
area_comparison_baselines.csv"""
import json
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import LABEL, OUTPUT, PRODUCTS, SOURCE, SSPS, YEARS, output_dir

OUT = output_dir("area_reconciliation")
calc = [json.loads(p.read_text()) for p in sorted((OUT / "maps").glob("*.json"))]
calc = {(c["product"], c["kind"], c["ssp"], c["year"]): c for c in calc}
gao_file = OUT / "gao_land_weighted.json"
gao = {(g["kind"], g["ssp"], g["year"]): g["derived_native_land_weighted_km2"] for g in json.loads(gao_file.read_text())} if gao_file.exists() else {}
reported = {}
li = SOURCE / "Li2021" / "Cntry_SSP_UrbanTS.xlsx"
if li.exists():
    import openpyxl
    wb = openpyxl.load_workbook(li, data_only=True, read_only=True)
    for s in SSPS:
        rows = list(wb[f"SSP{s}"].values)
        hdr = rows[0]
        assert len(rows[1:]) == 146
        for y in YEARS:
            c = hdr.index("y" + str(y))
            reported[("Li2021", s, y)] = sum(float(r[c]) for r in rows[1:])
g = pd.read_csv(OUTPUT / "measure" / "global_metrics.csv")
measured = {(r["product"], r.ssp, r.year): r.urban_area_km2 for _, r in g[g.threshold == .4].iterrows()}

rows = []
for p in PRODUCTS:
    for s in SSPS:
        for y in YEARS:
            c = calc[(p, "forecast", s, y)]
            full, common = c["raster_full_km2"], c["raster_common_km2"]
            assert abs(common - measured[(p, s, y)]) <= 1, (p, s, y, common, measured[(p, s, y)])
            derived = gao.get(("forecast", s, y)) if p == "Gao2021" else c["fraction_full_km2"] if p == "He2023" else None
            rows.append(dict(forecast=LABEL[p], product=p, ssp=s, year=y, reported_km2=reported.get((p, s, y)),
                             full_source_km2=full, common_domain_km2=common, outside_common_domain_km2=full - common,
                             excluded_percent=100 * (full - common) / full, fraction_weighted_km2=derived,
                             urban_rule=c["urban_rule"]))
f = pd.DataFrame(rows)
f.to_csv(OUT / "area_comparison_forecasts.csv", index=False)
base = []
for p in PRODUCTS:
    hit = [c for k, c in calc.items() if k[0] == p and k[1] == "baseline"]
    if not hit:
        base.append(dict(forecast=LABEL[p], product=p, base_year=None, reported_km2=None, full_source_km2=None,
                         common_domain_km2=None, excluded_percent=None, fraction_weighted_km2=None))
        continue
    c = hit[0]
    derived = gao.get(("baseline", None, c["year"])) if p == "Gao2021" else c["fraction_full_km2"] if p == "He2023" else None
    base.append(dict(forecast=LABEL[p], product=p, base_year=c["year"], reported_km2=610742 if p == "Gao2021" else None,
                     full_source_km2=c["raster_full_km2"], common_domain_km2=c["raster_common_km2"],
                     excluded_percent=100 * (c["raster_full_km2"] - c["raster_common_km2"]) / c["raster_full_km2"],
                     fraction_weighted_km2=derived))
b = pd.DataFrame(base)
b.to_csv(OUT / "area_comparison_baselines.csv", index=False)
top = f.loc[f.excluded_percent.idxmax()]
print(f"excluded share across the 81 maps: median {f.excluded_percent.median():.2f}%, maximum "
      f"{top.excluded_percent:.2f}% ({top.forecast}, SSP{top.ssp}, {top.year}); base-year maps: median "
      f"{b.excluded_percent.median():.2f}%, maximum {b.excluded_percent.max():.2f}%")
