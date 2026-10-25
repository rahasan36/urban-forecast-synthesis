"""Endpoint changes from 2020 to 2100 and ensemble summaries of boundary contact.

Each forecast is kept as its own member; the spread among forecasts is described, never treated as
a sample. Ensemble means of contact use four weightings: the nine forecasts with equal weight, the
three forecasts sharing the FLUS allocation model (Chen 2020, Chen 2022, Hou 2022) at one-third each
(Table S8), the original four forecasts, and SSP1 without Zhang 2023.

Run after measure_maps.py. Outputs (analysis/output/measure): endpoint_changes.csv,
landcover_endpoint_changes.csv (net overlap increase and its share, Table S4), landcover_rankings.csv
and ensemble_contact_changes.csv."""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import GRID, PRODUCTS, output_dir

OUT = output_dir("measure")
g = pd.read_csv(OUT / "global_metrics.csv")
r = pd.read_csv(OUT / "regional_metrics.csv")
lc = pd.read_csv(OUT / "landcover_metrics.csv")


def changes(df, values, extra=[]):
    keys = ["product", "ssp", "threshold", "region"] + extra
    b = df[df.year == 2020][keys + values].rename(columns={v: v + "_2020" for v in values})
    q = df[df.year == 2100].merge(b, on=keys, validate="one_to_one")
    for v in values:
        q[v + "_change"] = q[v] - q[v + "_2020"]
        q[v + "_growth_percent"] = 100 * q[v + "_change"] / q[v + "_2020"].replace(0, np.nan)
    return q


ep = changes(pd.concat([g, r]), ["urban_area_km2", "perimeter_km"])
ep.to_csv(OUT / "endpoint_changes.csv", index=False)
cl = changes(lc, ["overlap_km2", "contact_km", "contact_share_percent"], ["landcover"])
den = (cl[cl.landcover != "Baseline urban"].groupby(["product", "ssp", "threshold", "region"])
       .overlap_km2_change.sum().rename("nonurban_net_increment_km2").reset_index())
cl = cl.merge(den, on=["product", "ssp", "threshold", "region"])
cl["conversion_increment_share_percent"] = 100 * cl.overlap_km2_change / cl.nonurban_net_increment_km2.replace(0, np.nan)
cl.to_csv(OUT / "landcover_endpoint_changes.csv", index=False)
rank = lc[(lc.threshold == .4) & (lc.landcover != "Baseline urban")].copy()
rank["rank"] = rank.groupby(["product", "ssp", "year", "region"]).contact_km.rank(method="min", ascending=False)
rank.to_csv(OUT / "landcover_rankings.csv", index=False)
primary = cl[cl.threshold == .4]
rows = []
for (s, region, c), q in primary.groupby(["ssp", "region", "landcover"]):
    assert len(q) == 9
    v = q.contact_km.to_numpy()
    b = q.contact_km_2020.to_numpy()
    groups = {"Nine products": np.ones(9),
              "FLUS down-weighted": np.where(q["product"].isin(["Chen2020", "Chen2022", "Hou2022"]), 1 / 3, 1),
              "Original four": q["product"].isin(PRODUCTS[:4]).astype(float),
              "SSP1 without Zhang": (q["product"] != "Zhang2023").astype(float)}
    for name, w in groups.items():
        if name == "SSP1 without Zhang" and s != 1:
            continue
        base = np.average(b, weights=w)
        end = np.average(v, weights=w)
        rows.append(dict(ssp=s, region=region, landcover=c, ensemble=name, mean_contact_2020_km=base,
                         mean_contact_2100_km=end, mean_length_growth_percent=100 * (end / base - 1) if base else np.nan,
                         minimum_2100_km=v[w > 0].min(), maximum_2100_km=v[w > 0].max(),
                         spread_ratio=v[w > 0].max() / v[w > 0].min() if v[w > 0].min() > 0 else np.nan,
                         products_increasing=int(((v > b) & (w > 0)).sum()), n_products=int((w > 0).sum())))
pd.DataFrame(rows).to_csv(OUT / "ensemble_contact_changes.csv", index=False)
counts = dict(n_products=9, product_scenarios=27, focal_maps=81, landcover_rows=len(lc),
              common_valid_km2=json.loads(GRID.read_text())["valid_cells"])
(OUT / "analysis_counts.json").write_text(json.dumps(counts, indent=2))
print("COMPLETE summaries", counts)
