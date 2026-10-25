"""Tables S10 and S11 from the outputs of observed_2020.py and builtup_2020.py.

For each forecast and scenario: mapped urban area, intersection with the reference, Jaccard index
(intersection / union), share of the reference recovered (intersection / reference) and the ratio of
mapped to reference area. The built-up reference is thresholded at a built-up fraction of 0.13
(--threshold). The Spearman correlation compares the forecasts' Jaccard indices under the two
references in each scenario.
Outputs (analysis/output/benchmarks): table_observed_2020.csv and table_builtup_2020.csv"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import LABEL, PRODUCTS, SSPS, output_dir

ap = argparse.ArgumentParser()
ap.add_argument("--threshold", type=float, default=0.13)
a = ap.parse_args()
OUT = output_dir("benchmarks")
o = json.loads((OUT / "observed_2020_counts.json").read_text())
g = json.loads((OUT / "builtup_2020_histograms.json").read_text())
k0 = int(round(a.threshold * 100))
ref = int(np.array(g["hist"])[k0:].sum())


def row(p, s, area, inter, total):
    return dict(forecast=LABEL[p], ssp=s, mapped_km2=area, intersection_km2=inter,
                jaccard=round(inter / (area + total - inter), 3), reference_recovered_percent=round(100 * inter / total, 1),
                mapped_to_reference=round(area / total, 3))


obs_rows, ghs_rows = [], []
for s in SSPS:
    for p in PRODUCTS:
        c = o["c"][f"{p}|SSP{s}"]
        obs_rows.append(row(p, s, c["area"], c["inter"], o["obs"]))
        h = np.array(g["ph"][f"{p}|SSP{s}"])
        ghs_rows.append(row(p, s, int(h.sum()), int(h[k0:].sum()), ref))
t10, t11 = pd.DataFrame(obs_rows), pd.DataFrame(ghs_rows)
t10.to_csv(OUT / "table_observed_2020.csv", index=False)
t11.to_csv(OUT / "table_builtup_2020.csv", index=False)
print("observed urban land", o["obs"], "km2; built-up reference at", a.threshold, ":", ref, "km2")
def jaccard(t, total):
    return t.intersection_km2 / (t.mapped_km2 + total - t.intersection_km2)


for s in SSPS:
    rho = spearmanr(jaccard(t10[t10.ssp == s], o["obs"]), jaccard(t11[t11.ssp == s], ref))[0]
    print(f"SSP{s}: Spearman correlation of Jaccard indices under the two references {rho:.3f}")
print(t10.to_string(index=False))
print(t11.to_string(index=False))
