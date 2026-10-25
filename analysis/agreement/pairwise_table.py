"""Table S12: pairwise agreement and the quantity and allocation components of disagreement.

Pooled values sum intersections, unions and each disagreement component over the 36 forecast pairs
before dividing (Pontius and Millones, 2011); the pairwise columns summarize the individual pairs.
The note gives the amount share of the shortfall of the pooled index under SSP1 in 2100 for zones of
5 to 100 km and its regional range. Run after ensemble.py and amount_placement.py.
Outputs: analysis/output/agreement/table_pairwise_agreement.csv and table_pairwise_agreement_note.txt"""
import json
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import output_dir

OUT = output_dir("agreement")
d = pd.read_csv(OUT / "pairwise_agreement.csv")
rows = []
for (s, y), g in d.groupby(["ssp", "year"]):
    q = 100 * g.quantity_disagreement_km2.sum() / g.disagreement_km2.sum()
    j, qp = g.jaccard, g.quantity_share_percent
    rows.append(dict(ssp=s, year=y, pooled_jaccard=round(g.intersection_km2.sum() / g.union_km2.sum(), 3),
                     pairwise_jaccard_median=round(j.median(), 2), pairwise_jaccard_min=round(j.min(), 2),
                     pairwise_jaccard_max=round(j.max(), 2), quantity_percent=round(q, 1),
                     allocation_percent=round(100 - q, 1), pairwise_quantity_median_percent=round(qp.median(), 1),
                     pairwise_quantity_min_percent=round(qp.min(), 1), pairwise_quantity_max_percent=round(qp.max(), 1)))
pd.DataFrame(rows).to_csv(OUT / "table_pairwise_agreement.csv", index=False)

n3 = json.loads((OUT / "amount_placement_SSP1_2100.json").read_text())
zs = {k: 100 * v[0] for k, v in n3["zone_amount_placement_share"].items()}
reg = {r: 100 * (1 - c) / (1 - o) for r, (o, c, e) in n3["regions_J_Jmax_Jind"].items()}
lo, hi = min(reg, key=reg.get), max(reg, key=reg.get)
note = (f"Under SSP1 in 2100, differences in the amount of urban land each forecast assigns to a 20-km zone account for "
        f"{zs['20']:.1f} percent of the shortfall of the pooled index from full agreement. The amount share is "
        f"{zs['5']:.1f}, {zs['10']:.1f}, {zs['50']:.1f}, and {zs['100']:.1f} percent in zones of 5, 10, 50, and 100 km, "
        f"and ranges across regions from {reg[lo]:.1f} percent in {lo} to {reg[hi]:.1f} percent in {hi}.")
(OUT / "table_pairwise_agreement_note.txt").write_text(note + "\n", encoding="utf-8")
print(pd.DataFrame(rows).to_string(index=False))
print(note)
