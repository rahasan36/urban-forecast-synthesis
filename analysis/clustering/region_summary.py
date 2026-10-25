"""Regional shares of contested and consensus clusters under both specifications (main text and
Fig. S12) and the Spearman correlation of the regional rankings between them.

Regions with at least 15 eligible blocks under the primary specification are reported. Run after
primary_clusters.py and strict_clusters.py.
Output: analysis/output/clustering/cluster_regions.csv"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import REGIONS, REGIONS_MAIN, output_dir

OUT = output_dir("clustering")
G = np.load(OUT / "gistar_SSP1_2100.npz")
zD, BC, mask = G["zD"], G["BC"], G["mask"]
reg = np.load(OUT / "block_regions.npz")["region_all"]
HOT = (zD >= 1.96) & mask
COLD = (zD <= -1.96) & mask
S = np.load(OUT / "support_clusters.npz")
cat, smask, sreg = S["category"], S["mask"], S["region"]
rows = []
for c in REGIONS_MAIN:
    s = mask & (reg == c)
    t = smask & (sreg == c)
    if s.sum() < 15 or t.sum() == 0:
        continue
    rows.append(dict(region_id=c, region=REGIONS[c], primary_blocks=int(s.sum()),
                     primary_contested_percent=100 * (HOT & s).sum() / s.sum(),
                     primary_consensus_percent=100 * (COLD & s).sum() / s.sum(),
                     strict_blocks=int(t.sum()), strict_lower_support_percent=100 * np.sum(t & (cat == -1)) / t.sum(),
                     strict_higher_support_percent=100 * np.sum(t & (cat == 1)) / t.sum()))
d = pd.DataFrame(rows).sort_values("primary_contested_percent", ascending=False)
d.to_csv(OUT / "cluster_regions.csv", index=False)
print("eligible blocks %d, contested %.1f%%, consensus %.1f%%" % (mask.sum(), 100 * HOT.sum() / mask.sum(), 100 * COLD.sum() / mask.sum()))
print("Spearman, contested vs lower support: %.2f" % spearmanr(d.primary_contested_percent, d.strict_lower_support_percent)[0])
print("Spearman, consensus vs higher support: %.2f" % spearmanr(d.primary_consensus_percent, d.strict_higher_support_percent)[0])
print(d.round(1).to_string(index=False))
