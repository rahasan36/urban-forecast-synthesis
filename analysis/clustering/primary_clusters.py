"""Divergent urban boundaries and high spatial uncertainty in global urban land forecasts

Author: Rakibul Ahasan
"""
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import REGION_CODES, output_dir
from gistar import gistar

OUT = output_dir("clustering")
z1 = np.load(OUT / "block_sums_SSP1_2100.npz")
n20 = np.load(OUT / "block_sums_SSP1_2020.npz")
U, NS, BC, BF, BT = z1["union"], z1["nsum"], z1["bc"], z1["bf"], z1["bt"]
U20, NS20 = n20["union"], n20["nsum"]
reg = np.load(OUT / "block_regions.npz")["region_all"]
MINU = 10.0
mask = U >= MINU
A = np.where(mask, NS / np.maximum(9 * U, 1e-9), np.nan)
D = 1.0 - A
m20 = U20 >= MINU
A20 = np.where(m20, NS20 / np.maximum(9 * U20, 1e-9), np.nan)
print("blocks in the eligible set, 2100:", int(mask.sum()), "2020:", int(m20.sum()))
print("median support 2100 %.3f, 2020 %.3f" % (np.nanmedian(A), np.nanmedian(A20)))
zD = gistar(np.nan_to_num(D), mask, r=4)
zC = gistar(BC, mask, r=4)
np.savez_compressed(OUT / "gistar_SSP1_2100.npz", zD=zD, zC=zC, A=A, D=D, U=U, BC=BC, BF=BF, BT=BT, mask=mask)
HOT = (zD >= 1.96) & mask
COLD = (zD <= -1.96) & mask
CH = (zC >= 1.96) & mask
both = HOT & CH
print("contested blocks %d (%.1f%%), consensus blocks %d (%.1f%%)" % (
    HOT.sum(), 100 * HOT.sum() / mask.sum(), COLD.sum(), 100 * COLD.sum() / mask.sum()))
tot = BC[mask].sum()
print("cropland boundary of the union in the eligible set: %.0f thousand km" % (tot / 1e3))
print("  in contested clusters %.1f%%, in consensus clusters %.1f%%" % (100 * BC[HOT].sum() / tot, 100 * BC[COLD].sum() / tot))
codes = json.loads(REGION_CODES.read_text(encoding="utf-8"))["code"]
inv = {v: k for k, v in codes.items()}
rows = []
for c, nm in sorted(inv.items()):
    s = mask & (reg == c)
    if s.sum() < 30:
        continue
    rows.append(dict(region=nm, blocks=int(s.sum()), median_support=float(np.nanmedian(A[s])),
                     contested_percent=100 * float((HOT & s).sum()) / s.sum(),
                     consensus_percent=100 * float((COLD & s).sum()) / s.sum(),
                     cropland_boundary_thousand_km=float(BC[s].sum() / 1e3),
                     cropland_boundary_in_contested_cropland_hotspots_percent=100 * float(BC[both & s].sum()) / max(BC[s].sum(), 1e-9)))
rows.sort(key=lambda r: r["median_support"])
(OUT / "primary_clusters_regions.json").write_text(json.dumps(rows, indent=1))
for r in rows:
    print("%-32s %6d %7.3f %7.1f %7.1f" % (r["region"], r["blocks"], r["median_support"], r["contested_percent"], r["consensus_percent"]))
