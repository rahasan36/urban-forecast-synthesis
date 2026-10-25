"""Reduce the SSP1 agreement aggregation to the arrays drawn in Fig. 4.

The aggregation is stored in two parts, each covering half of the grid rows. The
20-km block histograms and the region and country tallies are summed over the two
parts. The output holds:
  hist2100            valid 1-km cells by the number of forecasts calling them urban
                      (0 to 9) in each 20-km block, 2100
  hist2020_global     the same counts summed over the globe, 2020
  inter<Y>_<unit>     pair intersection (36 pairs, itertools.combinations order) by
                      region or country, km2
  area<Y>_<unit>      urban area of each forecast by region or country, km2
"""
import argparse
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))

KEYS = [f"{kind}{year}_{unit}" for kind in ("inter", "area")
        for year in (2020, 2100) for unit in ("region", "country")]


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("parts", nargs=2, help="the two parts of the SSP1 agreement aggregation (.npz)")
    ap.add_argument("--out", default=os.path.join(REPO, "data", "figures", "pairwise_agreement",
                                                  "pairwise_agreement_SSP1.npz"))
    a = ap.parse_args()
    p0, p1 = (np.load(f) for f in a.parts)
    out = {"hist2100": (p0["hist2100"] + p1["hist2100"]).astype(np.int32),
           "hist2020_global": (p0["hist2020"] + p1["hist2020"]).sum(axis=(0, 1)).astype(np.int64)}
    for k in KEYS:
        out[k] = (p0[k] + p1[k]).astype(np.int64)
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    np.savez_compressed(a.out, **out)
    print("wrote", a.out)


if __name__ == "__main__":
    main()
