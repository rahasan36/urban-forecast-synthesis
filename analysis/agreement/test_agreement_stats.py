"""Synthetic checks of agreement_stats.py: the count-based pooled Jaccard index and the summed pair
intersection equal the direct pairwise calculation, and the location agreement is about 1 for
identical or nested forecasts, about 0 for independent placement (also under a regional gradient)
and negative for disjoint placement."""
import itertools
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from agreement_stats import jaccard_hist, kappa, strata

rng = np.random.default_rng(1)
H, W = 200, 300
key = (np.arange(H)[:, None] // 50 * 6 + np.arange(W)[None, :] // 50).astype(np.int64) * np.ones((H, W), np.int64)
K = 24
v = np.ones((H, W), bool)
v[:5] = False


def direct_J(masks):
    I = sum((a & b & v).sum() for a, b in itertools.combinations(masks, 2))
    Un = sum(((a | b) & v).sum() for a, b in itertools.combinations(masks, 2))
    return I / Un


def hist(masks):
    n = sum(m.astype(int) for m in masks)
    return np.bincount(n[v], minlength=10)


masks = [rng.random((H, W)) < p for p in np.linspace(0.05, 0.4, 9)]
print("pooled Jaccard from counts vs direct", round(jaccard_hist(hist(masks)), 6), round(direct_J(masks), 6))
s = strata(key, v, masks, K)
print("summed pair intersection", s["O"].sum(), sum((a & b & v).sum() for a, b in itertools.combinations(masks, 2)))
print("independent placement kappa", round(kappa(s), 4))
base = rng.random((H, W)) < 0.2
print("identical kappa", round(kappa(strata(key, v, [base] * 9, K)), 4))
z = rng.random((H, W))
print("nested kappa", round(kappa(strata(key, v, [z < p for p in np.linspace(0.05, 0.4, 9)], K)), 4))
r = rng.integers(0, 9, (H, W))
print("disjoint kappa", round(kappa(strata(key, v, [(r == i) & (rng.random((H, W)) < 0.9) for i in range(9)], K)), 4))
g = (np.arange(W)[None, :] // 50) / 6 * 0.5 * np.ones((H, W))
print("regional gradient, random placement kappa", round(kappa(strata(key, v, [rng.random((H, W)) < g for _ in range(9)], K)), 4))
