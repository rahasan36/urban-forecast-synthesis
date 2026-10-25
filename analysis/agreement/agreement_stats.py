"""Agreement statistics of nine binary forecast maps within strata (zones or zone x distance band).

For each stratum: A valid cells, S1 = sum of u_i and S2 = sum of u_i^2 (u_i = cells forecast i calls
urban), O = sum over cells of C(n, 2) = the summed intersection of the 36 forecast pairs, M = the sum
over pairs of min(u_i, u_j), the largest pair intersection the amounts allow, and Un = cells at least
one forecast calls urban. The pair intersection expected under independent placement inside a
stratum is E = (S1^2 - S2) / (2A), and the location agreement kappa = sum(O - E) / sum(M - E).
The pooled pairwise Jaccard index from the histogram h[n] of the number n of forecasts calling a cell
urban is J = sum h[n] C(n, 2) / sum h[n] (36 - C(9 - n, 2)) (Eqs. S6 to S10)."""
import numpy as np


def C2(n):
    return n * (n - 1) / 2.0


def strata(key, v, masks, K):
    kv = key[v]
    A = np.bincount(kv, minlength=K).astype(float)
    U = np.stack([np.bincount(key[m & v], minlength=K) for m in masks]).astype(float)
    n = np.zeros(key.shape, np.int64)
    for m in masks:
        n += m
    O = np.bincount(kv, weights=C2(n[v]), minlength=K)
    Un = np.bincount(key[v & (n > 0)], minlength=K).astype(float)
    S1 = U.sum(0)
    S2 = (U ** 2).sum(0)
    Us = np.sort(U, 0)
    k = len(masks)
    M = sum(Us[i] * (k - 1 - i) for i in range(k))
    return dict(A=A, S1=S1, S2=S2, O=O, M=M, Un=Un)


def kappa(s, idx=slice(None)):
    A = s["A"][idx]
    E = np.where(A > 0, (s["S1"][idx] ** 2 - s["S2"][idx]) / (2 * np.maximum(A, 1)), 0)
    return (s["O"][idx] - E).sum() / max((s["M"][idx] - E).sum(), 1e-12)


def jaccard_hist(h, axis=-1, k=9):
    n = np.arange(h.shape[axis])
    sh = [1] * h.ndim
    sh[axis] = -1
    n = n.reshape(sh)
    return (h * C2(n)).sum(axis) / np.maximum((h * (C2(k) - C2(k - n))).sum(axis), 1e-12)
