"""Divergent urban boundaries and high spatial uncertainty in global urban land forecasts

Author: Rakibul Ahasan
"""
import numpy as np
from scipy.signal import fftconvolve


def disk(r):
    y, x = np.mgrid[-r:r + 1, -r:r + 1]
    return ((x * x + y * y) <= r * r).astype(np.float64)


def gistar(x, mask, r=4):
    K = disk(r)
    m = mask.astype(np.float64)
    xv = np.where(mask, x, 0.0)
    S1 = fftconvolve(xv, K, mode="same")
    W = fftconvolve(m, K, mode="same")
    n = float(mask.sum())
    Xb = xv.sum() / n
    S = np.sqrt((xv ** 2).sum() / n - Xb ** 2)
    num = S1 - Xb * W
    den = S * np.sqrt(np.maximum((n * W - W ** 2) / (n - 1.0), 0.0))
    z = np.where((den > 0) & mask, num / np.maximum(den, 1e-12), np.nan)
    return np.where(mask, z, np.nan)


def brute(x, mask, r):
    K = disk(r)
    H, Wd = x.shape
    xv = np.where(mask, x, 0.0)
    n = float(mask.sum())
    Xb = xv.sum() / n
    S = np.sqrt((xv ** 2).sum() / n - Xb ** 2)
    z = np.full(x.shape, np.nan)
    for i in range(H):
        for j in range(Wd):
            if not mask[i, j]:
                continue
            s1 = 0.0
            w = 0.0
            for di in range(-r, r + 1):
                for dj in range(-r, r + 1):
                    if K[di + r, dj + r] == 0:
                        continue
                    ii, jj = i + di, j + dj
                    if 0 <= ii < H and 0 <= jj < Wd and mask[ii, jj]:
                        s1 += x[ii, jj]
                        w += 1
            d = S * np.sqrt((n * w - w * w) / (n - 1.0))
            z[i, j] = (s1 - Xb * w) / d if d > 0 else np.nan
    return z


if __name__ == "__main__":
    rng = np.random.default_rng(7)
    g = rng.normal(10, 1, (60, 80))
    mk = np.ones((60, 80), bool)
    g[12:20, 12:20] += 8.0
    g[40:48, 55:63] -= 8.0
    zf = gistar(g, mk, r=3)
    zb = brute(g, mk, r=3)
    print("max abs difference from brute force:", float(np.nanmax(np.abs(zf - zb))))
    print("z at the implanted high cluster center:", round(float(zf[16, 16]), 2))
    print("z at the implanted low cluster center:", round(float(zf[44, 59]), 2))
    print("z at a background point:", round(float(zf[5, 70]), 2))
    z2 = gistar(rng.normal(0, 1, (60, 80)), mk, r=3)
    print("share |z| > 2.58 in a pure-noise field:", round(float(np.nanmean(np.abs(z2) > 2.58)), 4))
