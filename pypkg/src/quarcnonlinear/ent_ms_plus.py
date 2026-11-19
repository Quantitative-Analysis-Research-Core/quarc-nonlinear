"""Multiscale entropies: refined composite, composite, standard, fuzzy and generalised."""

# Copyright (c) 2021-2026 Quantitative Analysis Research Core,
# Center for Human Movement Variability, University of Nebraska at Omaha.
# MIT licence. See LICENSE.txt.

import numpy as np


def ent_ms_plus(x, scale, dim, radius):
    """Multiscale entropies at scales 1..scale, as MATLAB's ent_ms_plus.

    Returns (RCMSE, CMSE, MSE, MSFE, GMSE), five arrays of length scale.
    MSE, CMSE, RCMSE and MSFE use the tolerance radius * std(x) at every
    scale. GMSE coarse-grains by variance, so its tolerance is radius times
    the SD of the scale-4 variance series, fixed across scales (Costa &
    Goldberger 2015). GMSE[0] (scale 1) is NaN: the variance of a single
    sample is undefined.
    """
    x = np.asarray(x, dtype=np.float64).ravel()
    R = radius * np.std(x, ddof=1)
    N = x.size
    nw4 = N // 4
    Rg = radius * np.std(np.var(x[:4 * nw4].reshape(nw4, 4), axis=1, ddof=1), ddof=1)
    RCMSE, CMSE, MSE, MSFE, GMSE = (np.zeros(scale) for _ in range(5))
    with np.errstate(divide="ignore", invalid="ignore"):
        for i in range(1, scale + 1):
            nj = N // i
            # Row k, column j: the window starting at j*i + k, length i.
            # Windows running past the end are NaN, as MATLAB's try/catch.
            o2 = np.full((i, nj), np.nan)
            y = np.full((i, nj), np.nan)
            for j in range(nj):
                for k in range(i):
                    lo = j * i + k
                    if lo + i <= N:
                        w = x[lo:lo + i]
                        o2[k, j] = np.var(w, ddof=1) if i > 1 else 0.0
                        y[k, j] = np.sum(w) / i
            GMSE[i - 1] = np.nan if i == 1 else _ms_samp_ent(o2[0], dim, Rg)[0]
            row1 = y[0, ~np.isnan(y[0])]
            MSE[i - 1] = _ms_samp_ent(row1, dim, R)[0]
            MSFE[i - 1] = _ms_fuzzy_ent(row1, dim, R, 2)
            nm = np.zeros(i)
            nm1 = np.zeros(i)
            for k in range(i):
                _, nm[k], nm1[k] = _ms_samp_ent(y[k, ~np.isnan(y[k])], dim, R)
                CMSE[i - 1] += 1 / i * -np.log(nm1[k] / nm[k])
            RCMSE[i - 1] = -np.log((1 / i * np.sum(nm1)) / (1 / i * np.sum(nm)))
    return RCMSE, CMSE, MSE, MSFE, GMSE


def _ms_samp_ent(data, m, r):
    """Sample entropy as ent_ms_plus computes it: N - m templates, self-matches
    removed. Returns (SE, sum_nm, sum_nm1), the last two the B and A totals."""
    N = data.size
    n = N - m
    Bm = np.empty(n)
    Am = np.empty(n)
    sum_nm = sum_nm1 = 0.0
    for i in range(n):
        dj = np.abs(data[:n] - data[i])
        for k in range(1, m):
            dj = np.maximum(dj, np.abs(data[k:n + k] - data[i + k]))
        dj1 = np.maximum(dj, np.abs(data[m:n + m] - data[i + m]))
        nmi = np.sum(dj <= r) - 1
        nm1i = np.sum(dj1 <= r) - 1
        sum_nm += nmi
        sum_nm1 += nm1i
        Bm[i] = nmi / n
        Am[i] = nm1i / n
    return -np.log((np.sum(Am) / n) / (np.sum(Bm) / n)), sum_nm, sum_nm1


def _ms_fuzzy_ent(series, dim, r, n):
    """Fuzzy entropy (Chen et al. 2007) as ent_ms_plus computes it; r is in
    the units of the series."""
    N = series.size
    phi = np.zeros(2)
    for j in range(2):
        m = dim + j
        nP = N - m + 1
        P = np.vstack([series[q:N - m + q + 1] for q in range(m)])
        P = P - P.mean(axis=0)
        aux = np.zeros(nP)
        for i in range(N - m):
            d = np.max(np.abs(P - P[:, [i]]), axis=0)
            simi = np.exp(-(d ** n) / r)
            aux[i] = (np.sum(simi) - 1) / (N - m - 1)
        phi[j] = np.sum(aux) / (N - m)
    return np.log(phi[0]) - np.log(phi[1])
