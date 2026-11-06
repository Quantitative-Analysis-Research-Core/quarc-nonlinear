"""Correlation dimension."""

# Copyright (c) 2021-2026 Quantitative Analysis Research Core,
# Center for Human Movement Variability, University of Nebraska at Omaha.
# MIT licence. See LICENSE.txt.

import numpy as np

from . import _core


def corr_dim(x, delay, dim):
    """Correlation dimension (Grassberger & Procaccia 1983), as MATLAB's corr_dim.

    The correlation sum is computed over 200 log-spaced radii, excluding
    pairs closer than dim*delay samples in time, and the slope is fitted
    over the band of ln C from the midpoint of its range to a quarter-range
    above it. Returns (CoD, extra); extra holds logEps, logC, idx (0-based
    indices of the fitted points), slope, intercept, r2 and bounds.
    """
    x = np.ascontiguousarray(x, dtype=np.float64).ravel()
    bins = 200
    n = x.size - (dim - 1) * delay
    k = dim * delay
    eps2, eps1 = _core.corr_dim_range(x, dim, delay)
    if eps2 == 0:
        eps2 = np.finfo(float).eps
    lo, hi = np.log(eps2), np.log(eps1)
    e = lo + np.arange(bins) * (hi - lo) / (bins - 1)
    e[-1] = hi
    epsilon = np.e ** e
    epsilon[0], epsilon[-1] = eps2, eps1      # exact end radii; see the MATLAB note
    CI = np.asarray(_core.corr_dim_counts(x, dim, delay, epsilon))
    with np.errstate(divide="ignore"):
        CI = np.log(np.cumsum(CI) / (n - k) ** 2)
    CI = CI + np.log(n - k - 1)
    logEps = np.log(epsilon)
    if np.isinf(CI.min()):
        i = np.flatnonzero(CI == -np.inf).max()
        minCI = CI[i + 1 :].min()
    else:
        minCI = CI.min()
    mid = (CI.max() - minCI) / 2 + minCI
    quarter = (CI.max() - minCI) / 4
    bounds = np.array([mid, mid + quarter])
    inside = np.flatnonzero((CI > bounds[0]) & (CI < bounds[1]))
    idx = np.arange(inside.min(), inside.max() + 1)
    A = np.column_stack([logEps[idx], np.ones(idx.size)])
    slope, intercept = np.linalg.lstsq(A, CI[idx], rcond=None)[0]
    yHat = slope * logEps[idx] + intercept
    ssTot = np.sum((CI[idx] - CI[idx].mean()) ** 2)
    r2 = 1 - np.sum((CI[idx] - yHat) ** 2) / ssTot if ssTot > 0 and idx.size > 2 else np.nan
    extra = dict(logEps=logEps, logC=CI, idx=idx, slope=slope, intercept=intercept,
                 r2=r2, bounds=bounds)
    return slope, extra
