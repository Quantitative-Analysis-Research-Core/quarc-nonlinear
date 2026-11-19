"""Delay embedding with arbitrary lags."""

# Copyright (c) 2021-2026 Quantitative Analysis Research Core,
# Center for Human Movement Variability, University of Nebraska at Omaha.
# MIT licence. See LICENSE.txt.

import numpy as np


def embed(z, v=(0, 1, 2), w=None):
    """Delay embedding with arbitrary lags, as MATLAB's embed (one output).

    v is the list of lags or, with w, the number of lags, which are then
    0, w, ..., w*(v-1). Returns a matrix with one row per lag, in increasing
    lag order, and one column per embedded time: row i is z delayed by the
    i-th smallest lag, so the lag-0 row holds the latest sample of each
    vector. This is the transpose of psr's layout.
    """
    if w is not None:
        v = np.arange(int(v)) * w
    z = np.asarray(z, dtype=np.float64).ravel()
    lags = np.sort(np.asarray(v).ravel()).astype(np.int64)
    dim = lags.size
    n = z.size
    if n < lags[-1]:
        raise ValueError("Vector is too small to be embedded with the given lags")
    m = n - (lags[-1] - lags[0])
    t = np.arange(1, m + 1) + lags[-1]                 # MATLAB's 1-based times
    x = np.zeros((dim, m))
    for i in range(dim):
        x[i] = z[t - lags[i] - 1]
    return x
