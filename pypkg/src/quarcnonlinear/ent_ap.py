"""Approximate entropy."""

# Copyright (c) 2021-2026 Quantitative Analysis Research Core,
# Center for Human Movement Variability, University of Nebraska at Omaha.
# MIT licence. See LICENSE.txt.

import numpy as np


def ent_ap(x, dim, radius):
    """Approximate entropy (Pincus 1991), in nats, as MATLAB's ent_ap.

    radius is a proportion of std(x) (N-1 normalisation). Self-matches are
    counted, as in the original definition.
    """
    x = np.asarray(x, dtype=np.float64).ravel()
    radius = radius * np.std(x, ddof=1)
    N = x.size
    phim = np.zeros(2)
    for j in range(2):
        m = dim + j
        n = N - m + 1
        phi = np.empty(n)
        for i in range(n):
            d = np.abs(x[:n] - x[i])
            for k in range(1, m):
                d = np.maximum(d, np.abs(x[k:n + k] - x[i + k]))
            phi[i] = np.sum(d <= radius) / n
        phim[j] = np.sum(np.log(phi)) / n
    return float(phim[0] - phim[1])
