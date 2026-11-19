"""Cross approximate entropy."""

# Copyright (c) 2021-2026 Quantitative Analysis Research Core,
# Center for Human Movement Variability, University of Nebraska at Omaha.
# MIT licence. See LICENSE.txt.

import numpy as np


def ent_xap(x, y, dim, radius, k):
    """Cross approximate entropy, in nats, as MATLAB's ent_xap.

    Both series are standardised. Templates are similar when every component
    differs by strictly less than radius (in standard deviations). k is the
    step between the two template lengths compared, dim and dim + k.
    """
    x = np.asarray(x, dtype=np.float64).ravel()
    y = np.asarray(y, dtype=np.float64).ravel()
    x = (x - np.mean(x)) / np.std(x, ddof=1)
    y = (y - np.mean(y)) / np.std(y, ddof=1)
    N = x.size
    radius = np.std(x, ddof=1) * radius
    Cm = np.zeros(2)
    for j in range(2):
        m = dim + j * k
        n = N - m + 1
        C = np.empty(n)
        for i in range(n):
            ok = np.abs(x[i] - y[:n]) < radius
            for q in range(1, m):
                ok &= np.abs(x[i + q] - y[q:n + q]) < radius
            C[i] = np.sum(ok) / n
        Cm[j] = np.sum(C) / n
    with np.errstate(divide="ignore", invalid="ignore"):
        return float(np.log(Cm[0] / Cm[1]))
