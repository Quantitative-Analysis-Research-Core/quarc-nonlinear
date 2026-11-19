"""Cross sample entropy."""

# Copyright (c) 2021-2026 Quantitative Analysis Research Core,
# Center for Human Movement Variability, University of Nebraska at Omaha.
# MIT licence. See LICENSE.txt.

import numpy as np


def ent_xsamp(x, y, dim, radius, norm):
    """Cross sample entropy (Richman & Moorman 2000), in nats, as MATLAB's ent_xsamp.

    norm 1 rescales both series to 0-1 and uses radius times the mean of
    their SDs; norm 2 standardises both and uses radius as is; anything else
    leaves the data alone and uses radius times the mean SD of the raw series.
    """
    x = np.asarray(x, dtype=np.float64).ravel()
    y = np.asarray(y, dtype=np.float64).ravel()
    if x.size != y.size:
        raise ValueError("The data series need to be the same length!")
    N = x.size
    if norm == 1:
        xn = (x - x.min()) / (x.max() - x.min())
        yn = (y - y.min()) / (y.max() - y.min())
        r = radius * ((np.std(xn, ddof=1) + np.std(yn, ddof=1)) / 2)
    elif norm == 2:
        xn = (x - np.mean(x)) / np.std(x, ddof=1)
        yn = (y - np.mean(y)) / np.std(y, ddof=1)
        r = radius
    else:
        xn, yn = x, y
        r = radius * ((np.std(x, ddof=1) + np.std(y, ddof=1)) / 2)
    n = N - dim
    Bm = np.empty(n)
    Am = np.empty(n)
    for i in range(n):
        dj = np.abs(xn[:n] - yn[i])
        for k in range(1, dim):
            dj = np.maximum(dj, np.abs(xn[k:n + k] - yn[i + k]))
        dj1 = np.maximum(dj, np.abs(xn[dim:n + dim] - yn[i + dim]))
        Bm[i] = np.sum(dj <= r) / n
        Am[i] = np.sum(dj1 <= r) / n
    with np.errstate(divide="ignore", invalid="ignore"):
        return float(-np.log((np.sum(Am) / n) / (np.sum(Bm) / n)))
