"""Detrended fluctuation analysis."""

# Copyright (c) 2021-2026 Quantitative Analysis Research Core,
# Center for Human Movement Variability, University of Nebraska at Omaha.
# MIT licence. See LICENSE.txt.

import numpy as np


def dfa(x, scales, order=1):
    """Detrended fluctuation analysis (Peng et al. 1995), as MATLAB's dfa.

    The profile (cumulative sum of x minus its mean) is cut into
    non-overlapping windows of each scale, a polynomial of the given order is
    removed from each, and the root mean square of the residuals is averaged
    over windows. Returns (scales, fluctuation, alpha): alpha is the slope of
    log fluctuation on log scale. Plotting is omitted.
    """
    x = np.asarray(x, dtype=np.float64).ravel()
    scales = np.asarray(scales, dtype=np.float64).ravel()
    profile = np.cumsum(x - x.mean())
    fluctuation = np.zeros(scales.size)
    for s, scale in enumerate(scales):
        scale = int(scale)
        chunks = x.size // scale
        idx = np.arange(1, scale + 1, dtype=np.float64)
        V = idx[:, None] ** np.arange(order, -1, -1)   # Vandermonde, as polyfit
        Q = np.linalg.qr(V)[0]
        segs = profile[: chunks * scale].reshape(chunks, scale)
        resid = segs - (segs @ Q) @ Q.T
        fluctuation[s] = np.sqrt(np.sum(np.mean(resid ** 2, axis=1)) / chunks)
    ls = np.log(scales)
    with np.errstate(divide="ignore", invalid="ignore"):
        lf = np.log(fluctuation)
    keep = ~np.isnan(lf)
    A = np.column_stack([ls[keep], np.ones(keep.sum())])
    alpha = np.linalg.lstsq(A, lf[keep], rcond=None)[0][0]
    return scales, fluctuation, alpha
