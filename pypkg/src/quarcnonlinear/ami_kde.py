"""Average mutual information by Gaussian kernel density estimation."""

# Copyright (c) 2021-2026 Quantitative Analysis Research Core,
# Center for Human Movement Variability, University of Nebraska at Omaha.
# MIT licence. See LICENSE.txt.

import numpy as np

from ._ami_first_minimum import ami_first_minimum
from ._check_ami_input import check_ami_input


def ami_kde(x, maxlag, fraction=0.2, chunk=512, bandwidth="fixed"):
    """AMI versus lag 0..maxlag from Gaussian kernel densities (Thomas et al. 2014).

    Lag k uses all N - k pairs; cost is O(maxlag * N^2). chunk is the number
    of rows of the pairwise kernel evaluated at once (memory only).
    bandwidth "adaptive" scales each kernel by sqrt(g / pilot_i), g the
    geometric mean of a fixed-bandwidth pilot (Abramson 1982). Returns
    (tau, curve, info): the lag of the first local minimum, a (maxlag+1, 2)
    array of lag and AMI in bits, and a dict with allMinima, usedFallback,
    fractionLag, bandwidth and estimator.
    """
    if bandwidth not in ("fixed", "adaptive"):
        raise ValueError('bandwidth must be "fixed" or "adaptive".')
    adaptive = bandwidth == "adaptive"
    x = check_ami_input(x, maxlag, "ami_kde")
    N = x.size
    chunk = int(chunk)
    curve = np.column_stack([np.arange(maxlag + 1, dtype=np.float64), np.zeros(maxlag + 1)])
    curve[0, 1] = _kde_mi_bits(x, x, chunk, adaptive)
    for lag in range(1, maxlag + 1):
        curve[lag, 1] = _kde_mi_bits(x[:N - lag], x[lag:], chunk, adaptive)
    tau, curve, info = ami_first_minimum(curve, fraction, "ami_kde")
    info.update(bandwidth=bandwidth, estimator="kde-adaptive" if adaptive else "kde")
    return tau, curve, info


def _kde_mi_bits(X, Y, chunk, adaptive=False):
    n = X.size
    hx = np.std(X, ddof=1) / n ** (1 / 6)
    hy = np.std(Y, ddof=1) / n ** (1 / 6)
    if hx == 0 or hy == 0:
        return 0.0
    pX = _kde_marginal(X, hx, chunk)
    pY = _kde_marginal(Y, hy, chunk)
    pXY = _kde_joint(X, Y, hx, hy, chunk)
    if adaptive:
        # The fixed-bandwidth estimates above are the pilots.
        pX = _kde_marginal(X, hx, chunk, _abramson(pX))
        pY = _kde_marginal(Y, hy, chunk, _abramson(pY))
        pXY = _kde_joint(X, Y, hx, hy, chunk, _abramson(pXY))
    ok = (pXY > 0) & (pX > 0) & (pY > 0)
    return float(np.sum(np.log2(pXY[ok] / (pX[ok] * pY[ok]))) / n)


def _kde_marginal(v, h, chunk, lam=None):
    n = v.size
    p = np.zeros(n)
    cst = 1 / (n * h * np.sqrt(2 * np.pi))
    for lo in range(0, n, chunk):
        rows = slice(lo, min(lo + chunk, n))
        z = (v[rows, None] - v[None, :]) / h
        if lam is None:
            p[rows] = cst * np.exp(-0.5 * z ** 2).sum(axis=1)
        else:
            p[rows] = cst * (np.exp(-0.5 * (z / lam) ** 2) / lam).sum(axis=1)
    return p


def _kde_joint(X, Y, hx, hy, chunk, lam=None):
    n = X.size
    dx = X - X.mean()
    dy = Y - Y.mean()
    r = np.sum(dx * dy) / (np.sqrt(np.sum(dx ** 2)) * np.sqrt(np.sum(dy ** 2)))
    r = min(max(r, -0.999999), 0.999999)
    omr2 = 1 - r ** 2
    cst = 1 / (n * 2 * np.pi * hx * hy * np.sqrt(omr2))
    p = np.zeros(n)
    for lo in range(0, n, chunk):
        rows = slice(lo, min(lo + chunk, n))
        u = (X[rows, None] - X[None, :]) / hx
        v = (Y[rows, None] - Y[None, :]) / hy
        q = (u ** 2 - 2 * r * (u * v) + v ** 2) / omr2
        if lam is None:
            p[rows] = cst * np.exp(-0.5 * q).sum(axis=1)
        else:
            L2 = lam ** 2
            p[rows] = cst * (np.exp(-0.5 * q / L2) / L2).sum(axis=1)
    return p


def _abramson(pilot):
    """Abramson's square-root law: lambda_i = sqrt(g / pilot_i)."""
    pilot = np.maximum(pilot, np.finfo(float).tiny)
    g = np.exp(np.mean(np.log(pilot)))
    return np.sqrt(g / pilot)
