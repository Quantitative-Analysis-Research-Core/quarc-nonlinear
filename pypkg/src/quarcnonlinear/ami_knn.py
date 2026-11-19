"""Average mutual information by the Kraskov-Stoegbauer-Grassberger estimator."""

# Copyright (c) 2021-2026 Quantitative Analysis Research Core,
# Center for Human Movement Variability, University of Nebraska at Omaha.
# MIT licence. See LICENSE.txt.

import numpy as np
from scipy.special import digamma

from ._ami_first_minimum import ami_first_minimum
from ._check_ami_input import check_ami_input


def ami_knn(x, maxlag, k=10, fraction=0.2, chunk=512):
    """AMI versus lag 0..maxlag by Kraskov-Stoegbauer-Grassberger, algorithm 1.

    The series is standardised and moved by 1e-10 * (((1..N) * g) % 1 - 0.5),
    g the golden-ratio conjugate, to break ties deterministically, as in
    MATLAB. Every lag uses N - maxlag pairs; cost O(maxlag * N^2). Can return
    small negative values where the true AMI is zero. Returns (tau, curve,
    info): the lag of the first local minimum, a (maxlag+1, 2) array of lag
    and AMI in bits, and a dict with allMinima, usedFallback, fractionLag,
    samplesPerLag, k and estimator.
    """
    x = check_ami_input(x, maxlag, "ami_knn")
    N = x.size
    m = N - maxlag
    if k >= m:
        raise ValueError(f"ami_knn: k ({int(k)}) must be smaller than the number of pairs "
                         f"per lag ({m}).")
    s = np.std(x, ddof=1)
    if s > 0:
        x = (x - x.mean()) / s
    x = x + 1e-10 * ((np.arange(1, N + 1) * 0.6180339887498949) % 1 - 0.5)
    curve = np.column_stack([np.arange(maxlag + 1, dtype=np.float64), np.zeros(maxlag + 1)])
    for lag in range(maxlag + 1):
        curve[lag, 1] = _ksg_bits(x[:m], x[lag:m + lag], int(k), int(chunk))
    tau, curve, info = ami_first_minimum(curve, fraction, "ami_knn")
    info.update(samplesPerLag=m, k=k, estimator="knn")
    return tau, curve, info


def _ksg_bits(a, b, k, chunk):
    """Kraskov-Stoegbauer-Grassberger algorithm 1, in bits."""
    n = a.size
    nx = np.zeros(n)
    ny = np.zeros(n)
    for lo in range(0, n, chunk):
        hi = min(lo + chunk, n)
        dx = np.abs(a[lo:hi, None] - a[None, :])
        dy = np.abs(b[lo:hi, None] - b[None, :])
        d = np.maximum(dx, dy)
        d[np.arange(hi - lo), np.arange(lo, hi)] = np.inf      # exclude self
        e = np.partition(d, k - 1, axis=1)[:, k - 1:k]
        # Strictly closer than eps_i; the self term is counted and removed.
        nx[lo:hi] = (dx < e).sum(axis=1) - 1
        ny[lo:hi] = (dy < e).sum(axis=1) - 1
    return float((digamma(k) + digamma(n) - np.mean(digamma(nx + 1) + digamma(ny + 1)))
                 / np.log(2))
