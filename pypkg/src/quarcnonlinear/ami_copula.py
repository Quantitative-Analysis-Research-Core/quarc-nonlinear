"""Average mutual information by the Gaussian copula estimator."""

# Copyright (c) 2021-2026 Quantitative Analysis Research Core,
# Center for Human Movement Variability, University of Nebraska at Omaha.
# MIT licence. See LICENSE.txt.

import numpy as np
from scipy.special import ndtri

from ._ami_first_minimum import ami_first_minimum
from ._check_ami_input import check_ami_input


def ami_copula(x, maxlag, fraction=0.2):
    """AMI versus lag 0..maxlag by the Gaussian copula (Ince et al. 2017).

    Ranks are mapped to normal scores, then AMI = -0.5 * log2(1 - r^2), r
    capped at 1 - 1e-6. Exact for a Gaussian copula and a lower bound
    otherwise; blind to non-monotonic dependence (near zero on a sine at a
    quarter-period lag). Every lag uses N - maxlag pairs. Returns (tau,
    curve, info): the lag of the first local minimum, a (maxlag+1, 2) array
    of lag and AMI in bits, and a dict with allMinima, usedFallback,
    fractionLag, samplesPerLag and estimator.
    """
    x = check_ami_input(x, maxlag, "ami_copula")
    N = x.size
    m = N - maxlag
    curve = np.column_stack([np.arange(maxlag + 1, dtype=np.float64), np.zeros(maxlag + 1)])
    for lag in range(maxlag + 1):
        a = _normal_scores(x[:m])
        b = _normal_scores(x[lag:m + lag])
        r = np.sum(a * b) / (np.sqrt(np.sum(a ** 2)) * np.sqrt(np.sum(b ** 2)))
        r = min(max(r, -0.999999), 0.999999)
        curve[lag, 1] = -0.5 * np.log2(1 - r ** 2)
    tau, curve, info = ami_first_minimum(curve, fraction, "ami_copula")
    info.update(samplesPerLag=m, estimator="copula")
    return tau, curve, info


def _normal_scores(v):
    """Ranks (ties by position) mapped to standard normal quantiles."""
    n = v.size
    u = np.empty(n)
    u[np.argsort(v, kind="stable")] = np.arange(1, n + 1) / (n + 1)
    return ndtri(u)
