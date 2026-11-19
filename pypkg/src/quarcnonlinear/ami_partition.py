"""Average mutual information by adaptive partition of the rank plane."""

# Copyright (c) 2021-2026 Quantitative Analysis Research Core,
# Center for Human Movement Variability, University of Nebraska at Omaha.
# MIT licence. See LICENSE.txt.

import numpy as np
from scipy.special import chdtri

from ._ami_first_minimum import ami_first_minimum
from ._check_ami_input import check_ami_input
from ._rank_of import rank_of


def ami_partition(x, maxlag, alpha=0.05, fraction=0.2):
    """AMI versus lag 0..maxlag by Darbellay-Vajda adaptive partition.

    A cell of the rank plane is cut into four at its median ranks when a
    chi-square test (3 df, level alpha) rejects a uniform spread; the root is
    always cut. Every lag uses N - maxlag pairs. Returns (tau, curve, info):
    the lag of the first local minimum, a (maxlag+1, 2) array of lag and AMI
    in bits, and a dict with allMinima, usedFallback, fractionLag,
    samplesPerLag, alpha and estimator.
    """
    x = check_ami_input(x, maxlag, "ami_partition")
    if not (0 < alpha < 1):
        raise ValueError("alpha must lie in (0, 1).")
    N = x.size
    m = N - maxlag
    critical = chdtri(3, alpha)          # qchisq(1 - alpha, 3)
    curve = np.column_stack([np.arange(maxlag + 1, dtype=np.float64), np.zeros(maxlag + 1)])
    for lag in range(maxlag + 1):
        curve[lag, 1] = _partition_bits(x[:m], x[lag:m + lag], critical)
    tau, curve, info = ami_first_minimum(curve, fraction, "ami_partition")
    info.update(samplesPerLag=m, alpha=alpha, estimator="adaptive")
    return tau, curve, info


def _partition_bits(a, b, critical):
    """Darbellay-Vajda estimate of I(A;B) in bits, on ranks."""
    n = a.size
    rx = rank_of(a)
    ry = rank_of(b)
    # Explicit stack of cells: bounds (xlo, xhi, ylo, yhi) in ranks, members.
    stack = [((1, n, 1, n), np.arange(n))]
    isRoot = True
    I = 0.0
    while stack:
        cb, idx = stack.pop()
        nc = idx.size
        if nc == 0:
            continue
        wx = cb[1] - cb[0] + 1
        wy = cb[3] - cb[2] + 1
        split = False
        if wx >= 2 and wy >= 2:
            xm = cb[0] + wx // 2 - 1
            ym = cb[2] + wy // 2 - 1
            lx = rx[idx] <= xm
            ly = ry[idx] <= ym
            q = np.array([np.sum(lx & ly), np.sum(lx & ~ly), np.sum(~lx & ly),
                          np.sum(~lx & ~ly)], dtype=np.float64)
            fx = (xm - cb[0] + 1) / wx
            fy = (ym - cb[2] + 1) / wy
            E = nc * np.array([fx * fy, fx * (1 - fy), (1 - fx) * fy, (1 - fx) * (1 - fy)])
            split = isRoot or np.sum((q - E) ** 2 / E) > critical
        isRoot = False
        if split:
            # Pushed in MATLAB's order so the summation order matches.
            stack += [((cb[0], xm, cb[2], ym), idx[lx & ly]),
                      ((cb[0], xm, ym + 1, cb[3]), idx[lx & ~ly]),
                      ((xm + 1, cb[1], cb[2], ym), idx[~lx & ly]),
                      ((xm + 1, cb[1], ym + 1, cb[3]), idx[~lx & ~ly])]
        else:
            I += (nc / n) * np.log2(nc * n / (wx * wy))
    return float(I)
