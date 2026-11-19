"""Weighted recurrence entropy."""

# Copyright (c) 2021-2026 Quantitative Analysis Research Core,
# Center for Human Movement Variability, University of Nebraska at Omaha.
# MIT licence. See LICENSE.txt.

import numpy as np


def ent_weighted(wrp):
    """Weighted recurrence entropy, in nats, as MATLAB's ent_weighted.

    Shannon entropy of the column sums of the square weighted recurrence
    matrix wrp, binned into 50 equal-width bins and weighted by strength.
    Returns 0 when every column sum is equal.
    """
    si = np.atleast_2d(np.asarray(wrp, dtype=np.float64)).sum(axis=0)
    lo, hi = si.min(), si.max()
    if hi == lo:
        return 0.0
    width = (hi - lo) / 49
    S = np.sum(si)
    # Bins 1-49 by index; bin 50 holds the maximum alone, as in MATLAB.
    b = np.minimum(np.floor((si - lo) / width), 48) + 1
    b[si == hi] = 50
    p1 = np.array([np.sum(si[b == k]) / S for k in range(1, 51)])
    with np.errstate(divide="ignore", invalid="ignore"):
        pp = p1 * np.log(p1)
    pp[np.isnan(pp)] = 0
    return float(-np.sum(pp))
