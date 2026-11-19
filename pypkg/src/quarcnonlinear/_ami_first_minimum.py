"""First local minimum of an AMI curve, with a fraction-of-AMI(0) fallback."""

# Copyright (c) 2021-2026 Quantitative Analysis Research Core,
# Center for Human Movement Variability, University of Nebraska at Omaha.
# MIT licence. See LICENSE.txt.

import warnings

import numpy as np


def ami_first_minimum(curve, fraction, caller):
    """Lag of the first strict local minimum of curve[:, 1].

    Falls back to the first lag at which AMI drops to fraction * AMI(0)
    (Abarbanel et al. 1993). Returns (tau, curve, info); tau is an int, or
    nan (with a warning) when neither rule finds a lag. info holds allMinima
    (rows of curve), usedFallback and fractionLag (nan if never reached).
    """
    v = curve[:, 1]
    k = v.size
    isMin = np.zeros(k, dtype=bool)
    if k >= 3:
        isMin[1:-1] = (v[1:-1] < v[:-2]) & (v[1:-1] <= v[2:])
    idx = np.flatnonzero(isMin)
    info = {"allMinima": curve[idx, :], "usedFallback": False, "fractionLag": np.nan}
    below = np.flatnonzero(v <= fraction * v[0])
    if below.size:
        info["fractionLag"] = int(curve[below[0], 0])
    if idx.size:
        tau = int(curve[idx[0], 0])
    else:
        tau = info["fractionLag"]
        info["usedFallback"] = True
        if np.isnan(tau):
            warnings.warn(f"{caller}: no local minimum in AMI over lags 0:{int(curve[-1, 0])} "
                          f"and the curve never fell to {100 * fraction:.0f}% of its lag-0 "
                          "value.", stacklevel=3)
    return tau, curve, info
