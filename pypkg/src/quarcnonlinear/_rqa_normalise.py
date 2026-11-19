"""Distance normalisation for the RQA functions."""

# Copyright (c) 2021-2026 Quantitative Analysis Research Core,
# Center for Human Movement Variability, University of Nebraska at Omaha.
# MIT licence. See LICENSE.txt.

import numpy as np


def rqa_normalise(a, norm, scaleDim):
    """Scale a negated distance matrix a by norm: "euc", "min", "max" or "none".

    "euc" divides by sqrt(mean distance^2 + 2*scaleDim), "min" by the
    smallest nonzero distance and "max" by the largest. Returns the scaled
    matrix (a itself for "none").
    """
    a = np.asarray(a, dtype=np.float64)
    neg = a[a < 0]
    if norm == "euc":
        b = neg.mean() if neg.size else np.nan
        b = -np.sqrt(abs(b ** 2 + 2 * scaleDim))
        return a / abs(b)
    if norm == "min":
        return a / abs(neg.max() if neg.size else -np.inf)
    if norm == "max":
        return a / abs(neg.min() if neg.size else np.inf)
    return a
