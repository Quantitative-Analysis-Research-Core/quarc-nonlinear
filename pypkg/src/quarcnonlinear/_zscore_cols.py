"""Standardise each column."""

# Copyright (c) 2021-2026 Quantitative Analysis Research Core,
# Center for Human Movement Variability, University of Nebraska at Omaha.
# MIT licence. See LICENSE.txt.

import numpy as np


def zscore_cols(x):
    """Each column minus its mean, over its sd (ddof=1), as R's zscore_cols.

    A vector is treated as one column. Returns a float 2-D array; a constant
    column becomes NaN, as in R and MATLAB.
    """
    x = np.asarray(x, dtype=np.float64)
    if x.ndim == 1:
        x = x[:, None]
    with np.errstate(invalid="ignore", divide="ignore"):
        return (x - x.mean(axis=0)) / x.std(axis=0, ddof=1)
