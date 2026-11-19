"""Normalised corrected Shannon entropy of a binary series."""

# Copyright (c) 2021-2026 Quantitative Analysis Research Core,
# Center for Human Movement Variability, University of Nebraska at Omaha.
# MIT licence. See LICENSE.txt.

import numpy as np


def ent_symbolic(x, dim):
    """Normalised corrected Shannon entropy, as MATLAB's ent_symbolic.

    x is a series of 0s and 1s. Words of dim consecutive symbols are read as
    binary numbers; their Shannon entropy, corrected for the number of
    distinct words seen, is divided by its maximum. Returns a value between 0
    and about 1.
    """
    x = np.asarray(x, dtype=np.float64).ravel()
    u = np.unique(x)
    if not np.all(np.isin(u, [0, 1])):
        raise ValueError(
            f"x must be a binary series of 0s and 1s; found {u.size} distinct values "
            f"in [{u.min():g}, {u.max():g}]. Threshold or symbolise the series first, "
            "for example x = (x > np.median(x)).astype(float).")
    nW = x.size - dim + 1
    words = np.zeros(nW)
    for q in range(dim):
        words = words * 2 + x[q:q + nW]
    Sm = 2 ** dim
    P = np.bincount(words.astype(np.int64), minlength=Sm) / nW
    with np.errstate(divide="ignore", invalid="ignore"):
        H = P * np.log2(P)
    H = -np.sum(H[~np.isnan(H)])
    So = np.unique(words).size
    CSE = H + (So - 1) / (2 * Sm * np.log(2))
    CSEm = -np.log2(1 / Sm) + (Sm - 1) / (2 * Sm * np.log(2))
    return float(CSE / CSEm)
