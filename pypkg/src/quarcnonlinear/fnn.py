"""False nearest neighbours."""

# Copyright (c) 2021-2026 Quantitative Analysis Research Core,
# Center for Human Movement Variability, University of Nebraska at Omaha.
# MIT licence. See LICENSE.txt.

import warnings

import numpy as np

from . import _core


def fnn(x, delay, maxdim, Rtol=15, Atol=2, speed=0):
    """Embedding dimension by false nearest neighbours, as MATLAB's fnn.

    Returns (dim, dE): the selected dimension and the fraction of false
    neighbours at each dimension 1..maxdim. With speed=1 the search stops at
    the first dimension that meets a stopping rule, and later entries of dE
    stay 0. The neighbour search is the compiled k-d tree that MATLAB's fnn
    uses, so ties resolve as in MATLAB.
    """
    x = np.ascontiguousarray(x, dtype=np.float64).ravel()
    n = x.size - delay * maxdim
    if n < 2:
        raise ValueError("Series too short for this delay and maxdim.")
    RA = np.std(x, ddof=1)
    dE = np.zeros(maxdim)
    dim = None
    for j in range(1, maxdim + 1):
        dE[j - 1] = _core.fnn_false_fraction(x, delay, j, n, Rtol, Atol, RA)
        if speed == 1:
            if dE[j - 1] <= 0:
                dim = j
                break
            if j >= 3 and dE[j - 3] > dE[j - 2] < dE[j - 1]:
                dim = j - 1
                break
            if j >= 2 and abs(dE[j - 1] - dE[j - 2]) <= 0.001:
                dim = j
                break
    if speed == 0 and maxdim >= 3:
        for i in range(2, maxdim):          # MATLAB's i, 1-based
            if dE[i - 1] <= 0 or dE[i - 2] > dE[i - 1] < dE[i]:
                dim = i
                break
            if dE[i - 2] - dE[i - 1] <= 0.001:
                dim = i - 1
                break
    if dim is None:
        dim = maxdim
        warnings.warn("no dimension found, dim set to maxdim")
    return dim, dE
