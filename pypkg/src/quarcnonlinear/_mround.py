"""MATLAB's round, which takes halves away from zero."""

# Copyright (c) 2021-2026 Quantitative Analysis Research Core,
# Center for Human Movement Variability, University of Nebraska at Omaha.
# MIT licence. See LICENSE.txt.

import numpy as np


def mround(v):
    """Round halves away from zero, elementwise (np.round goes to even).

    Returns a float array for array input and a float for scalar input; wrap
    in int() where an index is needed.
    """
    v = np.asarray(v, dtype=np.float64)
    out = np.sign(v) * np.floor(np.abs(v) + 0.5)
    return out.item() if out.ndim == 0 else out
