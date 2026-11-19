"""Ranks with ties broken by position."""

# Copyright (c) 2021-2026 Quantitative Analysis Research Core,
# Center for Human Movement Variability, University of Nebraska at Omaha.
# MIT licence. See LICENSE.txt.

import numpy as np


def rank_of(v):
    """Ranks 1..n, ties broken by position, as MATLAB's stable sort gives."""
    v = np.asarray(v).ravel()
    r = np.empty(v.size, dtype=np.int64)
    r[np.argsort(v, kind="stable")] = np.arange(1, v.size + 1)
    return r
