"""QUARC Nonlinear: Python bindings to the shared C++ core.

Function names, arguments and defaults follow the MATLAB functions of the
same name, which are the specification.
"""

# Copyright (c) 2021-2026 Quantitative Analysis Research Core,
# Center for Human Movement Variability, University of Nebraska at Omaha.
# MIT licence. See LICENSE.txt.

import numpy as np

from . import _core

__all__ = ["ent_samp"]


def ent_samp(x, dim, radius=0.2, flag="prop"):
    """Sample entropy (Richman & Moorman 2000), as MATLAB's ent_samp.

    If flag is "prop", radius is a proportion of std(x) (N-1 normalisation,
    as MATLAB); otherwise it is used as given. Returns inf when no pair of
    templates of length dim matches.
    """
    x = np.ascontiguousarray(x, dtype=np.float64).ravel()
    if x.size < dim + 1:
        raise ValueError("The time series must have at least m+1 data points.")
    if str(flag).lower() == "prop":
        radius = radius * np.std(x, ddof=1)
    return _core.ent_samp(x, int(dim), float(radius))
