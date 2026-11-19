"""Input checks shared by the AMI estimators."""

# Copyright (c) 2021-2026 Quantitative Analysis Research Core,
# Center for Human Movement Variability, University of Nebraska at Omaha.
# MIT licence. See LICENSE.txt.

import numpy as np


def check_ami_input(x, maxlag, caller):
    """Return x as a flat float array; raise ValueError on NaN or too large a maxlag."""
    x = np.asarray(x, dtype=np.float64).ravel()
    if np.isnan(x).any():
        raise ValueError(f"{caller}: x contains NA. Mutual information is undefined for "
                         "missing data; remove or impute first.")
    if maxlag >= x.size:
        raise ValueError(f"{caller}: maxlag ({int(maxlag)}) must be smaller than the "
                         f"number of samples ({x.size}).")
    return x
