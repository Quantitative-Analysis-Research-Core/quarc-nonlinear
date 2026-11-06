"""Phase space reconstruction."""

# Copyright (c) 2021-2026 Quantitative Analysis Research Core,
# Center for Human Movement Variability, University of Nebraska at Omaha.
# MIT licence. See LICENSE.txt.

import numpy as np


def psr(x, delay, dim):
    """Delay embedding, as MATLAB's psr.

    x is a vector or an N-by-D matrix (one signal per column). Returns an
    M-by-(dim*D) matrix, M = N - (dim-1)*delay, whose block i (columns
    i*D .. i*D+D-1) is x delayed by i*delay samples.
    """
    x = np.asarray(x, dtype=np.float64)
    if x.ndim == 1:
        x = x[:, None]
    N, D = x.shape
    M = N - (dim - 1) * delay
    if M < 1:
        raise ValueError(f"Series of {N} samples is too short for dim = {dim} at delay = {delay}.")
    return np.hstack([x[i * delay : i * delay + M] for i in range(dim)])
