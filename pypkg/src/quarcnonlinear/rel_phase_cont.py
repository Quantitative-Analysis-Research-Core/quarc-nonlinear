"""Continuous relative phase."""

# Copyright (c) 2021-2026 Quantitative Analysis Research Core,
# Center for Human Movement Variability, University of Nebraska at Omaha.
# MIT licence. See LICENSE.txt.

import numpy as np

from .hilbert import hilbert


def rel_phase_cont(x, y, fs):
    """Continuous relative phase between two oscillators, as MATLAB's
    rel_phase_cont.

    Each series is centred on the midpoint of its range. Returns
    (crp, crpH) in degrees: crp from position and velocity (one shorter than
    the input), crpH from the analytic signal (Hilbert transform).
    """
    x = np.asarray(x, dtype=np.float64).ravel()
    y = np.asarray(y, dtype=np.float64).ravel()
    x = x - x.min() - (x.max() - x.min()) / 2
    y = y - y.min() - (y.max() - y.min()) / 2
    x0, y0 = x[:-1], y[:-1]
    dx = np.diff(x) * fs
    dy = np.diff(y) * fs
    hx = np.imag(hilbert(x))
    hy = np.imag(hilbert(y))
    with np.errstate(divide="ignore", invalid="ignore"):
        crp = np.arctan((dx * y0 - dy * x0) / (x0 * y0 - dy * dx)) * 180 / np.pi
        crpH = np.arctan((hx * y - hy * x) / (x * y + hx * hy)) * 180 / np.pi
    return crp, crpH
