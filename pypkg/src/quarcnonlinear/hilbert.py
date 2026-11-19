"""Analytic signal by the FFT method."""

# Copyright (c) 2021-2026 Quantitative Analysis Research Core,
# Center for Human Movement Variability, University of Nebraska at Omaha.
# MIT licence. See LICENSE.txt.

import numpy as np


def hilbert(x):
    """Analytic signal, as MATLAB's hilbert: zero the negative frequencies,
    double the positive ones, invert.

    Returns a complex vector whose real part is x and imaginary part its
    Hilbert transform.
    """
    x = np.asarray(x, dtype=np.float64).ravel()
    n = x.size
    h = np.zeros(n)
    if n % 2 == 0:
        h[0] = 1
        h[n // 2] = 1
        h[1 : n // 2] = 2
    else:
        h[0] = 1
        h[1 : (n + 1) // 2] = 2
    return np.fft.ifft(np.fft.fft(x) * h)
