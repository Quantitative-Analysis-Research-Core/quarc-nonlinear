"""Fractional Gaussian noise."""

# Copyright (c) 2021-2026 Quantitative Analysis Research Core,
# Center for Human Movement Variability, University of Nebraska at Omaha.
# MIT licence. See LICENSE.txt.

import numpy as np


def fgn_sim(n, H, rng=None):
    """Fractional Gaussian noise by Davies-Harte circulant embedding, as
    MATLAB's fgn_sim.

    Returns a vector of length n with Hurst exponent H in (0, 1), zero mean
    and unit variance. rng is a numpy Generator or seed (MATLAB and R use
    their global generators, so values never match theirs).
    """
    rng = np.random.default_rng(rng)
    n = int(n)
    z = rng.standard_normal(2 * n)
    zr = z[:n].copy()
    zi = z[n:].copy()
    zic = -zi
    zi[0] = 0
    zr[0] = zr[0] * np.sqrt(2)
    zi[n - 1] = 0
    zr[n - 1] = zr[n - 1] * np.sqrt(2)
    back = np.arange(n - 2, 0, -1)                      # MATLAB's (n-1):-1:2, 0-based
    zr = np.concatenate([zr, zr[back]])
    zi = np.concatenate([zi, zic[back]])
    z = zr + 1j * zi
    k = np.arange(n, dtype=np.float64)
    gammak = (np.abs(k - 1) ** (2 * H) - 2 * np.abs(k) ** (2 * H) + np.abs(k + 1) ** (2 * H)) / 2
    ind = np.concatenate([np.arange(n - 1), [n - 1], np.arange(n - 2, 0, -1)])
    # MATLAB's ifft(v)*length(z): numpy's ifft is normalised, as MATLAB's.
    gk = np.real(np.fft.ifft(gammak[ind])) * z.size
    if not np.all(gk > 0):
        raise ValueError("Re(gk)-vector not positive")
    z = z * np.sqrt(gk)
    z = np.fft.ifft(z) * z.size
    z = 0.5 * (n - 1) ** (-0.5) * z
    return np.real(z[:n])
