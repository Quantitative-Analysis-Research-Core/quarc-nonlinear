"""Surrogate data, Theiler's algorithms 0, 1 and 2."""

# Copyright (c) 2021-2026 Quantitative Analysis Research Core,
# Center for Human Movement Variability, University of Nebraska at Omaha.
# MIT licence. See LICENSE.txt.

import numpy as np

from ._rank_of import rank_of


def surr_theiler(x, algorithm, rng=None):
    """Surrogate data (Theiler et al. 1992), as MATLAB's surr_theiler.

    Algorithm 0 shuffles the series (keeps the distribution), 1 randomises
    the Fourier phases (keeps the amplitude spectrum), 2 is the amplitude
    adjusted Fourier transform (keeps the distribution exactly and the
    spectrum approximately). x is a vector or a matrix with one series per
    column; rng is a numpy Generator or seed. Returns a surrogate of the
    same shape as x.
    """
    rng = np.random.default_rng(rng)
    x = np.asarray(x, dtype=np.float64)
    is_vec = x.ndim == 1
    X = x[:, None] if is_vec else x
    if algorithm == 0:
        z = np.column_stack([v[np.argsort(rng.standard_normal(v.size), kind="stable")]
                             for v in X.T])
    elif algorithm in (1, 2):
        z = _surr1(X, algorithm, rng)
    else:
        raise ValueError(f"Unknown algorithm {algorithm}. Supported: 0, 1, 2.")
    z = z.reshape(X.shape)
    return z.ravel() if is_vec else z


def _surr1(x, algorithm, rng):
    r, c = x.shape
    if algorithm == 2:
        # Gaussianise: each value becomes the same-rank value of a sorted
        # normal sample.
        sr = np.sort(rng.standard_normal((c, r)).T, axis=0)
        sx = np.sort(x, axis=0)
        y = np.column_stack([sr[rank_of(x[:, k]) - 1, k] for k in range(c)])
    else:
        y = x
    m = y.mean(axis=0)
    y = y - m
    fy = np.fft.fft(y, axis=0)
    nHalf = (r - 1) // 2
    rot = np.ones(r, dtype=complex)
    if nHalf > 0:
        ph = 2 * np.pi * rng.random(nHalf)
        rot[1 : nHalf + 1] = np.exp(1j * ph)
        rot[r - 1 : r - nHalf - 1 : -1] = np.conj(np.exp(1j * ph))
    yy = np.real(np.fft.ifft(fy * rot[:, None], axis=0))
    yy = yy + m
    if algorithm == 2:
        yy = np.column_stack([sx[rank_of(yy[:, k]) - 1, k] for k in range(c)])
    return yy
