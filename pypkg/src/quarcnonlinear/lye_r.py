"""Largest Lyapunov exponent, Rosenstein's method."""

# Copyright (c) 2021-2026 Quantitative Analysis Research Core,
# Center for Human Movement Variability, University of Nebraska at Omaha.
# MIT licence. See LICENSE.txt.

import numpy as np

from . import _core
from .psr import psr


def _round(v):
    """MATLAB's round: halves away from zero (Python's round goes to even)."""
    return int(np.sign(v) * np.floor(abs(v) + 0.5))


def lye_r(x, fs, delay, dim, theiler=None, slope=None, mean_period=None):
    """Largest Lyapunov exponent (Rosenstein, Collins & De Luca 1993), as MATLAB's lye_r.

    x is a vector, embedded with psr(x, delay, dim), or a phase space with
    one coordinate per column. Each point's nearest neighbour outside a
    Theiler window (default round(0.8*delay) samples) is followed forward,
    and the mean log divergence of the pairs is taken at each step.

    Returns out, an M-by-3 array of [point, neighbour, mean log divergence]
    (indices 0-based). With slope = [s0, s1, l0, l1] (in mean periods) and
    mean_period (seconds), returns (LyES, LyEL, out): the short- and
    long-range slopes of the divergence curve; NaN where a range runs past
    the curve.
    """
    if theiler is None:
        theiler = _round(delay * 0.8)
    if np.ndim(theiler) != 0 or not np.isfinite(theiler) or theiler < 0:
        raise ValueError("TheilerWindow must be a non-negative finite scalar, in samples.")
    theiler = _round(theiler)
    x = np.asarray(x, dtype=np.float64)
    if x.ndim == 1:
        x = x[:, None]
    if x.shape[1] > x.shape[0]:
        x = x.T
    Y = x if x.shape[1] > 1 else psr(x, delay, dim)
    M, D = Y.shape
    nn, avg, count = _core.lye_r_divergence(np.ascontiguousarray(Y.ravel(order="F")), M, D, theiler)
    out = np.column_stack([np.arange(M), nn, avg])
    if slope is None:
        return out

    ave = np.asarray(avg[:count])
    time = np.arange(count) / fs / mean_period
    short = (1 if slope[0] == 0 else _round(slope[0] * mean_period * fs),
             _round(slope[1] * mean_period * fs))
    long_ = (_round(slope[2] * mean_period * fs), _round(slope[3] * mean_period * fs))

    def fit(r):                              # r is MATLAB's 1-based, inclusive range
        if r[1] > count:
            return np.nan
        t, y = time[r[0] - 1 : r[1]], ave[r[0] - 1 : r[1]]
        return np.linalg.lstsq(np.column_stack([t, np.ones(t.size)]), y, rcond=None)[0][0]

    return fit(short), fit(long_), out
