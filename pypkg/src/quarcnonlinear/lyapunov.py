"""Largest Lyapunov exponent, one entry point for both methods."""

# Copyright (c) 2021-2026 Quantitative Analysis Research Core,
# Center for Human Movement Variability, University of Nebraska at Omaha.
# MIT licence. See LICENSE.txt.

import warnings

import numpy as np

from ._mround import mround
from .lye_r import lye_r
from .lye_w import lye_w
from .psr import psr


def lyapunov(x, fs, algorithm="rosenstein", delay=1, dim=3, theiler=None, evolve=None):
    """Largest Lyapunov exponent in nats per unit time, as MATLAB's lyapunov.

    x is a vector (embedded with delay and dim) or a phase space with one
    coordinate per column. algorithm is "rosenstein" (alias "r", default)
    or "wolf" (alias "w"). theiler is the Theiler window for "rosenstein"
    (default round(0.8*delay)); evolve, the propagation length in samples,
    is required for "wolf" (1 for a map, about 10 for a flow at ~40 samples
    per orbit).

    Returns (lambda, extra). extra holds phaseSpace, algorithm, units, fs,
    theilerWindow, tau, dim and, for "rosenstein", divergence, matchedPairs
    (0-based), scalingRegion (0-based indices into divergence), fitR2 and
    estimator; for "wolf", iterations (lye_w's out), bits, estimator and
    evolve.
    """
    x = np.asarray(x, dtype=np.float64)
    if np.isnan(x).any():
        raise ValueError("x contains NA. Remove or impute before estimating an exponent.")
    if x.ndim == 1 or x.shape[1] == 1:
        x = x.ravel()
        if x.size <= (dim - 1) * delay + 2:
            raise ValueError(f"Series of {x.size} samples is too short for dim = {int(dim)} "
                             f"at delay = {int(delay)}.")
        Y = psr(x, delay, dim)
    else:
        Y = x
        if theiler is None and delay == 1:
            warnings.warn("A phase space was supplied without theiler or delay, so the temporal "
                          "exclusion defaults to 1 sample. That is almost certainly too small "
                          "for an embedded series and will bias the exponent. Pass theiler (or "
                          "delay) to match the delay used to build the phase space.")
    if theiler is None:
        theiler = int(mround(delay * 0.8))
    algorithm = str(algorithm).lower()
    extra = {"phaseSpace": Y, "algorithm": algorithm, "units": "nats per unit time", "fs": fs,
             "theilerWindow": theiler, "tau": delay, "dim": Y.shape[1]}
    if algorithm in ("rosenstein", "r"):
        out = lye_r(Y, fs, delay, Y.shape[1], theiler=theiler)
        d = out[:, 2]
        slope, idx, r2 = _scaling_region(d, fs)
        extra["divergence"] = d
        extra["matchedPairs"] = out[:, :2]
        extra["scalingRegion"] = idx - 1
        extra["fitR2"] = r2
        extra["estimator"] = "rosenstein"
        lam = slope
    elif algorithm in ("wolf", "w"):
        if evolve is None:
            raise ValueError('algorithm="wolf" needs evolve, the propagation length in samples. '
                             "There is no default: use evolve = 1 for a map, and about 10 for a "
                             "flow sampled at ~40 samples per orbit.")
        out, L = lye_w(Y, fs, delay, Y.shape[1], evolve)
        lam = L * np.log(2)
        extra["iterations"] = out
        extra["bits"] = L
        extra["estimator"] = "wolf"
        extra["evolve"] = evolve
    else:
        raise ValueError(f'Unknown algorithm "{algorithm}". Supported: "rosenstein" (alias "r"), '
                         '"wolf" (alias "w").')
    return lam, extra


def _fit(xv, yv):
    """Least-squares line yv = a*xv + b; returns (a, b)."""
    A = np.column_stack([xv, np.ones(xv.size)])
    return np.linalg.lstsq(A, yv, rcond=None)[0]


def _line_r2(xv, yv):
    if xv.size < 2 or np.all(yv == yv[0]):
        return -np.inf
    ss = np.sum((yv - yv.mean()) ** 2)
    if ss <= 0:
        return -np.inf
    a, b = _fit(xv, yv)
    return 1 - np.sum((yv - (a * xv + b)) ** 2) / ss


def _scaling_region(d, fs):
    """Longest near-linear window of the divergence curve before it
    saturates. Returns (slope, idx, r2) with idx 1-based."""
    d = d[np.isfinite(d)]
    n = d.size
    empty = np.zeros(0, dtype=np.int64)
    if n < 6:
        return np.nan, empty, np.nan
    hi = max(6, min(_saturation_index(d), 200))
    best = {"r2": -np.inf, "a": 1, "b": min(hi, 5), "len": 0}
    cand = []
    for a in range(1, min(4, max(1, hi - 4)) + 1):
        if a + 3 > hi:
            continue
        for b in range(a + 3, hi + 1):
            xv = np.arange(a, b + 1, dtype=np.float64)
            c = {"r2": _line_r2(xv, d[a - 1 : b]), "a": a, "b": b, "len": b - a + 1}
            cand.append(c)
            if c["r2"] > best["r2"]:
                best = c
    if not cand:
        return np.nan, empty, np.nan
    pick = best
    for c in cand:
        if c["r2"] >= best["r2"] - 0.02 and c["len"] > pick["len"]:
            pick = c
    idx = np.arange(pick["a"], pick["b"] + 1)
    slope = _fit(idx / fs, d[idx - 1])[0]
    return slope, idx, pick["r2"]


def _saturation_index(d):
    """1-based index where the divergence curve stops rising."""
    n = d.size
    span = d.max() - d.min()
    if span <= 0:
        return n
    tol = 0.01 * span
    look = int(max(3, mround(0.02 * n)))
    runmax = d[0]
    for i in range(2, n - look + 1):
        runmax = max(runmax, d[i - 1])
        if d[i : i + look].max() <= runmax + tol:
            return i
    return n
