"""Largest Lyapunov exponent, Wolf's method."""

# Copyright (c) 2021-2026 Quantitative Analysis Research Core,
# Center for Human Movement Variability, University of Nebraska at Omaha.
# MIT licence. See LICENSE.txt.

import numpy as np

from .psr import psr


def _min_index(v):
    """1-based index of the smallest non-NaN entry (first on ties), 1 if none."""
    if np.all(np.isnan(v)):
        return 1
    return int(np.nanargmin(v)) + 1


def lye_w(x, fs, delay, dim, evolve, SCALEMX=None, ANGLMX=None, ZMULT=None):
    """Largest Lyapunov exponent (Wolf, Swift, Swinney & Vastano 1985), as
    MATLAB's lye_w.

    x is a vector, embedded with psr(x, delay, dim), or a phase space with
    one coordinate per column. A reference trajectory and a neighbour are
    followed for evolve samples, the log stretch recorded, and the neighbour
    replaced by one near the evolved point and close in direction.
    SCALEMX, ANGLMX and ZMULT (default a tenth of the data range, 30 degrees
    in radians, 1) must be given all together or not at all.

    Returns (out, LyE): out has one row per replacement step [iteration,
    point, neighbour, start distance, end distance, running exponent,
    current scale, best angle, maximum angle] with point and neighbour
    0-based and angles in degrees (rows never reached stay 0, as in MATLAB);
    LyE is the exponent in bits per unit time.
    """
    given = [a is not None for a in (SCALEMX, ANGLMX, ZMULT)]
    x = np.asarray(x, dtype=np.float64)
    if not any(given):
        SCALEMX = (x.max() - x.min()) / 10
        ANGLMX = 30 * np.pi / 180
        ZMULT = 1
    elif not all(given):
        raise ValueError("not enough input arguements")
    DT = 1 / fs
    if x.ndim == 1:
        x = x[:, None]
    if x.shape[1] == 1:
        Y = psr(x, delay, dim)
        NPT = x.shape[0] - (dim - 1) * delay - evolve
        Y = Y[: NPT + evolve]
    else:
        Y = x
        NPT = Y.shape[0] - evolve
    rows = [np.zeros(9) for _ in range(NPT // evolve)]
    ITS = 0
    distSUM = 0.0
    thbest = 0.0
    OUTMX = SCALEMX
    pair = None                                   # 1-based, as MATLAB
    LyE = np.nan
    for i in range(1, NPT + 1, evolve):
        cp = i                                    # 1-based
        if cp == 1:
            d = np.sqrt(np.sum((Y[:NPT] - Y[cp - 1]) ** 2, axis=1))
            d[d <= 0] = np.nan
            ex = np.arange(cp - 10, cp + 11)
            ex = ex[(ex >= 1) & (ex <= NPT)]
            d[ex - 1] = np.nan
            pair = _min_index(d)
        start_dist = np.sqrt(np.sum((Y[cp - 1] - Y[pair - 1]) ** 2))
        end_dist = np.sqrt(np.sum((Y[cp + evolve - 1] - Y[pair + evolve - 1]) ** 2))
        distSUM = distSUM + np.log2(end_dist / start_dist) / (evolve * DT)
        ITS += 1
        LyE = distSUM / ITS
        row = i // evolve + 1
        while row > len(rows):
            rows.append(np.zeros(9))
        rows[row - 1] = np.array([ITS, cp - 1, pair - 1, start_dist, end_dist, LyE, OUTMX,
                                  thbest * 180 / np.pi, ANGLMX * 180 / np.pi])
        ZMULT = 1
        if end_dist < SCALEMX:
            pair = pair + evolve
            if pair > NPT:
                pair = pair - evolve
                pair, ANGLMX, thbest, OUTMX = _next_point(1, Y, cp, pair, NPT, evolve,
                                                          SCALEMX, ANGLMX)
            continue
        pair, ANGLMX, thbest, OUTMX = _next_point(0, Y, cp, pair, NPT, evolve, SCALEMX, ANGLMX)
    out = np.vstack(rows) if rows else np.zeros((0, 9))
    return out, LyE


def _next_point(flag, Y, cp, pair, NPT, evolve, SCALEMX, ANGLMX):
    """Replacement neighbour near the evolved point, preferring a small angle
    to the current separation (MATLAB's GetNextPoint). Indices 1-based."""
    p1 = Y[cp + evolve - 1]
    Vnew = p1 - Y[:NPT]
    d = np.sqrt(np.sum(Vnew ** 2, axis=1))
    ex = np.arange(cp + evolve - 10, cp + evolve + 11)
    ex = ex[(ex >= 1) & (ex <= NPT)]
    d[ex - 1] = np.nan
    Vcurr = p1 - Y[pair + evolve - 1]
    end_dist = np.sqrt(np.sum(Vcurr ** 2))
    with np.errstate(divide="ignore", invalid="ignore"):
        cosTheta = np.abs((Vnew @ Vcurr) / (d * end_dist))
        # MATLAB's acos of a cosine rounded just above 1 is complex:
        # comparisons use its real part (0) and abs() its modulus (acosh).
        above = ~np.isnan(cosTheta) & (cosTheta > 1)
        thetaRe = np.where(above, 0.0, np.arccos(np.minimum(cosTheta, 1)))
        thetaAbs = np.where(above, np.arccosh(np.maximum(cosTheta, 1)), thetaRe)
        if flag == 0:
            pot = d.copy()
            pot[~np.isnan(d) & ((d <= 0) | (thetaRe >= ANGLMX))] = np.nan
            np_ = _min_index(pot)
            if not np.isnan(pot[np_ - 1]) and pot[np_ - 1] <= SCALEMX:
                return np_, 30 * np.pi / 180, abs(thetaAbs[np_ - 1]), SCALEMX
        pot = d.copy()
        pot[~np.isnan(d) & (d <= 0)] = np.nan
    return _min_index(pot), ANGLMX, ANGLMX, SCALEMX
