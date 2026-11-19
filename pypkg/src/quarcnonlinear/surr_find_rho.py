"""Noise radius for pseudo-periodic surrogates."""

# Copyright (c) 2021-2026 Quantitative Analysis Research Core,
# Center for Human Movement Variability, University of Nebraska at Omaha.
# MIT licence. See LICENSE.txt.

import warnings

import numpy as np

from .surr_pseudo_periodic import surr_pseudo_periodic


def _di(yi, n=2):
    """Number of runs copied from the original longer than n samples."""
    return int(np.sum(np.diff(np.flatnonzero(np.diff(yi) != 1)) > n))


def surr_find_rho(x, delay, dim, rhoL=0.1, rhoH=1, rng=None):
    """Noise radius rho for surr_pseudo_periodic, as MATLAB's surr_find_rho.

    Bisects between rhoL and rhoH for the rho that maximises the number of
    short segments copied from the original series (Small et al. 2001).
    rng is a numpy Generator or seed. Returns (rho, out): out has one row per
    trial with columns [trial, rho, di].
    """
    if np.size(rhoL) != 1 or np.size(rhoH) != 1 or not (rhoL > 0) or not (rhoH > rhoL):
        raise ValueError("Search bounds must satisfy 0 < rhoL < rhoH.")
    rng = np.random.default_rng(rng)

    def di_at(r):
        return _di(surr_pseudo_periodic(x, delay, dim, r, rng=rng)[1])

    diH = di_at(rhoH)
    diL = di_at(rhoL)
    out = [[1, rhoH, diH], [2, rhoL, diL]]
    if diH > diL:
        dmax, rho = diH, rhoH
    else:
        dmax, rho = diL, rhoL
    ind = 3
    while abs(rhoH - rhoL) / rhoL > 0.02:
        rhoi = (rhoH + rhoL) / 2
        di = di_at(rhoi)
        out.append([ind, rhoi, di])
        ind += 1
        if di > dmax:
            dmax, rho = di, rhoi
        if diL < diH:
            diL, rhoL = di, rhoi
        else:
            diH, rhoH = di, rhoi
    out = np.array(out, dtype=np.float64)
    lo, hi = out[1, 1], out[0, 1]
    margin = 0.15 * np.log(hi / lo)
    if np.log(rho / lo) < margin or np.log(hi / rho) < margin:
        warnings.warn(f"The best rho found ({rho:g}) lies against a search bound [{lo:g} {hi:g}], "
                      "so the optimum is probably not bracketed and this value may not be a "
                      "maximum of di(rho). Re-run with wider bounds, e.g. "
                      f"surr_find_rho(y, delay, dim, {lo / 10:g}, {hi * 2:g}).")
    return rho, out
