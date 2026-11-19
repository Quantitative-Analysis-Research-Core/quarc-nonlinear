"""Dormand-Prince 5(4) solver with MATLAB's step control."""

# Copyright (c) 2021-2026 Quantitative Analysis Research Core,
# Center for Human Movement Variability, University of Nebraska at Omaha.
# MIT licence. See LICENSE.txt.

import numpy as np

from ._rk_solve import rk_solve

_TAB = {
    "pow": 1 / 5, "refine": 4,
    "A": np.array([1 / 5, 3 / 10, 4 / 5, 8 / 9, 1, 1]),
    "B": np.array([[1 / 5, 3 / 40, 44 / 45, 19372 / 6561, 9017 / 3168, 35 / 384],
                   [0, 9 / 40, -56 / 15, -25360 / 2187, -355 / 33, 0],
                   [0, 0, 32 / 9, 64448 / 6561, 46732 / 5247, 500 / 1113],
                   [0, 0, 0, -212 / 729, 49 / 176, 125 / 192],
                   [0, 0, 0, 0, -5103 / 18656, -2187 / 6784],
                   [0, 0, 0, 0, 0, 11 / 84],
                   [0, 0, 0, 0, 0, 0]]),
    "E": np.array([71 / 57600, 0, -71 / 16695, 71 / 1920, -17253 / 339200, 22 / 525, -1 / 40]),
    "BI": np.array([[1, -183 / 64, 37 / 12, -145 / 128],
                    [0, 0, 0, 0],
                    [0, 1500 / 371, -1000 / 159, 1000 / 371],
                    [0, -125 / 32, 125 / 12, -375 / 64],
                    [0, 9477 / 3392, -729 / 106, 25515 / 6784],
                    [0, -11 / 7, 11 / 3, -55 / 28],
                    [0, 3 / 2, -4, 5 / 2]]),
}


def ode45(f, tspan, y0, rtol=1e-3, atol=1e-6):
    """Solve y' = f(t, y) with the Dormand-Prince 5(4) pair, following the
    error control, initial step, step update and continuous extension of
    MATLAB's ode45 (Shampine & Reichelt 1997).

    tspan holds the output times (increasing or decreasing); with two times,
    every step and three points between are reported (MATLAB's Refine 4).
    Returns (t, y), one row of y per time.
    """
    return rk_solve(f, tspan, y0, rtol, atol, _TAB)
