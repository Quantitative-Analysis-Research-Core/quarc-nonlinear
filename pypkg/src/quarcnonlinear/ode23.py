"""Bogacki-Shampine 3(2) solver with MATLAB's step control."""

# Copyright (c) 2021-2026 Quantitative Analysis Research Core,
# Center for Human Movement Variability, University of Nebraska at Omaha.
# MIT licence. See LICENSE.txt.

import numpy as np

from ._rk_solve import rk_solve

_TAB = {
    "pow": 1 / 3, "refine": 1,
    "A": np.array([1 / 2, 3 / 4, 1]),
    "B": np.array([[1 / 2, 0, 2 / 9],
                   [0, 3 / 4, 1 / 3],
                   [0, 0, 4 / 9],
                   [0, 0, 0]]),
    "E": np.array([-5 / 72, 1 / 12, 1 / 9, -1 / 8]),
    "BI": np.array([[1, -4 / 3, 5 / 9],
                    [0, 1, -2 / 3],
                    [0, 4 / 3, -8 / 9],
                    [0, -1, 1]]),
}


def ode23(f, tspan, y0, rtol=1e-3, atol=1e-6):
    """Solve y' = f(t, y) with the Bogacki-Shampine 3(2) pair, following
    MATLAB's ode23 step control. Same arguments as ode45; with two times in
    tspan every step is reported. Returns (t, y), one row of y per time.
    """
    return rk_solve(f, tspan, y0, rtol, atol, _TAB)
