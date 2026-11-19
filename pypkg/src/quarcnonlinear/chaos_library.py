"""Trajectories of standard chaotic systems."""

# Copyright (c) 2021-2026 Quantitative Analysis Research Core,
# Center for Human Movement Variability, University of Nebraska at Omaha.
# MIT licence. See LICENSE.txt.

import math

import numpy as np

from ._mround import mround
from .ode23 import ode23
from .ode45 import ode45


def chaos_library(S, t, IC, p):
    """Trajectory of a standard chaotic system, as MATLAB's chaos_library.

    S is "Lorenz", "Rossler", "Aizawa" (ode23), "DequanLi" (ode23),
    "NoseHoover", "QiChen", "YuWang", "LuChen", "Arneodo", "TSUCSI", or the
    maps "Hennon" and "Logistic". For a flow t holds the output times; for a
    map, the number of iterates or a vector whose length sets it. IC is the
    initial condition and p the parameters, in the MATLAB order.

    Returns (t, y), one row of y per time and one column per variable. Flows
    agree with MATLAB to rounding over a short horizon, then drift (chaos);
    maps agree exactly.
    """
    p = np.atleast_1d(np.asarray(p, dtype=np.float64)).tolist()
    IC = np.atleast_1d(np.asarray(IC, dtype=np.float64))
    flows = {"Rossler": (_rossler, ode45), "Lorenz": (_lorenz, ode45),
             "Aizawa": (_aizawa, ode23), "DequanLi": (_dequanli, ode23),
             "NoseHoover": (_nosehoover, ode45), "QiChen": (_qichen, ode45),
             "YuWang": (_yuwang, ode45), "LuChen": (_luchen, ode45),
             "Arneodo": (_arneodo, ode45), "TSUCSI": (_tsucsi, ode45)}
    if S in flows:
        rhs, solver = flows[S]
        return solver(lambda tt, y: rhs(tt, y, p), t, IC)
    if S == "Hennon":
        n = _map_length(t)
        x = np.zeros(n)
        z = np.zeros(n)
        x[0] = IC[0]
        z[0] = IC[1]
        for c in range(n - 1):
            x[c + 1] = 1 - p[0] * x[c] ** 2 + z[c]
            z[c + 1] = p[1] * x[c]
        return _map_time(t, n), np.column_stack([x, z])
    if S == "Logistic":
        n = _map_length(t)
        x = np.zeros(n)
        x[0] = IC[0]
        for c in range(n - 1):
            x[c + 1] = p[0] * x[c] * (1 - x[c])
        return _map_time(t, n), x[:, None]
    raise ValueError(f"Unknown system '{S}'.")


def _map_length(t):
    t = np.atleast_1d(t)
    n = int(mround(t[0])) if t.size == 1 else t.size
    if n < 1:
        raise ValueError("t must ask for at least one data point.")
    return n


def _map_time(t, n):
    t = np.atleast_1d(np.asarray(t, dtype=np.float64))
    return np.arange(n, dtype=np.float64) if t.size == 1 else t


def _rossler(t, y, p):
    a, b, c = y
    return np.array([-b - c, a + p[0] * b, p[1] + c * (a - p[2])])


def _lorenz(t, y, p):
    a, b, c = y
    return np.array([p[0] * (b - a), -a * c + p[1] * a - b, a * b - p[2] * c])


def _dequanli(t, y, p):
    a, b, c = y
    return np.array([p[0] * (b - a) + p[1] * (a * c), p[2] * a + p[3] * b - a * c,
                     p[4] * c + a * b - p[5] * (a * a)])


def _nosehoover(t, y, p):
    a, b, c = y
    return np.array([b, -a + b * c, p[0] - b * b])


def _qichen(t, y, p):
    a, b, c = y
    return np.array([p[0] * (b - a) + b * c, p[2] * a + b - a * c, a * b - p[1] * c])


def _yuwang(t, y, p):
    a, b, c = y
    return np.array([p[0] * (b - a), p[1] * a - p[2] * a * c, math.exp(a * b) - p[3] * c])


def _luchen(t, y, p):
    a, b, c = y
    return np.array([-1 * (p[0] * p[1] * a / (p[0] + p[1])) - b * c + p[2],
                     p[0] * b + a * c, p[1] * c + a * b])


def _arneodo(t, y, p):
    a, b, c = y
    return np.array([b, c, -p[0] * a - p[1] * b - c + p[2] * a ** 3])


def _tsucsi(t, y, p):
    a, b, c = y
    return np.array([p[0] * (b - a) + p[2] * a * c, p[4] * b - a * c,
                     p[1] * c + a * b - p[3] * a * a])


def _aizawa(t, y, p):
    a, b, c = y
    return np.array([(c - p[1]) * a - p[3] * b, p[3] * a + (c - p[1]) * b,
                     p[2] + p[0] * c - c ** 3 / 3 - a ** 2 + p[5] * c * a ** 3])
