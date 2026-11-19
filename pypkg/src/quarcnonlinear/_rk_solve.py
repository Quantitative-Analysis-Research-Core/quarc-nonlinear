"""Explicit Runge-Kutta stepping with MATLAB's error and step control."""

# Copyright (c) 2021-2026 Quantitative Analysis Research Core,
# Center for Human Movement Variability, University of Nebraska at Omaha.
# MIT licence. See LICENSE.txt.

import numpy as np


def rk_solve(f, tspan, y0, rtol, atol, tab):
    """Integrate y' = f(t, y) with the Butcher tableau tab (pow, refine, A,
    B, E, BI), as MATLAB's ode45/ode23 do (Shampine & Reichelt 1997).

    Returns (t, y): output times and one row of y per time. With two times
    in tspan, every step is reported plus refine-1 interpolated points;
    otherwise the solution at the times of tspan.
    """
    tspan = np.asarray(tspan, dtype=np.float64).ravel()
    y = np.asarray(y0, dtype=np.float64).ravel().copy()

    def F(t, v):
        return np.asarray(f(t, v), dtype=np.float64).ravel()

    A, B, E, BI = tab["A"], tab["B"], tab["E"], tab["BI"]
    t0 = tspan[0]
    tfinal = tspan[-1]
    tdir = np.sign(tfinal - t0)
    rtol = max(rtol, 100 * np.finfo(float).eps)
    threshold = atol / rtol
    pow_ = tab["pow"]
    nst = A.size + 1
    refine_mode = tspan.size == 2
    hmax = abs(0.1 * (tfinal - t0))
    t = t0
    f1 = F(t, y)
    K = np.zeros((y.size, nst))
    K[:, 0] = f1
    htspan = abs(tspan[1] - tspan[0])
    absh = min(hmax, htspan)
    rh = np.max(np.abs(f1 / np.maximum(np.abs(y), threshold))) / (0.8 * rtol ** pow_)
    if absh * rh > 1:
        absh = 1 / rh
    absh = max(absh, 16 * np.spacing(abs(t)))

    tout = [t0]
    yout = [y]
    next_idx = 1
    done = False
    while not done:
        hmin = 16 * np.spacing(abs(t))
        absh = min(hmax, max(hmin, absh))
        h = tdir * absh
        if 1.1 * absh >= abs(tfinal - t):
            h = tfinal - t
            absh = abs(h)
            done = True
        nofailed = True
        while True:
            hB = h * B
            for j in range(1, nst - 1):
                K[:, j] = F(t + h * A[j - 1], y + K[:, :j] @ hB[:j, j - 1])
            tnew = t + h * A[nst - 2]
            if done:
                tnew = tfinal
            h = tnew - t
            ynew = y + K[:, : nst - 1] @ hB[: nst - 1, nst - 2]
            K[:, nst - 1] = F(tnew, ynew)
            err = absh * np.max(np.abs(K @ E) / np.maximum(np.maximum(np.abs(y), np.abs(ynew)),
                                                         threshold))
            if err > rtol:
                if absh <= hmin:
                    raise ValueError(f"Unable to meet integration tolerances at t = {t:g}.")
                if nofailed:
                    nofailed = False
                    absh = max(hmin, absh * max(0.1 if pow_ == 1 / 5 else 0.5,
                                                0.8 * (rtol / err) ** pow_))
                else:
                    absh = max(hmin, 0.5 * absh)
                h = tdir * absh
                done = False
            else:
                break

        def interp(ti):
            s = (ti - t) / h
            return y + (K @ (h * BI)) @ np.cumprod(np.full(BI.shape[1], s))

        if refine_mode:
            r = tab["refine"]
            for ti in t + h * (np.arange(1, r) / r):
                tout.append(ti)
                yout.append(interp(ti))
            tout.append(tnew)
            yout.append(ynew)
        else:
            while next_idx < tspan.size and tdir * (tnew - tspan[next_idx]) >= 0:
                ti = tspan[next_idx]
                tout.append(ti)
                yout.append(ynew if ti == tnew else interp(ti))
                next_idx += 1
        if nofailed:
            temp = 1.25 * (err / rtol) ** pow_
            absh = absh / temp if temp > 0.2 else 5 * absh
        t = tnew
        y = ynew
        K[:, 0] = K[:, nst - 1]
    return np.array(tout), np.vstack(yout)
