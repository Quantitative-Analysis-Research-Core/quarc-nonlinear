"""Radius giving a target percent recurrence."""

# Copyright (c) 2021-2026 Quantitative Analysis Research Core,
# Center for Human Movement Variability, University of Nebraska at Omaha.
# MIT licence. See LICENSE.txt.

import warnings

from .line_hist import line_hist


def set_radius(x, a, radius_start, radius_end, threshold, type, iter, theiler=0):
    """Radius that gives threshold percent recurrence, as MATLAB's setRadius.

    Brackets the target from both sides (radius_start below, radius_end
    above), then bisects iter times. Arguments x, a, type and theiler are as
    in line_hist. Returns a dict with rec, diag_hist, vertical_hist,
    rad_final and A, all from the last bisection step. Reports the result
    with warnings.warn (R's message()).
    """
    maxAdjust = 200
    rec = line_hist(x, a, radius_start, type, theiler)["recurrence"]
    k = 0
    while rec == 0 or rec > threshold:
        k += 1
        if k > maxAdjust:
            raise ValueError(
                f"Could not bracket %REC = {threshold:g} from below after {maxAdjust} adjustments "
                f"(radius {radius_start:g} gives %REC = {rec:g}). The distance matrix takes too few "
                "distinct values for this target to be reachable, which is usual "
                'for binary or categorical data. Set the radius directly with '
                'param = "rad", or pick a target that lies on an achievable step.')
        radius_start = radius_start * 2 if rec == 0 else radius_start / 1.5
        rec = line_hist(x, a, radius_start, type, theiler)["recurrence"]
    rec = line_hist(x, a, radius_end, type, theiler)["recurrence"]
    k = 0
    while rec < threshold:
        k += 1
        if k > maxAdjust:
            raise ValueError(
                f"Could not bracket %REC = {threshold:g} from above after {maxAdjust} adjustments "
                f"(radius {radius_end:g} gives %REC = {rec:g}).")
        radius_end = radius_end * 2
        rec = line_hist(x, a, radius_end, type, theiler)["recurrence"]
    lv, hv = radius_start, radius_end
    lh, mid = None, None
    for _ in range(int(iter)):
        mid = (lv + hv) / 2
        lh = line_hist(x, a, mid, type, theiler)
        if lh["recurrence"] < threshold:
            lv = mid
        else:
            hv = mid
    if lh is None:
        raise ValueError("iter must be at least 1.")
    warnings.warn(f"% recurrence = {lh['recurrence']:.5g}, radius = {mid:.5g}", stacklevel=2)
    return dict(rec=lh["recurrence"], diag_hist=lh["diag_hist"],
                vertical_hist=lh["vertical_hist"], rad_final=mid, A=lh["a"])
