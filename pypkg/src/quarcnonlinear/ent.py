"""Entropy of a time series, by algorithm."""

# Copyright (c) 2021-2026 Quantitative Analysis Research Core,
# Center for Human Movement Variability, University of Nebraska at Omaha.
# MIT licence. See LICENSE.txt.

import numpy as np

from .ent_ap import ent_ap
from .ent_ms_plus import ent_ms_plus
from .ent_permu import ent_permu
from .ent_samp import ent_samp
from .ent_symbolic import ent_symbolic
from .ent_xap import ent_xap
from .ent_xsamp import ent_xsamp


def ent(x, algorithm="sample", y=None, dim=2, radius=0.2, delay=1, scale=20):
    """One entry point for the entropy family, as MATLAB's ent.

    algorithm is "sample" (default), "approximate", "permutation",
    "symbolic" (binary series only) or "multiscale" (refined composite
    multiscale sample entropy over scales 1..scale). Supplying y gives the
    cross entropy, with "sample" or "approximate" only.

    Returns (value, extra): the entropy (an array of length scale for
    "multiscale") and a dict with the algorithm, the arguments used and any
    secondary outputs (patternCounts; rcmse, cmse, mse, msfe, gmse).
    """
    x = np.asarray(x, dtype=np.float64).ravel()
    if np.any(np.isnan(x)):
        raise ValueError("x contains NaN. Entropy is undefined for missing data; "
                         "remove or impute first.")
    if x.size <= dim + 1:
        raise ValueError(f"Series of {x.size} samples is too short for dim = {int(dim)}.")
    isCross = y is not None
    if isCross:
        y = np.asarray(y, dtype=np.float64).ravel()
        if np.any(np.isnan(y)):
            raise ValueError("y contains NaN.")
        if y.size != x.size:
            raise ValueError(f"x has {x.size} samples and y has {y.size}; "
                             "a cross entropy needs equal lengths.")
    algorithm = str(algorithm).lower()

    def reject_cross(name):
        if isCross:
            raise ValueError(f'algorithm "{name}" has no cross-entropy form. Supply y only '
                             'with "sample" or "approximate".')

    extra = dict(algorithm=algorithm, dim=dim, cross=isCross)
    if algorithm == "sample":
        extra["radius"] = radius
        extra["estimator"] = "cross sample" if isCross else "sample"
        value = ent_xsamp(x, y, dim, radius, 1) if isCross else ent_samp(x, dim, radius)
    elif algorithm == "approximate":
        extra["radius"] = radius
        extra["estimator"] = "cross approximate" if isCross else "approximate"
        value = ent_xap(x, y, dim, radius, 1) if isCross else ent_ap(x, dim, radius)
    elif algorithm == "permutation":
        reject_cross("permutation")
        value, counts = ent_permu(x, dim, delay)
        extra["delay"] = delay
        extra["patternCounts"] = counts
        extra["estimator"] = "permutation"
    elif algorithm == "symbolic":
        reject_cross("symbolic")
        u = np.unique(x)
        if not np.all(np.isin(u, [0, 1])):
            raise ValueError(
                f'algorithm "symbolic" needs a binary series; x has {u.size} distinct '
                f"values in [{u.min():g}, {u.max():g}]. Threshold it first, for example "
                "x = (x > np.median(x)).astype(float).")
        extra["estimator"] = "symbolic"
        value = ent_symbolic(x, dim)
    elif algorithm == "multiscale":
        reject_cross("multiscale")
        rcmse, cmse, mse, msfe, gmse = ent_ms_plus(x, scale, dim, radius)
        extra.update(radius=radius, scale=scale, rcmse=rcmse, cmse=cmse, mse=mse,
                     msfe=msfe, gmse=gmse,
                     estimator="refined composite multiscale sample")
        value = rcmse
    else:
        raise ValueError(f'Unknown algorithm "{algorithm}". Supported: "sample", '
                         '"approximate", "permutation", "symbolic", "multiscale".')
    return value, extra
