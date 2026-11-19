"""Shared last step of the RQA functions: radius, line measures, plot."""

# Copyright (c) 2021-2026 Quantitative Analysis Research Core,
# Center for Human Movement Variability, University of Nebraska at Omaha.
# MIT licence. See LICENSE.txt.

import numpy as np

from .ent_weighted import ent_weighted
from .line_hist import line_hist
from .set_radius import set_radius


def _line_entropy(h):
    """Shannon entropy (bits) of the distribution of line lengths h."""
    if h.size == 0:
        return 0.0
    p = np.bincount(h.astype(np.int64))[1:].astype(np.float64)
    p = p[p > 0] / p.sum()
    return -np.sum(p * np.log2(p))


def rqa_finish(x, a, wrp, param, threshold, type, iter, theiler, dmin, vmin,
               RESULTS, withEntrW):
    """Pick the radius, count lines and fill RESULTS, as R's rqa_finish.

    Returns a dict with RP (the plot rotated as MATLAB's imrotate(1 - A, 90),
    so 0 marks a recurrence) and RESULTS.
    """
    if param == "rad":
        radius = threshold
        lh = line_hist(x, a, threshold, type, theiler)
        recurrence, dh, vh, A = lh["recurrence"], lh["diag_hist"], lh["vertical_hist"], lh["a"]
    elif param == "rec":
        sr = set_radius(x, a, 0.01, 0.5, threshold, type, iter, theiler)
        recurrence, dh, vh, A = sr["rec"], sr["diag_hist"], sr["vertical_hist"], sr["A"]
        radius = sr["rad_final"]
    else:
        raise ValueError('param must be "rad" or "rec".')
    RESULTS = dict(RESULTS)
    RESULTS["RADIUS"] = radius
    RESULTS["Size"] = A.shape[0]
    RESULTS["REC"] = recurrence
    if recurrence > 0:
        dh = np.asarray(dh, dtype=np.float64)
        vh = np.asarray(vh, dtype=np.float64)
        d = dh[dh >= dmin]
        v = vh[vh >= vmin]
        with np.errstate(invalid="ignore", divide="ignore"):
            RESULTS["DET"] = 100 * d.sum() / dh.sum()
            RESULTS["MeanL"] = d.mean() if d.size else np.nan
            RESULTS["MaxL"] = d.max() if d.size else -np.inf
            RESULTS["EntrL"] = _line_entropy(d)
            RESULTS["LAM"] = 100 * v.sum() / vh.sum()
            RESULTS["MeanV"] = v.mean() if v.size else np.nan
            RESULTS["MaxV"] = v.max() if v.size else -np.inf
            RESULTS["EntrV"] = _line_entropy(v)
        if withEntrW:
            RESULTS["EntrW"] = ent_weighted(wrp)
    else:
        for f in ("DET", "MeanL", "MaxL", "EntrL", "LAM", "MeanV", "MaxV", "EntrV"):
            RESULTS[f] = np.nan
        if withEntrW:
            RESULTS["EntrW"] = np.nan
    return dict(RP=np.rot90(1 - A), RESULTS=RESULTS)
