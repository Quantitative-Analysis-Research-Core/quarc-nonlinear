"""Multidimensional recurrence quantification analysis."""

# Copyright (c) 2021-2026 Quantitative Analysis Research Core,
# Center for Human Movement Variability, University of Nebraska at Omaha.
# MIT licence. See LICENSE.txt.

import warnings

import numpy as np

from ._pairwise_dist import pairwise_dist
from ._rqa_finish import rqa_finish
from ._rqa_normalise import rqa_normalise
from ._zscore_cols import zscore_cols
from .psr import psr


def mdrqa(x, delay=1, dim=1, param="rec", threshold=2.5, zscore=1, norm="none",
          dmin=2, vmin=2, iter=20, phasespace=False):
    """Recurrence of several channels in one phase space, as MATLAB's MDRQA.

    x is a matrix with one channel per column, embedded together with delay
    and dim (dim=1: unembedded), or with phasespace=True a phase space.
    param "rec" sets the radius to reach threshold percent recurrence
    (bisecting iter times; the result is reported with warnings.warn); "rad"
    uses threshold as the radius. zscore=1 standardises the data first. norm
    is "none", "euc", "min" or "max". dmin and vmin are the shortest diagonal
    and vertical lines counted. The plotting options of MATLAB are omitted.

    Returns a dict with RP, the recurrence plot as displayed (rotated as
    MATLAB's imrotate(1 - A, 90), so 0 marks a recurrence), and RESULTS, a
    dict with DIM, EMB, DEL, NORM, ZSCORE, RADIUS, Size, REC, DET, MeanL,
    MaxL, EntrL, LAM, MeanV, MaxV, EntrV and EntrW.
    """
    x = np.asarray(x, dtype=np.float64)
    if x.ndim == 1:
        x = x[:, None]
    if zscore:
        x = zscore_cols(x)
    if phasespace:
        if norm == "euc":
            warnings.warn("PhaseSpace is true, so DATA was not rebuilt from TAU and DIM and DIM "
            "cannot be confirmed to match its width. Norm=\"euc\" scales by DIM "
            "regardless; pass DIM equal to ncol(data) or use Norm=\"none\", \"min\" or \"max\".",
                          stacklevel=2)
    elif dim > 1:
        x = psr(x, delay, dim)
    a = rqa_normalise(-pairwise_dist(x, x), norm, 1 * dim)
    return rqa_finish(x, a, a, param, threshold, "mdrqa", iter, 0, dmin, vmin,
                      dict(DIM=1, EMB=dim, DEL=delay, NORM=norm, ZSCORE=zscore), True)
