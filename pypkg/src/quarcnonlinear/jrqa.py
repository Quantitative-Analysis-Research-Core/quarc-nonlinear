"""Joint recurrence quantification analysis."""

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


def jrqa(x, delay=(1, 1), dim=(1, 1), param="rec", threshold=2.5, zscore=1,
         norm="none", dmin=2, vmin=2, iter=20, phasespace=False):
    """Joint recurrence of several channels, as MATLAB's JRQA.

    x is a matrix with one channel per column, each embedded with its own
    delay[i] and dim[i] (dim 1: unembedded), or with phasespace=True a list
    of already-built phase spaces, one per channel. A point recurs when it
    recurs in every channel.
    param "rec" sets the radius to reach threshold percent recurrence
    (bisecting iter times; the result is reported with warnings.warn); "rad"
    uses threshold as the radius. zscore=1 standardises the data first. norm
    is "none", "euc", "min" or "max". dmin and vmin are the shortest diagonal
    and vertical lines counted. The plotting options of MATLAB are omitted.

    Returns a dict with RP, the recurrence plot as displayed (rotated as
    MATLAB's imrotate(1 - A, 90), so 0 marks a recurrence), and RESULTS, a
    dict with DIM, EMB, DEL, NORM, ZSCORE, RADIUS, Size, REC, DET, MeanL,
    MaxL, EntrL, LAM, MeanV, MaxV, EntrV (no EntrW).
    """
    delay = np.atleast_1d(delay)
    dim = np.atleast_1d(dim)
    if phasespace:
        if not isinstance(x, (list, tuple)):
            raise ValueError(
                "PhaseSpace is true, so DATA must be a list with one already-reconstructed "
                "phase space per channel, e.g. [Y1, Y2], rather than a raw multi-column series.")
        if len(x) < 2:
            raise ValueError("DATA must contain at least two channels.")
        x = [np.asarray(d, dtype=np.float64) for d in x]
        x = [d[:, None] if d.ndim == 1 else d for d in x]
    else:
        try:
            x = np.asarray(x, dtype=np.float64)
        except ValueError:
            raise ValueError("DATA must be a numeric matrix unless PhaseSpace is true.") from None
        if x.ndim == 1:
            x = x[:, None]
        if x.shape[1] < 2:
            raise ValueError("Data must be at least two column vectors.")
    if zscore:
        x = [zscore_cols(d) for d in x] if phasespace else zscore_cols(x)
    if phasespace:
        data2 = x
        if norm == "euc":
            warnings.warn("PhaseSpace is true, so DATA was not rebuilt from TAU and DIM and DIM "
            "cannot be confirmed to match its width. Norm=\"euc\" scales by DIM "
            "regardless; pass DIM to match each phase space's width or use "
                          "Norm=\"none\", \"min\" or \"max\".", stacklevel=2)
    else:
        data2 = [psr(x[:, i], delay[i], dim[i]) if dim[i] > 1 else x[:, i:i + 1]
                 for i in range(x.shape[1])]
    DIM = len(data2)
    if norm == "euc" and np.unique(dim).size > 1:
        raise ValueError('Norm="euc" scales by DIM, which is per channel here and differs between '
                         'channels. Use Norm="none", "min" or "max".')
    a = [rqa_normalise(-pairwise_dist(d, d), norm, DIM * dim[0]) for d in data2]
    return rqa_finish(x, a, a[0], param, threshold, "jrqa", iter, 0, dmin, vmin,
                      dict(DIM=1, EMB=dim, DEL=delay, NORM=norm, ZSCORE=zscore), False)
