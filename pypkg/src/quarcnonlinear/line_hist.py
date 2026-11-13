"""Line histograms of a recurrence plot."""

# Copyright (c) 2021-2026 Quantitative Analysis Research Core,
# Center for Human Movement Variability, University of Nebraska at Omaha.
# MIT licence. See LICENSE.txt.

import numpy as np

from . import _core


def line_hist(x, a, radius, type, theiler=0):
    """Line histograms and percent recurrence, as MATLAB's line_hist.

    a is a negated distance matrix, or a list of them (joint recurrence: the
    thresholded matrices are multiplied together). type is "rqa", "crqa",
    "jrqa" or "mdrqa"; cross recurrence counts every diagonal and ignores
    theiler. Diagonals within theiler of the main diagonal are excluded, and
    the main diagonal always is, except in cross recurrence. x is unused,
    kept for parity with MATLAB. Returns a dict with recurrence (percent),
    diag_hist, vertical_hist and a, the 0/1 recurrence matrix.
    """
    cross = "crqa" in type
    if cross:
        theiler = 0
    theiler = int(theiler)
    mats = a if isinstance(a, (list, tuple)) else [a]
    mats = [(np.asarray(m, dtype=np.float64) + radius >= 0).astype(np.float64) for m in mats]
    if len(mats) > 1:
        s = min(m.shape[1] for m in mats)
        A = mats[0][:s, :s]
        for m in mats[1:]:
            A = A * m[:s, :s]
    else:
        A = mats[0]
    n = A.shape[0]
    if A.shape[0] == A.shape[1]:
        flat = np.ascontiguousarray((A != 0).ravel(order="F"), dtype=np.uint8)
        diag, vert = _core.line_hist(flat, n, theiler, cross)
        diag, vert = np.asarray(diag), np.asarray(vert)
    else:
        diag, vert = _lines(A, theiler, cross)
    with np.errstate(invalid="ignore", divide="ignore"):   # n = 1: NaN, as MATLAB
        recurrence = _recurrence(A, n, theiler, cross)
    return dict(recurrence=recurrence, diag_hist=diag, vertical_hist=vert, a=A)


def _recurrence(A, n, theiler, cross):
    if not cross:
        band = sum(np.trace(A, k) for k in range(-theiler, theiler + 1))
        excluded = n + 2 * sum(n - k for k in range(1, theiler + 1))
        return 100 * (A.sum() - band) / (n ** 2 - excluded)
    return 100 * A.sum() / n ** 2


def _runs(v):
    """Lengths of the runs of 1s in v."""
    d = np.diff(np.concatenate([[0], (np.asarray(v) != 0).astype(int), [0]]))
    return np.flatnonzero(d == -1) - np.flatnonzero(d == 1)


def _lines(A, theiler, cross):
    """Line counts in plain numpy, as the MATLAB loops; used for non-square
    matrices and as the reference the tests hold the compiled kernel to."""
    n = A.shape[0]
    diag = [_runs(np.diagonal(A, k)) for k in range(-(A.shape[0] - 1), A.shape[1])
            if cross or abs(k) > theiler]
    diag = np.concatenate(diag) if diag else np.zeros(0)
    if not cross and diag.size == 0:
        diag = np.zeros(1)
    vert = []
    for j in range(A.shape[1]):
        v = A[:, j].copy()
        if theiler > 0:
            v[max(0, j - theiler) : min(n, j + theiler + 1)] = 0
        vert.append(_runs(v))
    vert = np.concatenate(vert) if vert else np.zeros(0)
    return diag.astype(np.float64), vert.astype(np.float64)
