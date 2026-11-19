"""Average mutual information by an equal-width joint histogram."""

# Copyright (c) 2021-2026 Quantitative Analysis Research Core,
# Center for Human Movement Variability, University of Nebraska at Omaha.
# MIT licence. See LICENSE.txt.

import numpy as np

from ._ami_first_minimum import ami_first_minimum
from ._check_ami_input import check_ami_input


def ami_histogram(x, maxlag, bins=None, fraction=0.2):
    """AMI versus lag 0..maxlag from an equal-width joint histogram.

    Every lag uses the same N - maxlag pairs. bins is None (Scott's rule), a
    positive integer count, or strictly increasing edges spanning the data
    (bins closed on the left, the last also on the right). Biased upward;
    adequate for locating a minimum. Returns (tau, curve, info): the lag of
    the first local minimum, a (maxlag+1, 2) array of lag and AMI in bits,
    and a dict with allMinima, usedFallback, fractionLag, bins, edges,
    samplesPerLag and estimator.
    """
    x = check_ami_input(x, maxlag, "ami_histogram")
    N = x.size
    _check_ami_bins(bins)
    if bins is None or np.size(bins) == 1:
        nBins = _scott_bins(x) if bins is None else int(np.ravel(bins)[0])
        edges = _bin_edges(x, nBins)
    else:
        edges = np.asarray(bins, dtype=np.float64).ravel()
        if x.min() < edges[0] or x.max() > edges[-1]:
            raise ValueError(f"ami_histogram: bins edges [{edges[0]:g}, {edges[-1]:g}] do not "
                             f"span the data [{x.min():g}, {x.max():g}].")
        nBins = edges.size - 1
    m = N - maxlag
    curve = np.column_stack([np.arange(maxlag + 1, dtype=np.float64), np.zeros(maxlag + 1)])
    for lag in range(maxlag + 1):
        curve[lag, 1] = _histogram_mi_bits(x[:m], x[lag:m + lag], edges)
    tau, curve, info = ami_first_minimum(curve, fraction, "ami_histogram")
    info.update(bins=nBins, edges=edges, samplesPerLag=m, estimator="histogram")
    return tau, curve, info


def _check_ami_bins(bins):
    if bins is None:
        return
    b = np.asarray(bins)
    if b.size == 1:
        ok = np.issubdtype(b.dtype, np.number) and b.item() >= 1 and b.item() == np.floor(b.item())
    else:
        ok = (np.issubdtype(b.dtype, np.number) and np.isfinite(b).all()
              and (np.diff(b.ravel()) > 0).all())
    if not ok:
        raise ValueError("bins must be a positive integer count or strictly increasing "
                         "finite edges.")


def _scott_bins(x):
    """Scott's rule: bin width h = 3.49 * sd * N^(-1/3)."""
    s = np.std(x, ddof=1)
    if s == 0:
        return 1
    h = 3.49 * s * x.size ** (-1 / 3)
    return int(max(2, np.ceil((x.max() - x.min()) / h)))


def _bin_edges(x, nBins):
    """Equal-width edges spanning the data, as MATLAB's linspace computes them."""
    lo, hi = x.min(), x.max()
    if hi == lo:
        hi = lo + 1
    e = lo + np.arange(nBins + 1) * (hi - lo) / nBins
    e[nBins] = hi
    return e


def _histogram_mi_bits(a, b, edges):
    n = a.size
    nb = edges.size - 1
    # Left-closed bins with the last one closed on the right, as histcounts.
    ia = np.minimum(np.searchsorted(edges, a, side="right"), nb) - 1
    ib = np.minimum(np.searchsorted(edges, b, side="right"), nb) - 1
    # P[j, i] is the joint probability of a in bin i and b in bin j, so C
    # order walks it as R and MATLAB walk their column-major matrix.
    P = np.bincount(ib * nb + ia, minlength=nb * nb).reshape(nb, nb) / n
    nz = P > 0
    if not nz.any():
        return 0.0
    outer_ = np.outer(P.sum(axis=1), P.sum(axis=0))
    return float(np.sum(P[nz] * np.log2(P[nz] / outer_[nz])))
