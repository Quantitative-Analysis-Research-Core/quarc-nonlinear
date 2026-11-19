"""Permutation entropy."""

# Copyright (c) 2021-2026 Quantitative Analysis Research Core,
# Center for Human Movement Variability, University of Nebraska at Omaha.
# MIT licence. See LICENSE.txt.

from itertools import permutations

import numpy as np


def ent_permu(x, dim, delay):
    """Permutation entropy (Bandt & Pompe 2002), in bits, as MATLAB's ent_permu.

    Returns (permEnt, hist): the entropy and the count of each ordinal
    pattern. hist is in the order of MATLAB's perms(1:dim) (reverse
    lexicographic), so it compares entry by entry with MATLAB's output. Ties
    within a window resolve by position, as MATLAB's stable sort.
    """
    x = np.asarray(x, dtype=np.float64).ravel()
    N = x.size
    nWin = N - delay * (dim - 1)
    # MATLAB's perms(1:dim): permutations of dim..1 in lexicographic order.
    index = {p: i for i, p in enumerate(permutations(range(dim - 1, -1, -1)))}
    hist = np.zeros(len(index))
    for s in range(nWin):
        w = x[s:s + delay * (dim - 1) + 1:delay]
        hist[index[tuple(np.argsort(w, kind="stable"))]] += 1
    p = hist[hist != 0] / np.sum(hist)
    return float(-np.sum(p * np.log2(p))), hist
