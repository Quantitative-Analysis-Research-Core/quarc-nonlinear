"""Euclidean distances between the rows of two matrices."""

# Copyright (c) 2021-2026 Quantitative Analysis Research Core,
# Center for Human Movement Variability, University of Nebraska at Omaha.
# MIT licence. See LICENSE.txt.

import numpy as np


def pairwise_dist(X, Y):
    """Matrix of Euclidean distances, row i of X to row j of Y.

    Squares are summed column by column, in the order R's pairwise_dist does.
    """
    X = np.asarray(X, dtype=np.float64)
    Y = np.asarray(Y, dtype=np.float64)
    if X.ndim == 1:
        X = X[:, None]
    if Y.ndim == 1:
        Y = Y[:, None]
    s = np.zeros((X.shape[0], Y.shape[0]))
    for k in range(X.shape[1]):
        s += np.subtract.outer(X[:, k], Y[:, k]) ** 2
    return np.sqrt(s)
