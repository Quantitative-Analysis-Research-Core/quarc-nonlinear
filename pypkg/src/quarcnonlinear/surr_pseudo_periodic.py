"""Pseudo-periodic surrogate."""

# Copyright (c) 2021-2026 Quantitative Analysis Research Core,
# Center for Human Movement Variability, University of Nebraska at Omaha.
# MIT licence. See LICENSE.txt.

import numpy as np

from .psr import psr


def surr_pseudo_periodic(x, delay, dim, rho, rng=None):
    """Pseudo-periodic surrogate (Small, Yu & Harrison 2001), as MATLAB's
    surr_pseudo_periodic.

    A new series is built by hopping through the reconstructed attractor:
    from the current point the next is chosen with probability falling off as
    exp(-distance / rho), and its successor is emitted. rng is a numpy
    Generator or seed. Returns (ys, yi): the surrogate and the 0-based index
    of x each sample was taken from.
    """
    rng = np.random.default_rng(rng)
    x = np.asarray(x, dtype=np.float64).ravel()
    Y = psr(x, delay, dim)
    nY = Y.shape[0]
    xi = int(np.floor(rng.random() * nY))           # 0-based
    ys = np.zeros(nY)
    yi = np.zeros(nY, dtype=np.int64)
    ys[0] = x[xi]
    yi[0] = xi
    M = nY - 2
    for i in range(1, nY):
        prob = np.exp(-np.sqrt(np.sum((Y[:M] - Y[xi]) ** 2, axis=1)) / rho)
        sum3 = np.cumsum(prob)
        xi_n = None
        ind = 1
        while xi_n is None or xi_n > nY - 1 or xi_n == xi:
            hit = np.flatnonzero(sum3 < rng.random() * sum3[M - 1])
            xi_n = int(hit.max()) + 2 if hit.size else None
            ind += 1
            if ind == 100:
                raise ValueError("a new value of xi could not be found, check that rho is not "
                                 "too low")
        xi = xi_n
        ys[i] = x[xi]
        yi[i] = xi
    return ys, yi
