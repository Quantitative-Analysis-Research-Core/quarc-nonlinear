"""Average mutual information versus lag, dispatching to an estimator."""

# Copyright (c) 2021-2026 Quantitative Analysis Research Core,
# Center for Human Movement Variability, University of Nebraska at Omaha.
# MIT licence. See LICENSE.txt.

import warnings

from .ami_copula import ami_copula
from .ami_histogram import ami_histogram
from .ami_kde import ami_kde
from .ami_knn import ami_knn
from .ami_partition import ami_partition


def ami(x, maxlag, algorithm="histogram", bins=None, k=None, fraction=0.2):
    """Embedding delay: the first local minimum of AMI over lags 0..maxlag.

    algorithm is "histogram" (alias "stergiou"), "kde" (alias "thomas"),
    "kde-adaptive", "adaptive" (alias "partition"), "knn" (alias "ksg") or
    "copula" (alias "gcmi"). bins is used only by "histogram" (None selects
    Scott's rule), k only by "knn" (None uses 10). fraction is the fallback
    threshold, as a fraction of AMI(0), when the curve has no local minimum.
    Returns (tau, curve, info): the selected lag (nan if none), a
    (maxlag+1, 2) array of lag and AMI in bits, and a dict with allMinima,
    usedFallback, fractionLag, estimator and algorithm plus
    estimator-specific fields.
    """
    algorithm = algorithm.lower()
    if bins is not None and algorithm not in ("histogram", "stergiou"):
        warnings.warn(f'bins applies to algorithm="histogram" and is ignored for '
                      f'"{algorithm}".', stacklevel=2)
    if k is not None and algorithm not in ("knn", "ksg"):
        warnings.warn(f'k applies to algorithm="knn" and is ignored for "{algorithm}".',
                      stacklevel=2)
    if algorithm in ("histogram", "stergiou"):
        out = ami_histogram(x, maxlag, bins=bins, fraction=fraction)
    elif algorithm in ("kde", "thomas"):
        out = ami_kde(x, maxlag, fraction=fraction)
    elif algorithm == "kde-adaptive":
        out = ami_kde(x, maxlag, fraction=fraction, bandwidth="adaptive")
    elif algorithm in ("adaptive", "partition"):
        out = ami_partition(x, maxlag, fraction=fraction)
    elif algorithm in ("knn", "ksg"):
        out = ami_knn(x, maxlag, k=10 if k is None else k, fraction=fraction)
    elif algorithm in ("copula", "gcmi"):
        out = ami_copula(x, maxlag, fraction=fraction)
    else:
        raise ValueError(f'Unknown algorithm "{algorithm}". Supported: "histogram" '
                         '(alias "stergiou"), "kde" (alias "thomas"), "kde-adaptive", '
                         '"adaptive" (alias "partition"), "knn" (alias "ksg"), '
                         '"copula" (alias "gcmi").')
    out[2]["algorithm"] = algorithm
    return out
