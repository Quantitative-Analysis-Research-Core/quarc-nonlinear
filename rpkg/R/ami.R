#' Average mutual information versus lag
#'
#' Returns the lag of the first local minimum of the average mutual
#' information of `x` with itself over lags `0:maxlag`: the embedding delay
#' for phase space reconstruction (Fraser and Swinney 1986).
#'
#' @param x Numeric vector, no `NA`.
#' @param maxlag Maximum lag, positive integer, less than `length(x)`.
#' @param algorithm The estimator:
#'   * `"histogram"` (default; alias `"stergiou"`): equal-width joint
#'     histogram. Fast, biased upward.
#'   * `"kde"` (alias `"thomas"`): bivariate Gaussian kernel density.
#'   * `"kde-adaptive"`: as `"kde"`, kernel widths set by local density
#'     (Abramson's square-root law).
#'   * `"adaptive"` (alias `"partition"`): Darbellay-Vajda adaptive partition.
#'   * `"knn"` (alias `"ksg"`): Kraskov-Stoegbauer-Grassberger.
#'   * `"copula"` (alias `"gcmi"`): Gaussian copula. Blind to non-monotonic
#'     dependence; see [ami_copula()] before using it to choose a delay.
#' @param bins Histogram bins per axis. Default `NULL` selects by Scott's
#'   rule. Used only by `"histogram"`.
#' @param k Number of neighbours. Default `NULL` uses 10. Used only by
#'   `"knn"`.
#' @param fraction Fallback threshold, as a fraction of AMI at lag 0, used
#'   when the curve has no local minimum. Default 0.2.
#' @return A list with `tau`, the selected lag; `curve`, a
#'   `(maxlag + 1)`-by-2 matrix of lag and AMI in bits; and `info`, with
#'   `allMinima`, `usedFallback`, `fractionLag`, `estimator` and `algorithm`
#'   plus estimator-specific fields.
#' @references Fraser, A. M. and Swinney, H. L. (1986). Independent
#'   coordinates for strange attractors from mutual information. Physical
#'   Review A, 33(2), 1134-1140.
#' @seealso [ami_histogram()], [ami_kde()], [ami_partition()], [ami_knn()],
#'   [ami_copula()]
#' @export
ami <- function(x, maxlag, algorithm = "histogram", bins = NULL, k = NULL,
                fraction = 0.2) {
  algorithm <- tolower(algorithm)
  if (!is.null(bins) && !algorithm %in% c("histogram", "stergiou")) {
    warning(sprintf('bins applies to algorithm="histogram" and is ignored for "%s".',
                    algorithm), call. = FALSE)
  }
  if (!is.null(k) && !algorithm %in% c("knn", "ksg")) {
    warning(sprintf('k applies to algorithm="knn" and is ignored for "%s".',
                    algorithm), call. = FALSE)
  }
  out <- switch(algorithm,
    histogram = , stergiou = ami_histogram(x, maxlag, bins = bins, fraction = fraction),
    kde = , thomas = ami_kde(x, maxlag, fraction = fraction),
    "kde-adaptive" = ami_kde(x, maxlag, fraction = fraction, bandwidth = "adaptive"),
    adaptive = , partition = ami_partition(x, maxlag, fraction = fraction),
    knn = , ksg = ami_knn(x, maxlag, k = if (is.null(k)) 10 else k, fraction = fraction),
    copula = , gcmi = ami_copula(x, maxlag, fraction = fraction),
    stop(sprintf(paste0('Unknown algorithm "%s". Supported: "histogram" ',
                        '(alias "stergiou"), "kde" (alias "thomas"), "kde-adaptive", ',
                        '"adaptive" (alias "partition"), "knn" (alias "ksg"), ',
                        '"copula" (alias "gcmi").'), algorithm),
         call. = FALSE))
  out$info$algorithm <- algorithm
  out
}

#' @describeIn ami Equal-width joint histogram estimator. Every lag uses the
#'   same number of pairs, `N - maxlag`, so the curve is not confounded by a
#'   changing sample size. Biased upward; adequate for locating a minimum.
#' @export
ami_histogram <- function(x, maxlag, bins = NULL, fraction = 0.2) {
  x <- check_ami_input(x, maxlag, "ami_histogram")
  N <- length(x)
  nBins <- if (is.null(bins)) scott_bins(x) else bins
  m <- N - maxlag
  edges <- bin_edges(x, nBins)
  curve <- cbind(lag = 0:maxlag, ami = 0)
  for (lag in 0:maxlag) {
    curve[lag + 1, 2] <- histogram_mi_bits(x[1:m], x[(1 + lag):(m + lag)], edges)
  }
  out <- ami_first_minimum(curve, fraction, "ami_histogram")
  out$info$bins <- nBins
  out$info$samplesPerLag <- m
  out$info$estimator <- "histogram"
  out
}

#' @describeIn ami Gaussian kernel density estimator (Thomas et al. 2014).
#'   Lag `k` uses all `N - k` pairs. Cost is O(maxlag * N^2).
#' @param chunk Rows of the pairwise kernel or distance matrix evaluated at
#'   once. Changes peak memory, not the result. Default 512.
#' @param bandwidth `"fixed"` (default) or `"adaptive"`: each kernel scaled
#'   by `sqrt(g / pilot_i)`, `g` the geometric mean of a fixed-bandwidth
#'   pilot (Abramson 1982). About twice the cost.
#' @export
ami_kde <- function(x, maxlag, fraction = 0.2, chunk = 512,
                    bandwidth = c("fixed", "adaptive")) {
  bandwidth <- match.arg(bandwidth)
  adaptive <- bandwidth == "adaptive"
  x <- check_ami_input(x, maxlag, "ami_kde")
  N <- length(x)
  curve <- cbind(lag = 0:maxlag, ami = 0)
  curve[1, 2] <- kde_mi_bits(x, x, chunk, adaptive)
  for (lag in seq_len(maxlag)) {
    curve[lag + 1, 2] <- kde_mi_bits(x[1:(N - lag)], x[(lag + 1):N], chunk, adaptive)
  }
  out <- ami_first_minimum(curve, fraction, "ami_kde")
  out$info$bandwidth <- bandwidth
  out$info$estimator <- if (adaptive) "kde-adaptive" else "kde"
  out
}

#' @describeIn ami Adaptive partition of the rank plane (Fraser and Swinney
#'   1986; Darbellay and Vajda 1999). A cell is cut into four at its median
#'   ranks when a chi-square test (3 df) rejects a uniform spread; the root
#'   is always cut. Every lag uses `N - maxlag` pairs.
#' @param alpha Significance level of the splitting test. Default 0.05.
#' @export
ami_partition <- function(x, maxlag, alpha = 0.05, fraction = 0.2) {
  x <- check_ami_input(x, maxlag, "ami_partition")
  if (!(alpha > 0 && alpha < 1)) stop("alpha must lie in (0, 1).", call. = FALSE)
  N <- length(x)
  m <- N - maxlag
  critical <- stats::qchisq(1 - alpha, 3)
  curve <- cbind(lag = 0:maxlag, ami = 0)
  for (lag in 0:maxlag) {
    curve[lag + 1, 2] <- partition_bits(x[1:m], x[(1 + lag):(m + lag)], critical)
  }
  out <- ami_first_minimum(curve, fraction, "ami_partition")
  out$info$samplesPerLag <- m
  out$info$alpha <- alpha
  out$info$estimator <- "adaptive"
  out
}

#' @describeIn ami Kraskov-Stoegbauer-Grassberger estimator, algorithm 1.
#'   The series is standardised and moved by `1e-10 * (((1:N) * g) %% 1 -
#'   0.5)`, `g` the golden-ratio conjugate, to break ties deterministically,
#'   exactly as in MATLAB. Every lag uses `N - maxlag` pairs. Cost
#'   O(maxlag * N^2). Can return small negative values where the true AMI
#'   is zero.
#' @export
ami_knn <- function(x, maxlag, k = 10, fraction = 0.2, chunk = 512) {
  x <- check_ami_input(x, maxlag, "ami_knn")
  N <- length(x)
  m <- N - maxlag
  if (k >= m) {
    stop(sprintf("ami_knn: k (%d) must be smaller than the number of pairs per lag (%d).",
                 as.integer(k), m), call. = FALSE)
  }
  s <- stats::sd(x)
  if (s > 0) x <- (x - mean(x)) / s
  x <- x + 1e-10 * (((1:N) * 0.6180339887498949) %% 1 - 0.5)
  curve <- cbind(lag = 0:maxlag, ami = 0)
  for (lag in 0:maxlag) {
    curve[lag + 1, 2] <- ksg_bits(x[1:m], x[(1 + lag):(m + lag)], k, chunk)
  }
  out <- ami_first_minimum(curve, fraction, "ami_knn")
  out$info$samplesPerLag <- m
  out$info$k <- k
  out$info$estimator <- "knn"
  out
}

#' @describeIn ami Gaussian copula estimator (Ince et al. 2017): ranks
#'   mapped to normal scores, then `-0.5 * log2(1 - r^2)`. Exact for a
#'   Gaussian copula and a lower bound otherwise. It sees only dependence a
#'   correlation of normal scores can see: on a sine at a quarter-period lag
#'   it is near zero although the true AMI is large. `r` is capped at
#'   `1 - 1e-6`. Every lag uses `N - maxlag` pairs.
#' @export
ami_copula <- function(x, maxlag, fraction = 0.2) {
  x <- check_ami_input(x, maxlag, "ami_copula")
  N <- length(x)
  m <- N - maxlag
  curve <- cbind(lag = 0:maxlag, ami = 0)
  for (lag in 0:maxlag) {
    a <- normal_scores(x[1:m])
    b <- normal_scores(x[(1 + lag):(m + lag)])
    r <- sum(a * b) / (sqrt(sum(a^2)) * sqrt(sum(b^2)))
    r <- min(max(r, -0.999999), 0.999999)
    curve[lag + 1, 2] <- -0.5 * log2(1 - r^2)
  }
  out <- ami_first_minimum(curve, fraction, "ami_copula")
  out$info$samplesPerLag <- m
  out$info$estimator <- "copula"
  out
}

# ------------------------------------------------------------------ helpers

check_ami_input <- function(x, maxlag, caller) {
  x <- as.numeric(x)
  if (anyNA(x)) {
    stop(sprintf(paste0("%s: x contains NA. Mutual information is undefined for ",
                        "missing data; remove or impute first."), caller), call. = FALSE)
  }
  if (maxlag >= length(x)) {
    stop(sprintf("%s: maxlag (%d) must be smaller than the number of samples (%d).",
                 caller, as.integer(maxlag), length(x)), call. = FALSE)
  }
  x
}

# Scott's rule: bin width h = 3.49 * sd * N^(-1/3).
scott_bins <- function(x) {
  s <- stats::sd(x)
  if (s == 0) return(1)
  h <- 3.49 * s * length(x)^(-1 / 3)
  max(2, ceiling((max(x) - min(x)) / h))
}

# Equal-width edges spanning the data, computed as MATLAB's linspace does,
# so a sample that falls exactly on an edge lands in the same bin.
bin_edges <- function(x, nBins) {
  lo <- min(x)
  hi <- max(x)
  if (hi == lo) hi <- lo + 1
  e <- lo + (0:nBins) * (hi - lo) / nBins
  e[nBins + 1] <- hi
  e
}

histogram_mi_bits <- function(a, b, edges) {
  n <- length(a)
  nb <- length(edges) - 1
  # Left-closed bins with the last one closed on the right, as histcounts.
  ia <- findInterval(a, edges, rightmost.closed = TRUE)
  ib <- findInterval(b, edges, rightmost.closed = TRUE)
  pAB <- matrix(tabulate((ib - 1) * nb + ia, nb * nb), nb, nb) / n
  nz <- pAB > 0
  if (!any(nz)) return(0)
  outer_ <- outer(rowSums(pAB), colSums(pAB))
  sum(pAB[nz] * log2(pAB[nz] / outer_[nz]))
}

kde_mi_bits <- function(X, Y, chunk, adaptive = FALSE) {
  n <- length(X)
  hx <- stats::sd(X) / n^(1 / 6)
  hy <- stats::sd(Y) / n^(1 / 6)
  if (hx == 0 || hy == 0) return(0)
  pX <- kde_marginal(X, hx, chunk)
  pY <- kde_marginal(Y, hy, chunk)
  pXY <- kde_joint(X, Y, hx, hy, chunk)
  if (adaptive) {
    # The fixed-bandwidth estimates above are the pilots.
    pX <- kde_marginal(X, hx, chunk, abramson(pX))
    pY <- kde_marginal(Y, hy, chunk, abramson(pY))
    pXY <- kde_joint(X, Y, hx, hy, chunk, abramson(pXY))
  }
  ok <- pXY > 0 & pX > 0 & pY > 0
  sum(log2(pXY[ok] / (pX[ok] * pY[ok]))) / n
}

kde_marginal <- function(v, h, chunk, lambda = NULL) {
  n <- length(v)
  p <- numeric(n)
  cst <- 1 / (n * h * sqrt(2 * pi))
  for (lo in seq(1, n, by = chunk)) {
    rows <- lo:min(lo + chunk - 1, n)
    z <- outer(v[rows], v, `-`) / h
    if (is.null(lambda)) {
      p[rows] <- cst * rowSums(exp(-0.5 * z^2))
    } else {
      L <- matrix(lambda, length(rows), n, byrow = TRUE)
      p[rows] <- cst * rowSums(exp(-0.5 * (z / L)^2) / L)
    }
  }
  p
}

kde_joint <- function(X, Y, hx, hy, chunk, lambda = NULL) {
  n <- length(X)
  r <- sum((X - mean(X)) * (Y - mean(Y))) /
    (sqrt(sum((X - mean(X))^2)) * sqrt(sum((Y - mean(Y))^2)))
  r <- min(max(r, -0.999999), 0.999999)
  omr2 <- 1 - r^2
  cst <- 1 / (n * 2 * pi * hx * hy * sqrt(omr2))
  p <- numeric(n)
  for (lo in seq(1, n, by = chunk)) {
    rows <- lo:min(lo + chunk - 1, n)
    u <- outer(X[rows], X, `-`) / hx
    v <- outer(Y[rows], Y, `-`) / hy
    q <- (u^2 - 2 * r * (u * v) + v^2) / omr2
    if (is.null(lambda)) {
      p[rows] <- cst * rowSums(exp(-0.5 * q))
    } else {
      L2 <- matrix(lambda^2, length(rows), n, byrow = TRUE)
      p[rows] <- cst * rowSums(exp(-0.5 * q / L2) / L2)
    }
  }
  p
}

# Abramson's square-root law: lambda_i = sqrt(g / pilot_i).
abramson <- function(pilot) {
  pilot <- pmax(pilot, .Machine$double.xmin)
  g <- exp(mean(log(pilot)))
  sqrt(g / pilot)
}

# Darbellay-Vajda estimate of I(A;B) in bits, on ranks.
partition_bits <- function(a, b, critical) {
  n <- length(a)
  rx <- rank_of(a)
  ry <- rank_of(b)
  # Explicit stack of cells: bounds [xlo xhi ylo yhi] in ranks, members.
  stackB <- list(c(1, n, 1, n))
  stackM <- list(seq_len(n))
  isRoot <- TRUE
  I <- 0
  while (length(stackB)) {
    top <- length(stackB)
    cb <- stackB[[top]]
    idx <- stackM[[top]]
    stackB[[top]] <- NULL
    stackM[[top]] <- NULL
    nc <- length(idx)
    if (nc == 0) next
    wx <- cb[2] - cb[1] + 1
    wy <- cb[4] - cb[3] + 1
    split <- FALSE
    if (wx >= 2 && wy >= 2) {
      xm <- cb[1] + floor(wx / 2) - 1
      ym <- cb[3] + floor(wy / 2) - 1
      lx <- rx[idx] <= xm
      ly <- ry[idx] <= ym
      q <- c(sum(lx & ly), sum(lx & !ly), sum(!lx & ly), sum(!lx & !ly))
      fx <- (xm - cb[1] + 1) / wx
      fy <- (ym - cb[3] + 1) / wy
      E <- nc * c(fx * fy, fx * (1 - fy), (1 - fx) * fy, (1 - fx) * (1 - fy))
      split <- isRoot || sum((q - E)^2 / E) > critical
    }
    isRoot <- FALSE
    if (split) {
      # Pushed in MATLAB's order so the summation order matches.
      stackB <- c(stackB, list(c(cb[1], xm, cb[3], ym), c(cb[1], xm, ym + 1, cb[4]),
                               c(xm + 1, cb[2], cb[3], ym), c(xm + 1, cb[2], ym + 1, cb[4])))
      stackM <- c(stackM, list(idx[lx & ly], idx[lx & !ly], idx[!lx & ly], idx[!lx & !ly]))
    } else {
      I <- I + (nc / n) * log2(nc * n / (wx * wy))
    }
  }
  I
}

# Ranks 1..n, ties broken by position (order() is stable, like MATLAB sort).
rank_of <- function(v) {
  r <- integer(length(v))
  r[order(v)] <- seq_along(v)
  r
}

# Ranks mapped to standard normal quantiles.
normal_scores <- function(v) {
  n <- length(v)
  u <- numeric(n)
  u[order(v)] <- seq_len(n) / (n + 1)
  stats::qnorm(u)
}

# Kraskov-Stoegbauer-Grassberger algorithm 1, in bits.
ksg_bits <- function(a, b, k, chunk) {
  n <- length(a)
  nx <- numeric(n)
  ny <- numeric(n)
  for (lo in seq(1, n, by = chunk)) {
    rows <- lo:min(lo + chunk - 1, n)
    dx <- abs(outer(a[rows], a, `-`))
    dy <- abs(outer(b[rows], b, `-`))
    d <- pmax(dx, dy)
    d[cbind(seq_along(rows), rows)] <- Inf          # exclude self
    e <- apply(d, 1, function(r) sort.int(r, partial = k)[k])
    # Strictly closer than eps_i; the self term is counted and removed.
    nx[rows] <- rowSums(dx < e) - 1
    ny[rows] <- rowSums(dy < e) - 1
  }
  (digamma(k) + digamma(n) - mean(digamma(nx + 1) + digamma(ny + 1))) / log(2)
}

# Lag of the first strict local minimum, falling back to the first lag at
# which AMI drops to `fraction` of AMI(0) (Abarbanel et al. 1993).
ami_first_minimum <- function(curve, fraction, caller) {
  v <- curve[, 2]
  k <- length(v)
  isMin <- c(FALSE, v[2:(k - 1)] < v[1:(k - 2)] & v[2:(k - 1)] <= v[3:k], FALSE)
  if (k < 3) isMin <- rep(FALSE, k)
  idx <- which(isMin)
  info <- list(allMinima = curve[idx, , drop = FALSE], usedFallback = FALSE,
               fractionLag = NA_real_)
  below <- which(v <= fraction * v[1])
  if (length(below)) info$fractionLag <- curve[below[1], 1]
  if (length(idx)) {
    tau <- curve[idx[1], 1]
  } else {
    tau <- info$fractionLag
    info$usedFallback <- TRUE
    if (is.na(tau)) {
      warning(sprintf(paste0("%s: no local minimum in AMI over lags 0:%d and the curve ",
                             "never fell to %.0f%% of its lag-0 value."),
                      caller, as.integer(curve[k, 1]), 100 * fraction), call. = FALSE)
    }
  }
  list(tau = unname(tau), curve = curve, info = info)
}
