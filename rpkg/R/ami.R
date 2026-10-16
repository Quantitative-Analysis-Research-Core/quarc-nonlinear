#' Average mutual information versus lag
#'
#' Returns the lag of the first local minimum of the average mutual
#' information of `x` with itself over lags `0:maxlag`: the embedding delay
#' for phase space reconstruction (Fraser and Swinney 1986).
#'
#' @param x Numeric vector, no `NA`.
#' @param maxlag Maximum lag, positive integer, less than `length(x)`.
#' @param algorithm `"histogram"` (default; alias `"stergiou"`), an
#'   equal-width joint histogram, or `"kde"` (alias `"thomas"`), a bivariate
#'   Gaussian kernel density estimate: slower, less biased.
#' @param bins Histogram bins per axis. Default `NULL` selects by Scott's
#'   rule. Used only by `"histogram"`.
#' @param fraction Fallback threshold, as a fraction of AMI at lag 0, used
#'   when the curve has no local minimum. Default 0.2.
#' @return A list with `tau`, the selected lag; `curve`, a
#'   `(maxlag + 1)`-by-2 matrix of lag and AMI in bits; and `info`, with
#'   `allMinima`, `usedFallback`, `fractionLag`, `estimator` and `algorithm`
#'   (plus `bins` and `samplesPerLag` for `"histogram"`).
#' @references Fraser, A. M. and Swinney, H. L. (1986). Independent
#'   coordinates for strange attractors from mutual information. Physical
#'   Review A, 33(2), 1134-1140.
#' @seealso [ami_histogram()], [ami_kde()]
#' @export
ami <- function(x, maxlag, algorithm = "histogram", bins = NULL, fraction = 0.2) {
  algorithm <- tolower(algorithm)
  if (!is.null(bins) && !algorithm %in% c("histogram", "stergiou")) {
    warning(sprintf('bins applies to algorithm="histogram" and is ignored for "%s".',
                    algorithm), call. = FALSE)
  }
  out <- switch(algorithm,
    histogram = , stergiou = ami_histogram(x, maxlag, bins = bins, fraction = fraction),
    kde = , thomas = ami_kde(x, maxlag, fraction = fraction),
    stop(sprintf(paste0('Unknown algorithm "%s". Supported: "histogram" ',
                        '(alias "stergiou"), "kde" (alias "thomas").'), algorithm),
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
#' @param chunk Rows of the pairwise kernel evaluated at once. Changes peak
#'   memory, not the result. Default 512.
#' @export
ami_kde <- function(x, maxlag, fraction = 0.2, chunk = 512) {
  x <- check_ami_input(x, maxlag, "ami_kde")
  N <- length(x)
  curve <- cbind(lag = 0:maxlag, ami = 0)
  curve[1, 2] <- kde_mi_bits(x, x, chunk)
  for (lag in seq_len(maxlag)) {
    curve[lag + 1, 2] <- kde_mi_bits(x[1:(N - lag)], x[(lag + 1):N], chunk)
  }
  out <- ami_first_minimum(curve, fraction, "ami_kde")
  out$info$estimator <- "kde"
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

kde_mi_bits <- function(X, Y, chunk) {
  n <- length(X)
  hx <- stats::sd(X) / n^(1 / 6)
  hy <- stats::sd(Y) / n^(1 / 6)
  if (hx == 0 || hy == 0) return(0)
  pX <- kde_marginal(X, hx, chunk)
  pY <- kde_marginal(Y, hy, chunk)
  pXY <- kde_joint(X, Y, hx, hy, chunk)
  ok <- pXY > 0 & pX > 0 & pY > 0
  sum(log2(pXY[ok] / (pX[ok] * pY[ok]))) / n
}

kde_marginal <- function(v, h, chunk) {
  n <- length(v)
  p <- numeric(n)
  cst <- 1 / (n * h * sqrt(2 * pi))
  for (lo in seq(1, n, by = chunk)) {
    rows <- lo:min(lo + chunk - 1, n)
    z <- outer(v[rows], v, `-`) / h
    p[rows] <- cst * rowSums(exp(-0.5 * z^2))
  }
  p
}

kde_joint <- function(X, Y, hx, hy, chunk) {
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
    p[rows] <- cst * rowSums(exp(-0.5 * q))
  }
  p
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
