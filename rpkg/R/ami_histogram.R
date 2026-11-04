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
