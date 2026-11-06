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
#' @param bins Histogram bins per axis: a positive integer count, or a
#'   vector of strictly increasing bin edges spanning the data (each bin
#'   closed on the left, the last also on the right). Default `NULL` selects
#'   the count by Scott's rule. Used only by `"histogram"`.
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
