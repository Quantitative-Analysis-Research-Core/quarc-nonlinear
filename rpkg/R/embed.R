#' Delay embedding with arbitrary lags
#'
#' @param z Numeric vector.
#' @param v Lags, or, with `w`, the number of lags.
#' @param w Spacing between lags: when given, the lags are `0:(v-1) * w`.
#' @return A matrix with one row per lag, in increasing lag order, and one
#'   column per embedded time. Row `i` is `z` delayed by the `i`-th smallest
#'   lag, so the lag-0 row holds the latest sample of each vector. This is the
#'   transpose of [psr()]'s layout, as in MATLAB.
#' @export
embed <- function(z, v = c(0, 1, 2), w = NULL) {
  if (!is.null(w)) v <- seq(0, w * (v - 1), by = w)
  z <- as.numeric(z)
  lags <- sort(v)
  dim <- length(lags)
  n <- length(z)
  if (n < lags[dim]) stop("Vector is too small to be embedded with the given lags", call. = FALSE)
  m <- n - (lags[dim] - lags[1])
  t <- (1:m) + lags[dim]
  x <- matrix(0, dim, m)
  for (i in seq_len(dim)) x[i, ] <- z[t - lags[i]]
  x
}
