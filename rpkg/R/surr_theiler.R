#' Surrogate data, Theiler's algorithms 0, 1 and 2
#'
#' Theiler et al. (1992). Algorithm 0 shuffles the series (destroys all
#' temporal structure, keeps the distribution). Algorithm 1 randomises the
#' Fourier phases (keeps the amplitude spectrum exactly, so the
#' autocorrelation; not the distribution). Algorithm 2 is the amplitude
#' adjusted Fourier transform (keeps the distribution exactly and the
#' spectrum approximately).
#'
#' Draws come from R's generator, so a surrogate differs from MATLAB's even
#' at the same seed; what each algorithm preserves is the same, and that is
#' what the tests check.
#'
#' @param x Numeric vector, or a matrix with one series per column.
#' @param algorithm 0, 1 or 2.
#' @return A surrogate of the same shape as `x`.
#' @references Theiler, J., Eubank, S., Longtin, A., Galdrikian, B. and
#'   Farmer, J. D. (1992). Testing for nonlinearity in time series: the
#'   method of surrogate data. Physica D, 58(1-4), 77-94.
#' @export
surr_theiler <- function(x, algorithm) {
  isVec <- is.null(dim(x))
  x <- as.matrix(x)
  z <- switch(as.character(algorithm),
    "0" = apply(x, 2, function(v) v[order(stats::rnorm(length(v)))]),
    "1" = surr1(x, 1),
    "2" = surr1(x, 2),
    stop(sprintf("Unknown algorithm %s. Supported: 0, 1, 2.", format(algorithm)), call. = FALSE))
  z <- matrix(z, nrow(x))
  if (isVec) as.numeric(z) else z
}

surr1 <- function(x, algorithm) {
  r <- nrow(x)
  c <- ncol(x)
  if (algorithm == 2) {
    # Gaussianise: replace each value by the same-rank value of a sorted
    # normal sample.
    sr <- apply(matrix(stats::rnorm(r * c), r, c), 2, sort)
    sx <- apply(x, 2, sort)
    y <- sapply(seq_len(c), function(k) sr[rank_of(x[, k]), k])
    y <- matrix(y, r, c)
  } else {
    sx <- NULL
    y <- x
  }
  m <- colMeans(y)
  y <- sweep(y, 2, m)
  fy <- stats::mvfft(y)
  nHalf <- floor((r - 1) / 2)
  rot <- rep(1 + 0i, r)
  if (nHalf > 0) {
    ph <- 2 * pi * stats::runif(nHalf)
    rot[2:(nHalf + 1)] <- exp(1i * ph)
    rot[r:(r - nHalf + 1)] <- Conj(exp(1i * ph))
  }
  yy <- Re(stats::mvfft(fy * rot, inverse = TRUE)) / r
  yy <- sweep(yy, 2, m, "+")
  if (algorithm == 2) {
    yy <- sapply(seq_len(c), function(k) sx[rank_of(yy[, k]), k])
  }
  matrix(yy, r, c)
}
