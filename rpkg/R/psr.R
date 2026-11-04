#' Phase space reconstruction by time-delay embedding
#'
#' Builds the delay-coordinate matrix of a time series. A matrix input is
#' embedded column by column, so each delay block holds every input channel.
#'
#' @param x Numeric vector, or a matrix with one channel per column.
#' @param delay Embedding delay in samples, positive integer.
#' @param dim Embedding dimension, positive integer.
#' @return A matrix with `N - (dim - 1) * delay` rows and `dim * ncol(x)`
#'   columns. Columns `1:ncol(x)` are the undelayed channels.
#' @examples
#' x <- sin(seq(0, 20 * pi, length.out = 2000))
#' Y <- psr(x, delay = 25, dim = 3)
#' dim(Y)
#' @export
psr <- function(x, delay, dim) {
  x <- as.matrix(x)
  N <- nrow(x)
  D <- ncol(x)
  M <- N - (dim - 1) * delay
  if (M < 1) {
    stop(sprintf("Series of %d samples is too short for dim = %d at delay = %d.",
                 N, dim, delay), call. = FALSE)
  }
  out <- matrix(0, M, dim * D)
  for (i in seq_len(dim)) {
    rows <- (1 + (i - 1) * delay):(N - (dim - i) * delay)
    out[, (1:D) + D * (i - 1)] <- x[rows, , drop = FALSE]
  }
  out
}
