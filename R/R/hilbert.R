#' Analytic signal by the FFT method
#'
#' The same construction as MATLAB's `hilbert`: zero the negative
#' frequencies, double the positive ones, invert.
#'
#' @param x Numeric vector.
#' @return Complex vector; `Re` is `x`, `Im` its Hilbert transform.
#' @export
hilbert <- function(x) {
  n <- length(x)
  h <- numeric(n)
  if (n %% 2 == 0) {
    h[c(1, n / 2 + 1)] <- 1
    h[2:(n / 2)] <- 2
  } else {
    h[1] <- 1
    h[2:((n + 1) / 2)] <- 2
  }
  stats::fft(stats::fft(x) * h, inverse = TRUE) / n
}
