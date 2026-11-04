#' Fractional Gaussian noise
#'
#' Davies-Harte circulant embedding: exact fGn of length `n` with Hurst
#' exponent `H`, zero mean and unit variance.
#'
#' @param n Length.
#' @param H Hurst exponent, in (0, 1).
#' @return Numeric vector of length `n`.
#' @references Davies, R. B. and Harte, D. S. (1987). Tests for Hurst
#'   effect. Biometrika, 74(1), 95-101.
#' @export
fgn_sim <- function(n, H) {
  z <- stats::rnorm(2 * n)
  zr <- z[1:n]
  zi <- z[(n + 1):(2 * n)]
  zic <- -zi
  zi[1] <- 0
  zr[1] <- zr[1] * sqrt(2)
  zi[n] <- 0
  zr[n] <- zr[n] * sqrt(2)
  back <- if (n > 2) (n - 1):2 else integer(0)
  zr <- c(zr[1:n], zr[back])
  zi <- c(zi[1:n], zic[back])
  z <- complex(real = zr, imaginary = zi)
  k <- 0:(n - 1)
  gammak <- (abs(k - 1)^(2 * H) - 2 * abs(k)^(2 * H) + abs(k + 1)^(2 * H)) / 2
  ind <- c(0:(n - 2), n - 1, if (n > 2) (n - 2):1 else integer(0))
  # MATLAB's ifft(v)*length(z) is an unnormalised inverse transform.
  gk <- Re(stats::fft(gammak[ind + 1], inverse = TRUE)) * length(z) / length(ind)
  if (!all(gk > 0)) stop("Re(gk)-vector not positive", call. = FALSE)
  z <- z * sqrt(gk)
  z <- stats::fft(z, inverse = TRUE) / length(z) * length(z)
  z <- 0.5 * (n - 1)^(-0.5) * z
  Re(z[1:n])
}
