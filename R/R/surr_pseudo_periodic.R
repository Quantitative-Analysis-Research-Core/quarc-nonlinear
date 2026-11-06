#' Pseudo-periodic surrogate
#'
#' Small, Yu and Harrison (2001). A new series is built by hopping through
#' the reconstructed attractor: from the current point, the next is chosen
#' with probability falling off as `exp(-distance / rho)`, and its successor
#' is emitted.
#'
#' @param x Numeric vector.
#' @param delay,dim Embedding delay and dimension.
#' @param rho Noise radius. See [surr_find_rho()].
#' @return A list with `ys`, the surrogate, and `yi`, the index of `x` that
#'   each sample was taken from.
#' @references Small, M., Yu, D. and Harrison, R. G. (2001). Surrogate test
#'   for pseudoperiodic time series data. Physical Review Letters, 87(18),
#'   188101.
#' @export
surr_pseudo_periodic <- function(x, delay, dim, rho) {
  x <- as.numeric(x)
  Y <- psr(x, delay, dim)
  nY <- nrow(Y)
  xi <- floor(stats::runif(1) * nY) + 1
  ys <- numeric(nY)
  yi <- numeric(nY)
  ys[1] <- x[xi]
  yi[1] <- xi
  M <- nY - 2
  for (i in 2:nY) {
    prob <- exp(-sqrt(rowSums((Y[1:M, , drop = FALSE] -
                               matrix(Y[xi, ], M, ncol(Y), byrow = TRUE))^2)) / rho)
    sum3 <- cumsum(prob)
    xi_n <- NA
    ind <- 1
    while (is.na(xi_n) || xi_n > nY || xi_n == xi) {
      hit <- which(sum3 < stats::runif(1) * sum3[M])
      xi_n <- if (length(hit)) max(hit) + 2 else NA
      ind <- ind + 1
      if (ind == 100) {
        stop("a new value of xi could not be found, check that rho is not too low", call. = FALSE)
      }
    }
    xi <- xi_n
    ys[i] <- x[xi]
    yi[i] <- xi
  }
  list(ys = ys, yi = yi)
}
