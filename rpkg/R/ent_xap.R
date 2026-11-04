#' Cross approximate entropy
#'
#' Both series are standardised. Templates are similar when every component
#' differs by strictly less than `radius`.
#'
#' @param x,y Numeric vectors of equal length.
#' @param dim Template length.
#' @param radius Tolerance, in standard deviations.
#' @param k Step between the two template lengths compared, `dim` and
#'   `dim + k`.
#' @return Cross approximate entropy in nats.
#' @export
ent_xap <- function(x, y, dim, radius, k) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  x <- (x - mean(x)) / stats::sd(x)
  y <- (y - mean(y)) / stats::sd(y)
  N <- length(x)
  radius <- stats::sd(x) * radius
  Cm <- numeric(2)
  for (j in 1:2) {
    m <- dim + (j - 1) * k
    n <- N - m + 1
    C <- numeric(n)
    for (i in seq_len(n)) {
      ok <- abs(x[i] - y[seq_len(n)]) < radius
      for (q in seq_len(m - 1)) ok <- ok & (abs(x[i + q] - y[seq_len(n) + q]) < radius)
      C[i] <- sum(ok) / n
    }
    Cm[j] <- sum(C) / n
  }
  log(Cm[1] / Cm[2])
}
