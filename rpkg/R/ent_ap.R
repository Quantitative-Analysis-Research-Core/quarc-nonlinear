#' Approximate entropy
#'
#' Pincus (1991). Self-matches are counted, as in the original definition.
#'
#' @param x Numeric vector.
#' @param dim Template length (`m` in the paper).
#' @param radius Tolerance as a proportion of `sd(x)` (`r` in the paper).
#' @return Approximate entropy in nats.
#' @references Pincus, S. M. (1991). Approximate entropy as a measure of
#'   system complexity. PNAS, 88(6), 2297-2301.
#' @export
ent_ap <- function(x, dim, radius) {
  x <- as.numeric(x)
  radius <- radius * stats::sd(x)
  N <- length(x)
  phim <- numeric(2)
  for (j in 1:2) {
    m <- dim + j - 1
    n <- N - m + 1
    phi <- numeric(n)
    for (i in seq_len(n)) {
      d <- abs(x[seq_len(n)] - x[i])
      for (k in seq_len(m - 1)) {
        d <- pmax(d, abs(x[seq_len(n) + k] - x[i + k]))
      }
      phi[i] <- sum(d <= radius) / n
    }
    phim[j] <- sum(log(phi)) / n
  }
  phim[1] - phim[2]
}
