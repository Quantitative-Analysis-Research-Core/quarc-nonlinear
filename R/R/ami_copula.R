#' @describeIn ami Gaussian copula estimator (Ince et al. 2017): ranks
#'   mapped to normal scores, then `-0.5 * log2(1 - r^2)`. Exact for a
#'   Gaussian copula and a lower bound otherwise. It sees only dependence a
#'   correlation of normal scores can see: on a sine at a quarter-period lag
#'   it is near zero although the true AMI is large. `r` is capped at
#'   `1 - 1e-6`. Every lag uses `N - maxlag` pairs.
#' @export
ami_copula <- function(x, maxlag, fraction = 0.2) {
  x <- check_ami_input(x, maxlag, "ami_copula")
  N <- length(x)
  m <- N - maxlag
  curve <- cbind(lag = 0:maxlag, ami = 0)
  for (lag in 0:maxlag) {
    a <- normal_scores(x[1:m])
    b <- normal_scores(x[(1 + lag):(m + lag)])
    r <- sum(a * b) / (sqrt(sum(a^2)) * sqrt(sum(b^2)))
    r <- min(max(r, -0.999999), 0.999999)
    curve[lag + 1, 2] <- -0.5 * log2(1 - r^2)
  }
  out <- ami_first_minimum(curve, fraction, "ami_copula")
  out$info$samplesPerLag <- m
  out$info$estimator <- "copula"
  out
}

# Ranks mapped to standard normal quantiles.
normal_scores <- function(v) {
  n <- length(v)
  u <- numeric(n)
  u[order(v)] <- seq_len(n) / (n + 1)
  stats::qnorm(u)
}
