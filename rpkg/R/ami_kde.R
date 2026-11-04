#' @describeIn ami Gaussian kernel density estimator (Thomas et al. 2014).
#'   Lag `k` uses all `N - k` pairs. Cost is O(maxlag * N^2).
#' @param chunk Rows of the pairwise kernel or distance matrix evaluated at
#'   once. Changes peak memory, not the result. Default 512.
#' @param bandwidth `"fixed"` (default) or `"adaptive"`: each kernel scaled
#'   by `sqrt(g / pilot_i)`, `g` the geometric mean of a fixed-bandwidth
#'   pilot (Abramson 1982). About twice the cost.
#' @export
ami_kde <- function(x, maxlag, fraction = 0.2, chunk = 512,
                    bandwidth = c("fixed", "adaptive")) {
  bandwidth <- match.arg(bandwidth)
  adaptive <- bandwidth == "adaptive"
  x <- check_ami_input(x, maxlag, "ami_kde")
  N <- length(x)
  curve <- cbind(lag = 0:maxlag, ami = 0)
  curve[1, 2] <- kde_mi_bits(x, x, chunk, adaptive)
  for (lag in seq_len(maxlag)) {
    curve[lag + 1, 2] <- kde_mi_bits(x[1:(N - lag)], x[(lag + 1):N], chunk, adaptive)
  }
  out <- ami_first_minimum(curve, fraction, "ami_kde")
  out$info$bandwidth <- bandwidth
  out$info$estimator <- if (adaptive) "kde-adaptive" else "kde"
  out
}

kde_mi_bits <- function(X, Y, chunk, adaptive = FALSE) {
  n <- length(X)
  hx <- stats::sd(X) / n^(1 / 6)
  hy <- stats::sd(Y) / n^(1 / 6)
  if (hx == 0 || hy == 0) return(0)
  pX <- kde_marginal(X, hx, chunk)
  pY <- kde_marginal(Y, hy, chunk)
  pXY <- kde_joint(X, Y, hx, hy, chunk)
  if (adaptive) {
    # The fixed-bandwidth estimates above are the pilots.
    pX <- kde_marginal(X, hx, chunk, abramson(pX))
    pY <- kde_marginal(Y, hy, chunk, abramson(pY))
    pXY <- kde_joint(X, Y, hx, hy, chunk, abramson(pXY))
  }
  ok <- pXY > 0 & pX > 0 & pY > 0
  sum(log2(pXY[ok] / (pX[ok] * pY[ok]))) / n
}

kde_marginal <- function(v, h, chunk, lambda = NULL) {
  n <- length(v)
  p <- numeric(n)
  cst <- 1 / (n * h * sqrt(2 * pi))
  for (lo in seq(1, n, by = chunk)) {
    rows <- lo:min(lo + chunk - 1, n)
    z <- outer(v[rows], v, `-`) / h
    if (is.null(lambda)) {
      p[rows] <- cst * rowSums(exp(-0.5 * z^2))
    } else {
      L <- matrix(lambda, length(rows), n, byrow = TRUE)
      p[rows] <- cst * rowSums(exp(-0.5 * (z / L)^2) / L)
    }
  }
  p
}

kde_joint <- function(X, Y, hx, hy, chunk, lambda = NULL) {
  n <- length(X)
  r <- sum((X - mean(X)) * (Y - mean(Y))) /
    (sqrt(sum((X - mean(X))^2)) * sqrt(sum((Y - mean(Y))^2)))
  r <- min(max(r, -0.999999), 0.999999)
  omr2 <- 1 - r^2
  cst <- 1 / (n * 2 * pi * hx * hy * sqrt(omr2))
  p <- numeric(n)
  for (lo in seq(1, n, by = chunk)) {
    rows <- lo:min(lo + chunk - 1, n)
    u <- outer(X[rows], X, `-`) / hx
    v <- outer(Y[rows], Y, `-`) / hy
    q <- (u^2 - 2 * r * (u * v) + v^2) / omr2
    if (is.null(lambda)) {
      p[rows] <- cst * rowSums(exp(-0.5 * q))
    } else {
      L2 <- matrix(lambda^2, length(rows), n, byrow = TRUE)
      p[rows] <- cst * rowSums(exp(-0.5 * q / L2) / L2)
    }
  }
  p
}

# Abramson's square-root law: lambda_i = sqrt(g / pilot_i).
abramson <- function(pilot) {
  pilot <- pmax(pilot, .Machine$double.xmin)
  g <- exp(mean(log(pilot)))
  sqrt(g / pilot)
}
