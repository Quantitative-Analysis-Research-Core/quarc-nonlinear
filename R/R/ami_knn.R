#' @describeIn ami Kraskov-Stoegbauer-Grassberger estimator, algorithm 1.
#'   The series is standardised and moved by `1e-10 * (((1:N) * g) %% 1 -
#'   0.5)`, `g` the golden-ratio conjugate, to break ties deterministically,
#'   exactly as in MATLAB. Every lag uses `N - maxlag` pairs. Cost
#'   O(maxlag * N^2). Can return small negative values where the true AMI
#'   is zero.
#' @export
ami_knn <- function(x, maxlag, k = 10, fraction = 0.2, chunk = 512) {
  x <- check_ami_input(x, maxlag, "ami_knn")
  N <- length(x)
  m <- N - maxlag
  if (k >= m) {
    stop(sprintf("ami_knn: k (%d) must be smaller than the number of pairs per lag (%d).",
                 as.integer(k), m), call. = FALSE)
  }
  s <- stats::sd(x)
  if (s > 0) x <- (x - mean(x)) / s
  x <- x + 1e-10 * (((1:N) * 0.6180339887498949) %% 1 - 0.5)
  curve <- cbind(lag = 0:maxlag, ami = 0)
  for (lag in 0:maxlag) {
    curve[lag + 1, 2] <- ksg_bits(x[1:m], x[(1 + lag):(m + lag)], k, chunk)
  }
  out <- ami_first_minimum(curve, fraction, "ami_knn")
  out$info$samplesPerLag <- m
  out$info$k <- k
  out$info$estimator <- "knn"
  out
}

# Kraskov-Stoegbauer-Grassberger algorithm 1, in bits.
ksg_bits <- function(a, b, k, chunk) {
  n <- length(a)
  nx <- numeric(n)
  ny <- numeric(n)
  for (lo in seq(1, n, by = chunk)) {
    rows <- lo:min(lo + chunk - 1, n)
    dx <- abs(outer(a[rows], a, `-`))
    dy <- abs(outer(b[rows], b, `-`))
    d <- pmax(dx, dy)
    d[cbind(seq_along(rows), rows)] <- Inf          # exclude self
    e <- apply(d, 1, function(r) sort.int(r, partial = k)[k])
    # Strictly closer than eps_i; the self term is counted and removed.
    nx[rows] <- rowSums(dx < e) - 1
    ny[rows] <- rowSums(dy < e) - 1
  }
  (digamma(k) + digamma(n) - mean(digamma(nx + 1) + digamma(ny + 1))) / log(2)
}
