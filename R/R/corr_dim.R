#' Correlation dimension
#'
#' Grassberger and Procaccia (1983). The correlation sum is computed over
#' 200 log-spaced radii, excluding pairs closer than `dim * delay` samples in
#' time, and the slope is fitted over the band of `ln C` from the midpoint of
#' its range to a quarter-range above it.
#'
#' The region is selected by height, not linearity: check `extra$r2`. A low
#' value means the number is a slope through a bend, not a scaling exponent.
#'
#' @param x Numeric vector.
#' @param delay Embedding delay in samples.
#' @param dim Embedding dimension.
#' @return A list with `CoD`, the correlation dimension, and `extra`: `logEps`,
#'   `logC`, `idx` (the points fitted), `slope`, `intercept`, `r2` and
#'   `bounds`.
#' @references Grassberger, P. and Procaccia, I. (1983). Measuring the
#'   strangeness of strange attractors. Physica D, 9(1-2), 189-208.
#' @export
corr_dim <- function(x, delay, dim) {
  x <- as.numeric(x)
  bins <- 200
  n <- length(x) - (dim - 1) * delay
  y <- embed(x, dim, delay)
  k <- dim * delay
  pair_dist <- function(i) {
    j <- (i + k + 1):n
    sqrt(colSums((y[, j, drop = FALSE] - y[, i])^2))
  }
  eps1 <- 0
  eps2 <- Inf
  for (i in seq_len(n - k - 1)) {
    d <- pair_dist(i)
    eps1 <- max(max(d), eps1)
    eps2 <- min(min(d), eps2)
  }
  if (eps2 == 0) eps2 <- .Machine$double.eps
  lo <- log(eps2)
  hi <- log(eps1)
  e <- lo + (0:(bins - 1)) * (hi - lo) / (bins - 1)
  e[bins] <- hi
  epsilon <- exp(1)^e
  epsilon[c(1, bins)] <- c(eps2, eps1)   # exact end radii; see the MATLAB note
  CI <- numeric(bins)
  for (i in seq_len(n - k - 1)) {
    d <- pair_dist(i)
    # histc: edges[b] <= d < edges[b+1]; the last bin counts d == edges[end].
    b <- findInterval(d, epsilon)
    b[d == epsilon[bins]] <- bins
    b <- b[b >= 1 & (b < bins | d == epsilon[bins])]
    CI <- CI + tabulate(b, bins)
  }
  CI <- log(cumsum(CI) / ((n - k)^2))
  CI <- CI + log(n - k - 1)
  epsilon <- log(epsilon)
  if (is.infinite(min(CI))) {
    i <- max(which(CI == -Inf))
    minCI <- min(CI[(i + 1):bins])
  } else {
    minCI <- min(CI)
  }
  mid <- (max(CI) - minCI) / 2 + minCI
  quarter <- (max(CI) - minCI) / 4
  bounds <- c(mid, mid + quarter)
  intervals <- which(CI > bounds[1] & CI < bounds[2])
  idx <- min(intervals):max(intervals)
  fit <- stats::lm.fit(cbind(epsilon[idx], 1), CI[idx])
  cf <- unname(fit$coefficients)
  yHat <- cf[1] * epsilon[idx] + cf[2]
  ssTot <- sum((CI[idx] - mean(CI[idx]))^2)
  r2 <- if (ssTot > 0 && length(idx) > 2) 1 - sum((CI[idx] - yHat)^2) / ssTot else NA_real_
  list(CoD = cf[1],
       extra = list(logEps = epsilon, logC = CI, idx = idx, slope = cf[1], intercept = cf[2],
                    r2 = r2, bounds = bounds))
}
