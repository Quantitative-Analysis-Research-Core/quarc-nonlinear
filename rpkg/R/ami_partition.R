#' @describeIn ami Adaptive partition of the rank plane (Fraser and Swinney
#'   1986; Darbellay and Vajda 1999). A cell is cut into four at its median
#'   ranks when a chi-square test (3 df) rejects a uniform spread; the root
#'   is always cut. Every lag uses `N - maxlag` pairs.
#' @param alpha Significance level of the splitting test. Default 0.05.
#' @export
ami_partition <- function(x, maxlag, alpha = 0.05, fraction = 0.2) {
  x <- check_ami_input(x, maxlag, "ami_partition")
  if (!(alpha > 0 && alpha < 1)) stop("alpha must lie in (0, 1).", call. = FALSE)
  N <- length(x)
  m <- N - maxlag
  critical <- stats::qchisq(1 - alpha, 3)
  curve <- cbind(lag = 0:maxlag, ami = 0)
  for (lag in 0:maxlag) {
    curve[lag + 1, 2] <- partition_bits(x[1:m], x[(1 + lag):(m + lag)], critical)
  }
  out <- ami_first_minimum(curve, fraction, "ami_partition")
  out$info$samplesPerLag <- m
  out$info$alpha <- alpha
  out$info$estimator <- "adaptive"
  out
}

# Darbellay-Vajda estimate of I(A;B) in bits, on ranks.
partition_bits <- function(a, b, critical) {
  n <- length(a)
  rx <- rank_of(a)
  ry <- rank_of(b)
  # Explicit stack of cells: bounds [xlo xhi ylo yhi] in ranks, members.
  stackB <- list(c(1, n, 1, n))
  stackM <- list(seq_len(n))
  isRoot <- TRUE
  I <- 0
  while (length(stackB)) {
    top <- length(stackB)
    cb <- stackB[[top]]
    idx <- stackM[[top]]
    stackB[[top]] <- NULL
    stackM[[top]] <- NULL
    nc <- length(idx)
    if (nc == 0) next
    wx <- cb[2] - cb[1] + 1
    wy <- cb[4] - cb[3] + 1
    split <- FALSE
    if (wx >= 2 && wy >= 2) {
      xm <- cb[1] + floor(wx / 2) - 1
      ym <- cb[3] + floor(wy / 2) - 1
      lx <- rx[idx] <= xm
      ly <- ry[idx] <= ym
      q <- c(sum(lx & ly), sum(lx & !ly), sum(!lx & ly), sum(!lx & !ly))
      fx <- (xm - cb[1] + 1) / wx
      fy <- (ym - cb[3] + 1) / wy
      E <- nc * c(fx * fy, fx * (1 - fy), (1 - fx) * fy, (1 - fx) * (1 - fy))
      split <- isRoot || sum((q - E)^2 / E) > critical
    }
    isRoot <- FALSE
    if (split) {
      # Pushed in MATLAB's order so the summation order matches.
      stackB <- c(stackB, list(c(cb[1], xm, cb[3], ym), c(cb[1], xm, ym + 1, cb[4]),
                               c(xm + 1, cb[2], cb[3], ym), c(xm + 1, cb[2], ym + 1, cb[4])))
      stackM <- c(stackM, list(idx[lx & ly], idx[lx & !ly], idx[!lx & ly], idx[!lx & !ly]))
    } else {
      I <- I + (nc / n) * log2(nc * n / (wx * wy))
    }
  }
  I
}
