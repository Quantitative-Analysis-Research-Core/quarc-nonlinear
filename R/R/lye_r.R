#' Largest Lyapunov exponent, Rosenstein's method
#'
#' Rosenstein, Collins and De Luca (1993): each point's nearest neighbour
#' outside a Theiler window, and the mean log divergence of the pairs over
#' time.
#'
#' @param x Numeric vector, or a phase space (one coordinate per column).
#' @param fs Sampling frequency, Hz.
#' @param delay,dim Embedding delay and dimension, used when `x` is a vector.
#' @param theiler Temporal exclusion in samples. Default `round(0.8 * delay)`.
#' @param slope Optional `c(shortStart, shortEnd, longStart, longEnd)` in
#'   periods. With `MeanPeriod`, fits short- and long-range slopes.
#' @param MeanPeriod Mean period in seconds, used with `slope`.
#' @return Without `slope`, the `M`-by-3 matrix `out`: index, matched
#'   neighbour and mean log divergence at each step. With `slope`, a list
#'   with `LyES`, `LyEL` and `out`.
#' @references Rosenstein, M. T., Collins, J. J. and De Luca, C. J. (1993). A
#'   practical method for calculating largest Lyapunov exponents from small
#'   data sets. Physica D, 65(1-2), 117-134.
#' @export
lye_r <- function(x, fs, delay, dim, theiler = NULL, slope = NULL, MeanPeriod = NULL) {
  if (is.null(theiler)) theiler <- round(delay * 0.8)
  if (length(theiler) != 1 || !is.finite(theiler) || theiler < 0) {
    stop("TheilerWindow must be a non-negative finite scalar, in samples.", call. = FALSE)
  }
  theiler <- round(theiler)
  x <- as.matrix(x)
  if (ncol(x) > nrow(x)) x <- t(x)
  if (ncol(x) > 1) {
    Y <- x
  } else {
    Y <- psr(x, delay, dim)
  }
  M <- nrow(Y)
  core <- .lye_r_core(Y, theiler)
  out <- cbind(seq_len(M), core$nn, core$avg)
  AveLnDiv <- core$avg[seq_len(core$count)]
  dimnames(out) <- NULL
  if (is.null(slope)) return(out)

  time <- (0:(length(AveLnDiv) - 1)) / fs / MeanPeriod
  short <- c(if (slope[1] == 0) 1 else round(slope[1] * MeanPeriod * fs),
             round(slope[2] * MeanPeriod * fs))
  long <- c(round(slope[3] * MeanPeriod * fs), round(slope[4] * MeanPeriod * fs))
  fit_slope <- function(r) {
    if (r[2] > length(AveLnDiv)) return(NaN)
    idx <- r[1]:r[2]
    unname(stats::lm.fit(cbind(time[idx], 1), AveLnDiv[idx])$coefficients[1])
  }
  list(LyES = fit_slope(short), LyEL = fit_slope(long), out = out)
}

# The R code for the two loops, kept as the readable reference for the
# compiled kernel (cpp/src/lye_r.cpp); the tests check that the two agree.
# Returns out (index, neighbour, mean log divergence) and AveLnDiv.
lye_r_divergence_r <- function(Y, theiler) {
  M <- nrow(Y)
  IND2 <- integer(M)
  for (i in seq_len(M)) {
    d <- sqrt(rowSums((Y - matrix(Y[i, ], M, ncol(Y), byrow = TRUE))^2))
    ex <- (i - theiler):(i + theiler)
    d[ex[ex >= 1 & ex <= M]] <- Inf
    IND2[i] <- which.min(d)
  }
  out <- cbind(seq_len(M), IND2, 0)
  sumLog <- numeric(M - 1)
  nDist <- numeric(M - 1)
  for (i in seq_len(M - 1)) {
    endITL <- min(M - IND2[i], M - i)
    if (endITL < 1) next
    d <- sqrt(rowSums((Y[(i + 1):(endITL + i), , drop = FALSE] -
                       Y[(IND2[i] + 1):(endITL + IND2[i]), , drop = FALSE])^2))
    pos <- which(d > 0)
    sumLog[pos] <- sumLog[pos] + log(d[pos])
    nDist[pos] <- nDist[pos] + 1
  }
  AveLnDiv <- sumLog / nDist
  ok <- nDist > 0
  out[which(ok), 3] <- AveLnDiv[ok]
  # MATLAB grows AveLnDiv only up to the last step with a pair.
  AveLnDiv <- AveLnDiv[seq_len(max(c(0, which(ok))))]
  dimnames(out) <- NULL
  list(out = out, AveLnDiv = AveLnDiv)
}
