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

#' Largest Lyapunov exponent, Wolf's method
#'
#' Wolf, Swift, Swinney and Vastano (1985): follow a reference trajectory
#' and a neighbour for `evolve` samples, record the log stretch, then replace
#' the neighbour with one near the evolved point and close in direction.
#'
#' @param x Numeric vector, or a phase space.
#' @param fs Sampling frequency, Hz.
#' @param delay,dim Embedding delay and dimension, used when `x` is a vector.
#' @param evolve Propagation length in samples. No default: use 1 for a map.
#' @param SCALEMX,ANGLMX,ZMULT Maximum replacement distance (default a tenth
#'   of the data range), maximum replacement angle in radians (default 30
#'   degrees) and Wolf's scale multiplier. Give all three or none.
#' @return A list with `out`, one row per replacement step (iteration,
#'   point, neighbour, start and end distance, running exponent, current
#'   scale, best angle and maximum angle in degrees), and `LyE`, the exponent
#'   in bits per unit time.
#' @references Wolf, A., Swift, J. B., Swinney, H. L. and Vastano, J. A.
#'   (1985). Determining Lyapunov exponents from a time series. Physica D,
#'   16(3), 285-317.
#' @export
lye_w <- function(x, fs, delay, dim, evolve, SCALEMX = NULL, ANGLMX = NULL, ZMULT = NULL) {
  given <- !c(is.null(SCALEMX), is.null(ANGLMX), is.null(ZMULT))
  if (!any(given)) {
    SCALEMX <- (max(x) - min(x)) / 10
    ANGLMX <- 30 * pi / 180
    ZMULT <- 1
  } else if (!all(given)) {
    stop("not enough input arguements", call. = FALSE)
  }
  DT <- 1 / fs
  x <- as.matrix(x)
  if (ncol(x) == 1) {
    Y <- psr(x, delay, dim)
    NPT <- nrow(x) - (dim - 1) * delay - evolve
    Y <- Y[1:(NPT + evolve), , drop = FALSE]
  } else {
    Y <- x
    NPT <- nrow(Y) - evolve
  }
  out <- matrix(0, floor(NPT / evolve), 9)
  ITS <- 0
  distSUM <- 0
  thbest <- 0
  OUTMX <- SCALEMX
  pair <- NA_integer_
  dist_to <- function(p) sqrt(rowSums((Y[1:NPT, , drop = FALSE] -
                                       matrix(Y[p, ], NPT, ncol(Y), byrow = TRUE))^2))
  minIndex <- function(v) { w <- which.min(v); if (length(w)) w else 1L }
  LyE <- NA_real_
  for (i in seq(1, NPT, by = evolve)) {
    cp <- i
    if (cp == 1) {
      d <- dist_to(cp)
      ex <- (cp - 10):(cp + 10)
      d[d <= 0] <- NA
      d[ex[ex >= 1 & ex <= NPT]] <- NA
      pair <- minIndex(d)
    }
    start_dist <- sqrt(sum((Y[cp, ] - Y[pair, ])^2))
    end_dist <- sqrt(sum((Y[cp + evolve, ] - Y[pair + evolve, ])^2))
    distSUM <- distSUM + log2(end_dist / start_dist) / (evolve * DT)
    ITS <- ITS + 1
    LyE <- distSUM / ITS
    row <- floor(i / evolve) + 1
    if (row > nrow(out)) out <- rbind(out, matrix(0, row - nrow(out), 9))
    out[row, ] <- c(ITS, cp, pair, start_dist, end_dist, LyE, OUTMX, thbest * 180 / pi,
                    ANGLMX * 180 / pi)
    ZMULT <- 1
    if (end_dist < SCALEMX) {
      pair <- pair + evolve
      if (pair > NPT) {
        pair <- pair - evolve
        nx <- wolf_next_point(1, Y, cp, pair, NPT, evolve, SCALEMX, ANGLMX)
        pair <- nx$next_point; ANGLMX <- nx$ANGLMX; thbest <- nx$thbest; OUTMX <- nx$SCALEMX
      }
      next
    }
    nx <- wolf_next_point(0, Y, cp, pair, NPT, evolve, SCALEMX, ANGLMX)
    pair <- nx$next_point; ANGLMX <- nx$ANGLMX; thbest <- nx$thbest; OUTMX <- nx$SCALEMX
  }
  list(out = out, LyE = LyE)
}

# Replacement neighbour near the evolved point, preferring a small angle to
# the current separation (GetNextPoint and find_next_point in MATLAB).
wolf_next_point <- function(flag, Y, cp, pair, NPT, evolve, SCALEMX, ANGLMX) {
  p1 <- Y[cp + evolve, ]
  Vnew <- matrix(p1, NPT, ncol(Y), byrow = TRUE) - Y[1:NPT, , drop = FALSE]
  d <- sqrt(rowSums(Vnew^2))
  ex <- (cp + evolve - 10):(cp + evolve + 10)
  d[ex[ex >= 1 & ex <= NPT]] <- NA
  Vcurr <- p1 - Y[pair + evolve, ]
  end_dist <- sqrt(sum(Vcurr^2))
  cosTheta <- abs(as.numeric(Vnew %*% Vcurr) / (d * end_dist))
  # MATLAB's acos of a cosine rounded just above 1 is complex: comparisons
  # use its real part (0) and abs() its modulus (acosh).
  above <- !is.na(cosTheta) & cosTheta > 1
  thetaRe <- ifelse(above, 0, acos(pmin(cosTheta, 1)))
  thetaAbs <- ifelse(above, acosh(pmax(cosTheta, 1)), thetaRe)
  minIndex <- function(v) { w <- which.min(v); if (length(w)) w else 1L }
  if (flag == 0) {
    pot <- d
    pot[!is.na(d) & (d <= 0 | thetaRe >= ANGLMX)] <- NA
    np <- minIndex(pot)
    if (!is.na(pot[np]) && pot[np] <= SCALEMX) {
      return(list(next_point = np, ANGLMX = 30 * pi / 180, thbest = abs(thetaAbs[np]),
                  SCALEMX = SCALEMX))
    }
  }
  pot <- d
  pot[!is.na(d) & d <= 0] <- NA
  list(next_point = minIndex(pot), ANGLMX = ANGLMX, thbest = ANGLMX, SCALEMX = SCALEMX)
}

#' Largest Lyapunov exponent
#'
#' One entry point for both methods, returning nats per unit time.
#'
#' @param x Numeric vector (embedded with `delay` and `dim`) or a phase space.
#' @param fs Sampling frequency, Hz. Use 1 for maps.
#' @param algorithm `"rosenstein"` (default; alias `"r"`) or `"wolf"`
#'   (alias `"w"`).
#' @param delay,dim Embedding delay and dimension. Defaults 1 and 3.
#' @param theiler Theiler window, `"rosenstein"` only. Default
#'   `round(0.8 * delay)`.
#' @param evolve Propagation length in samples, required for `"wolf"`: 1 for
#'   a map, about 10 for a flow at ~40 samples per orbit.
#' @return A list with `lambda`, in nats per unit time, and `extra`: the
#'   phase space and, for `"rosenstein"`, the divergence curve, matched
#'   pairs, fitted scaling region and its R^2; for `"wolf"`, the iteration
#'   table, the exponent in bits and `evolve`.
#' @export
lyapunov <- function(x, fs, algorithm = "rosenstein", delay = 1, dim = 3, theiler = NULL,
                     evolve = NULL) {
  if (anyNA(x)) stop("x contains NA. Remove or impute before estimating an exponent.", call. = FALSE)
  if (is.null(dim(x)) || ncol(as.matrix(x)) == 1) {
    x <- as.numeric(x)
    if (length(x) <= (dim - 1) * delay + 2) {
      stop(sprintf("Series of %d samples is too short for dim = %d at delay = %d.",
                   length(x), as.integer(dim), as.integer(delay)), call. = FALSE)
    }
    Y <- psr(x, delay, dim)
  } else {
    Y <- as.matrix(x)
    if (is.null(theiler) && delay == 1) {
      warning(paste0("A phase space was supplied without theiler or delay, so the temporal ",
                     "exclusion defaults to 1 sample. That is almost certainly too small ",
                     "for an embedded series and will bias the exponent. Pass theiler (or ",
                     "delay) to match the delay used to build the phase space."), call. = FALSE)
    }
  }
  if (is.null(theiler)) theiler <- round(delay * 0.8)
  algorithm <- tolower(algorithm)
  extra <- list(phaseSpace = Y, algorithm = algorithm, units = "nats per unit time", fs = fs,
                theilerWindow = theiler, tau = delay, dim = ncol(Y))
  if (algorithm %in% c("rosenstein", "r")) {
    out <- lye_r(Y, fs, delay, ncol(Y), theiler = theiler)
    d <- out[, 3]
    sr <- scaling_region(d, fs)
    extra$divergence <- d
    extra$matchedPairs <- out[, 1:2]
    extra$scalingRegion <- sr$idx
    extra$fitR2 <- sr$r2
    extra$estimator <- "rosenstein"
    lambda <- sr$slope
  } else if (algorithm %in% c("wolf", "w")) {
    if (is.null(evolve)) {
      stop(paste0('algorithm="wolf" needs evolve, the propagation length in samples. There ',
                  "is no default: use evolve = 1 for a map, and about 10 for a flow sampled ",
                  "at ~40 samples per orbit."), call. = FALSE)
    }
    w <- lye_w(Y, fs, delay, ncol(Y), evolve)
    lambda <- w$LyE * log(2)
    extra$iterations <- w$out
    extra$bits <- w$LyE
    extra$estimator <- "wolf"
    extra$evolve <- evolve
  } else {
    stop(sprintf('Unknown algorithm "%s". Supported: "rosenstein" (alias "r"), "wolf" (alias "w").',
                 algorithm), call. = FALSE)
  }
  list(lambda = lambda, extra = extra)
}

# Longest near-linear window of the divergence curve before it saturates.
scaling_region <- function(d, fs) {
  d <- d[is.finite(d)]
  n <- length(d)
  if (n < 6) return(list(slope = NaN, idx = integer(0), r2 = NaN))
  sat <- saturation_index(d)
  hi <- max(6, min(sat, 200))
  line_r2 <- function(xv, yv) {
    if (length(xv) < 2 || all(yv == yv[1])) return(-Inf)
    ss <- sum((yv - mean(yv))^2)
    if (ss <= 0) return(-Inf)
    cf <- stats::lm.fit(cbind(xv, 1), yv)$coefficients
    1 - sum((yv - (cf[1] * xv + cf[2]))^2) / ss
  }
  best <- list(r2 = -Inf, a = 1, b = min(hi, 5), len = 0)
  cand <- list()
  for (a in seq_len(min(4, max(1, hi - 4)))) {
    if (a + 3 > hi) next
    for (b in (a + 3):hi) {
      cnd <- list(r2 = unname(line_r2(a:b, d[a:b])), a = a, b = b, len = b - a + 1)
      cand[[length(cand) + 1]] <- cnd
      if (cnd$r2 > best$r2) best <- cnd
    }
  }
  if (!length(cand)) return(list(slope = NaN, idx = integer(0), r2 = NaN))
  pick <- best
  for (cnd in cand) if (cnd$r2 >= best$r2 - 0.02 && cnd$len > pick$len) pick <- cnd
  idx <- pick$a:pick$b
  slope <- unname(stats::lm.fit(cbind(idx / fs, 1), d[idx])$coefficients[1])
  list(slope = slope, idx = idx, r2 = pick$r2)
}

saturation_index <- function(d) {
  n <- length(d)
  span <- max(d) - min(d)
  if (span <= 0) return(n)
  tol <- 0.01 * span
  look <- max(3, round(0.02 * n))
  runmax <- d[1]
  if (n - look >= 2) {
    for (i in 2:(n - look)) {
      runmax <- max(runmax, d[i])
      if (max(d[(i + 1):(i + look)]) <= runmax + tol) return(i)
    }
  }
  n
}
