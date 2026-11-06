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
