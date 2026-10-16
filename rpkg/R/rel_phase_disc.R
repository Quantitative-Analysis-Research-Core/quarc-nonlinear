#' Discrete relative phase, frequency ratios and phase coupling
#'
#' Event-based relative phase between two oscillators (by default their
#' peaks), the distribution of frequency ratios it implies, and the
#' strength of phase coupling at the dominant ratio. This is the
#' computation of MATLAB's `rel_phase_disc` without its figures.
#'
#' @param x,y Numeric vectors.
#' @param S1_events,S2_events Event indices for `x` and `y`. Default: their
#'   peaks, found as MATLAB's `findpeaks` does.
#' @return A list of three data frames: `FR` (`Ratio`, `Occurances`,
#'   `Percentage_of_Total`), `PC` (`n`, `PC`) and `DRP` (`t`, `n`, `T`,
#'   `DRP`, one row per event, DRP in degrees).
#' @export
rel_phase_disc <- function(x, y, S1_events = NULL, S2_events = NULL) {
  if (is.null(S1_events)) S1_events <- find_peaks(x)
  if (is.null(S2_events)) S2_events <- find_peaks(y)
  S1 <- S1_events
  S2 <- S2_events
  if (length(S1) < 2) stop("x needs at least two events.", call. = FALSE)
  # MATLAB grows an array by assignment past its end, filling with zeros.
  put <- function(v, i, val) {
    if (i > length(v)) v <- c(v, numeric(i - length(v)))
    v[i] <- val
    v
  }
  mround <- function(v) sign(v) * floor(abs(v) + 0.5)

  i <- 2
  back <- 0
  tiRp <- numeric(0)
  for (i1 in seq_along(S2)) {
    if (S2[i1] < S1[1]) {
      back <- back + 1
    } else if (S2[i1] < S1[i]) {
      tiRp <- put(tiRp, i1 - back, S1[i] - S2[i1])
    } else if (i < length(S1)) {
      i <- i + 1
      tiRp <- put(tiRp, i1 - back, S1[i] - S2[i1])
    }
  }
  Ti <- diff(S2)
  N <- vapply(seq_len(length(S1) - 1), function(k) {
    f <- S2[seq_len(length(S2) - 1)]
    sum(f >= S1[k] & f < S1[k + 1]) - 1
  }, numeric(1))
  L <- min(length(tiRp), length(Ti))
  DRP <- tiRp[seq_len(L)] / Ti[seq_len(L)] * 360
  tcol <- tiRp[seq_len(L)]
  Tcol <- Ti[seq_len(L)]
  n <- unlist(lapply(N, function(Nk) if (Nk >= 0) Nk - (0:Nk) else numeric(0)))
  if (length(n) != L) {
    stop(sprintf("All table variables must have the same number of rows (n has %d, DRP %d).",
                 length(n), L), call. = FALSE)
  }
  Raw <- data.frame(t = tcol, n = n, T = Tcol, DRP = DRP)
  maxDRP <- mround(max(DRP) / 360) + 1

  DRPfirst <- DRP[1]
  if (nrow(Raw) - 1 >= 3) {
    for (k in 3:(nrow(Raw) - 1)) if (Raw$n[k] == 0) DRPfirst <- c(DRPfirst, Raw$DRP[k + 1])
  }
  freq <- numeric(0)
  for (k in seq_along(DRPfirst)) {
    for (i1 in seq_len(maxDRP)) {
      if (DRPfirst[k] > (i1 - 1) * 360 && DRPfirst[k] <= i1 * 360) freq <- put(freq, k, i1)
    }
  }
  if (length(freq) < 4) stop("Too few cycles to classify frequency ratios.", call. = FALSE)
  half <- numeric(0)
  half[1] <- if (freq[1] == freq[3] && freq[2] == freq[4] && abs(freq[1] - freq[2]) == 1) {
    (freq[1] + freq[2]) / 2
  } else {
    freq[1]
  }
  last <- length(freq)
  counter <- 0
  for (k in 2:(last - 1)) {
    if (freq[k - 1] == freq[k + 1] && abs(freq[k] - freq[k + 1]) == 1) {
      half <- put(half, k, (freq[k] + freq[k + 1]) / 2)
      counter <- counter + 1
    } else {
      if (counter %% 2 == 1 && counter < k && counter > 1) {
        half <- put(half, k - counter - 1, half[k - counter])
      }
      half <- put(half, k, freq[k])
      counter <- 0
    }
  }
  half <- put(half, last, freq[last])
  nh <- length(half)
  withNC <- numeric(nh)
  withNC[1] <- if (half[1] == half[2]) half[1] else 0
  for (k in 2:(nh - 1)) {
    withNC[k] <- if (half[k] == half[k + 1] || half[k] == half[k - 1]) half[k] else 0
  }
  withNC[nh] <- if (half[nh] == half[nh - 1]) half[nh] else 0

  maxn <- max(withNC)
  Ratio <- 0
  Occ <- sum(withNC == 0)
  Pct <- Occ / nh * 100
  if (maxn >= 1) {
    for (r in seq(1, maxn, by = 0.5)) {
      Ratio <- put(Ratio, r * 2, r)
      Occ <- put(Occ, r * 2, sum(withNC == r))
      Pct <- put(Pct, r * 2, Occ[r * 2] / nh * 100)
    }
  }
  lag <- which.max(Occ)
  if (lag %% 2 == 0) lag <- lag / 2
  original <- DRP[seq_len(length(DRP) - lag)]
  lagged <- DRP[(1 + lag):length(DRP)]
  kmax <- mround(Ratio[Pct == max(Pct)][1])
  wdn <- numeric(0)
  nPC <- PC <- numeric(0)
  for (k in seq_len(kmax)) {
    sel <- original <= 360 * k & lagged <= 360 * k
    dn <- (original[sel] - lagged[sel]) * cos(pi / 4)
    for (q in seq_along(dn)) wdn <- put(wdn, q, if (abs(dn[q]) <= 40) 1 - abs(dn[q]) / 40 else 0)
    nPC[k] <- kmax - k
    PC[k] <- sum(wdn) / length(wdn) * 100
  }
  list(FR = data.frame(Ratio = Ratio, Occurances = Occ, Percentage_of_Total = Pct),
       PC = data.frame(n = nPC, PC = PC),
       DRP = Raw)
}

# Indices of local maxima, as MATLAB's findpeaks: strictly above the
# preceding sample, and for a flat top the first sample of the plateau,
# provided it is followed by a fall. Endpoints are never peaks.
find_peaks <- function(x) {
  x <- as.numeric(x)
  n <- length(x)
  peaks <- integer(0)
  i <- 2
  while (i < n) {
    if (x[i] > x[i - 1]) {
      j <- i
      while (j < n && x[j + 1] == x[i]) j <- j + 1
      if (j < n && x[j + 1] < x[i]) peaks <- c(peaks, i)
      i <- j + 1
    } else {
      i <- i + 1
    }
  }
  peaks
}
