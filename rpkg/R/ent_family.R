#' Cross sample entropy
#'
#' Richman and Moorman (2000): how often a template of `y` is matched by
#' templates of `x`, for lengths `dim` and `dim + 1`.
#'
#' @param x,y Numeric vectors of equal length.
#' @param dim Template length.
#' @param radius Tolerance.
#' @param norm Normalisation: `1` rescales both series to the range 0-1 and
#'   sets the tolerance to `radius` times the mean of their SDs; `2`
#'   standardises both and uses `radius` as is; anything else leaves the data
#'   alone and uses `radius` times the mean SD of the raw series.
#' @return Cross sample entropy in nats.
#' @export
ent_xsamp <- function(x, y, dim, radius, norm) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  if (length(x) != length(y)) stop("The data series need to be the same length!", call. = FALSE)
  N <- length(x)
  if (norm == 1) {
    xn <- (x - min(x)) / (max(x) - min(x))
    yn <- (y - min(y)) / (max(y) - min(y))
    r <- radius * ((stats::sd(xn) + stats::sd(yn)) / 2)
  } else if (norm == 2) {
    xn <- (x - mean(x)) / stats::sd(x)
    yn <- (y - mean(y)) / stats::sd(y)
    r <- radius
  } else {
    xn <- x
    yn <- y
    r <- radius * ((stats::sd(x) + stats::sd(y)) / 2)
  }
  n <- N - dim
  Bm <- numeric(n)
  Am <- numeric(n)
  for (i in seq_len(n)) {
    dj <- abs(xn[seq_len(n)] - yn[i])
    for (k in seq_len(dim - 1)) dj <- pmax(dj, abs(xn[seq_len(n) + k] - yn[i + k]))
    dj1 <- pmax(dj, abs(xn[seq_len(n) + dim] - yn[i + dim]))
    Bm[i] <- sum(dj <= r) / n
    Am[i] <- sum(dj1 <= r) / n
  }
  -log((sum(Am) / n) / (sum(Bm) / n))
}

#' Cross approximate entropy
#'
#' Both series are standardised. Templates are similar when every component
#' differs by strictly less than `radius`.
#'
#' @param x,y Numeric vectors of equal length.
#' @param dim Template length.
#' @param radius Tolerance, in standard deviations.
#' @param k Step between the two template lengths compared, `dim` and
#'   `dim + k`.
#' @return Cross approximate entropy in nats.
#' @export
ent_xap <- function(x, y, dim, radius, k) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  x <- (x - mean(x)) / stats::sd(x)
  y <- (y - mean(y)) / stats::sd(y)
  N <- length(x)
  radius <- stats::sd(x) * radius
  Cm <- numeric(2)
  for (j in 1:2) {
    m <- dim + (j - 1) * k
    n <- N - m + 1
    C <- numeric(n)
    for (i in seq_len(n)) {
      ok <- abs(x[i] - y[seq_len(n)]) < radius
      for (q in seq_len(m - 1)) ok <- ok & (abs(x[i + q] - y[seq_len(n) + q]) < radius)
      C[i] <- sum(ok) / n
    }
    Cm[j] <- sum(C) / n
  }
  log(Cm[1] / Cm[2])
}

#' Normalised corrected Shannon entropy of a binary series
#'
#' Words of `dim` consecutive symbols are read as binary numbers; their
#' Shannon entropy, corrected for the number of distinct words seen, is
#' divided by its maximum.
#'
#' @param x Numeric vector of 0s and 1s.
#' @param dim Word length.
#' @return Normalised corrected entropy, between 0 and about 1.
#' @export
ent_symbolic <- function(x, dim) {
  x <- as.numeric(x)
  u <- unique(x)
  if (!all(u %in% c(0, 1))) {
    stop(sprintf(paste0("x must be a binary series of 0s and 1s; found %d distinct values ",
                        "in [%g, %g]. Threshold or symbolise the series first, for example ",
                        "x = as.numeric(x > median(x))."), length(u), min(u), max(u)),
         call. = FALSE)
  }
  nW <- length(x) - dim + 1
  words <- numeric(nW)
  for (q in seq_len(dim)) words <- words * 2 + x[q:(q + nW - 1)]
  Sm <- 2^dim
  P <- tabulate(words + 1, Sm) / nW
  H <- P * log2(P)
  H <- -sum(H[!is.nan(H)])
  So <- length(unique(words))
  CSE <- H + (So - 1) / (2 * Sm * log(2))
  CSEm <- -log2(1 / Sm) + (Sm - 1) / (2 * Sm * log(2))
  CSE / CSEm
}

#' Weighted recurrence entropy
#'
#' Shannon entropy of the column sums of a weighted recurrence plot, binned
#' into 50 equal-width bins and weighted by strength.
#'
#' @param wrp Square weighted recurrence matrix.
#' @return Weighted entropy in nats. 0 when every column sum is equal.
#' @export
ent_weighted <- function(wrp) {
  si <- colSums(as.matrix(wrp))
  lo <- min(si)
  hi <- max(si)
  if (hi == lo) return(0)
  width <- (hi - lo) / 49
  S <- sum(si)
  # Bins 1-49 by index; bin 50 holds the maximum alone, as in MATLAB.
  bin <- pmin(floor((si - lo) / width), 48) + 1
  bin[si == hi] <- 50
  p1 <- vapply(1:50, function(k) sum(si[bin == k]) / S, numeric(1))
  pp <- p1 * log(p1)
  pp[is.nan(pp)] <- 0
  -sum(pp)
}

#' Multiscale entropies: refined composite, composite, standard, fuzzy and
#' generalised
#'
#' At each scale `1:scale` the series is coarse-grained by non-overlapping
#' means (and, for GMSE, variances), and sample entropy is computed in
#' several ways.
#'
#' Two details follow the MATLAB implementation exactly and may be worth
#' revisiting there: GMSE uses `radius` as an absolute tolerance on the
#' coarse-grained variances, not scaled by any SD; and the fuzzy entropy
#' scales the already-scaled tolerance `radius * sd(x)` again by the SD of
#' the coarse-grained series.
#'
#' @param x Numeric vector.
#' @param scale Largest scale.
#' @param dim Template length.
#' @param radius Tolerance, as a proportion of `sd(x)`.
#' @return A list of five vectors of length `scale`: `RCMSE`, `CMSE`, `MSE`,
#'   `MSFE` and `GMSE`.
#' @references Wu, S.-D., Wu, C.-W., Lin, S.-G., Lee, K.-Y. and Peng, C.-K.
#'   (2014). Analysis of complex time series using refined composite
#'   multiscale entropy. Physics Letters A, 378(20), 1369-1374.
#' @export
ent_ms_plus <- function(x, scale, dim, radius) {
  x <- as.numeric(x)
  R <- radius * stats::sd(x)
  N <- length(x)
  RCMSE <- CMSE <- MSE <- MSFE <- GMSE <- numeric(scale)
  for (i in seq_len(scale)) {
    nj <- floor(N / i)
    # Row k, column j: the window starting at (j-1)*i + k, length i. Windows
    # running past the end are NaN, as MATLAB's try/catch makes them.
    win_stat <- function(f) {
      out <- matrix(NaN, i, nj)
      for (j in seq_len(nj)) for (k in seq_len(i)) {
        hiIdx <- j * i + k - 1
        if (hiIdx <= N) out[k, j] <- f(x[((j - 1) * i + k):hiIdx])
      }
      out
    }
    o2 <- win_stat(function(w) if (length(w) > 1) stats::var(w) else 0)
    GMSE[i] <- ms_samp_ent(o2[1, ], dim, radius)$SE
    y <- win_stat(function(w) sum(w) / i)
    row1 <- y[1, !is.nan(y[1, ])]
    MSE[i] <- ms_samp_ent(row1, dim, R)$SE
    MSFE[i] <- ms_fuzzy_ent(row1, dim, R, 2)
    nm <- nm1 <- numeric(i)
    for (k in seq_len(i)) {
      s <- ms_samp_ent(y[k, !is.nan(y[k, ])], dim, R)
      nm[k] <- s$sum_nm
      nm1[k] <- s$sum_nm1
      CMSE[i] <- CMSE[i] + 1 / i * -log(nm1[k] / nm[k])
    }
    RCMSE[i] <- -log((1 / i * sum(nm1)) / (1 / i * sum(nm)))
  }
  list(RCMSE = RCMSE, CMSE = CMSE, MSE = MSE, MSFE = MSFE, GMSE = GMSE)
}

# Sample entropy as ent_ms_plus computes it: N - m templates, self-matches
# removed, and the A and B match totals returned for the composite forms.
ms_samp_ent <- function(data, m, r) {
  N <- length(data)
  n <- N - m
  Bm <- Am <- numeric(n)
  sum_nm <- sum_nm1 <- 0
  for (i in seq_len(n)) {
    dj <- abs(data[seq_len(n)] - data[i])
    for (k in seq_len(m - 1)) dj <- pmax(dj, abs(data[seq_len(n) + k] - data[i + k]))
    dj1 <- pmax(dj, abs(data[seq_len(n) + m] - data[i + m]))
    nmi <- sum(dj <= r) - 1
    nm1i <- sum(dj1 <= r) - 1
    sum_nm <- sum_nm + nmi
    sum_nm1 <- sum_nm1 + nm1i
    Bm[i] <- nmi / n
    Am[i] <- nm1i / n
  }
  list(SE = -log((sum(Am) / n) / (sum(Bm) / n)), sum_nm = sum_nm, sum_nm1 = sum_nm1)
}

# Fuzzy entropy (Chen et al. 2007), as ent_ms_plus computes it.
ms_fuzzy_ent <- function(series, dim, r, n) {
  N <- length(series)
  phi <- numeric(2)
  r <- r * stats::sd(series)
  for (j in 1:2) {
    m <- dim + j - 1
    nP <- N - m + 1
    P <- matrix(0, m, nP)
    for (q in seq_len(m)) P[q, ] <- series[q:(N - m + q)]
    P <- sweep(P, 2, colMeans(P))
    aux <- numeric(nP)
    for (i in seq_len(N - m)) {
      d <- apply(abs(P - P[, i]), 2, max)
      simi <- exp(-(d^n) / r)
      aux[i] <- (sum(simi) - 1) / (N - m - 1)
    }
    phi[j] <- sum(aux) / (N - m)
  }
  log(phi[1]) - log(phi[2])
}

#' Entropy of a time series, by algorithm
#'
#' One entry point for the entropy family.
#'
#' @param x Numeric vector, no `NA`.
#' @param algorithm `"sample"` (default), `"approximate"`, `"permutation"`,
#'   `"symbolic"` (binary series only) or `"multiscale"` (refined composite
#'   multiscale sample entropy).
#' @param y Optional second series: gives the cross entropy. Only with
#'   `"sample"` or `"approximate"`.
#' @param dim Template length or pattern order. Default 2.
#' @param radius Tolerance, as a proportion of the SD. Default 0.2.
#' @param delay Delay for `"permutation"`. Default 1.
#' @param scale Largest scale for `"multiscale"`. Default 20.
#' @return A list with `value` and `extra`, a list of diagnostics.
#' @export
ent <- function(x, algorithm = "sample", y = NULL, dim = 2, radius = 0.2, delay = 1,
                scale = 20) {
  x <- as.numeric(x)
  if (anyNA(x)) stop("x contains NA. Entropy is undefined for missing data; remove or impute first.",
                     call. = FALSE)
  if (length(x) <= dim + 1) {
    stop(sprintf("Series of %d samples is too short for dim = %d.", length(x), as.integer(dim)),
         call. = FALSE)
  }
  isCross <- !is.null(y)
  if (isCross) {
    y <- as.numeric(y)
    if (anyNA(y)) stop("y contains NA.", call. = FALSE)
    if (length(y) != length(x)) {
      stop(sprintf("x has %d samples and y has %d; a cross entropy needs equal lengths.",
                   length(x), length(y)), call. = FALSE)
    }
  }
  algorithm <- tolower(algorithm)
  rejectCross <- function(name) {
    if (isCross) {
      stop(sprintf(paste0('algorithm "%s" has no cross-entropy form. Supply y only with ',
                          '"sample" or "approximate".'), name), call. = FALSE)
    }
  }
  extra <- list(algorithm = algorithm, dim = dim, cross = isCross)
  value <- switch(algorithm,
    sample = {
      extra$radius <- radius
      extra$estimator <- if (isCross) "cross sample" else "sample"
      if (isCross) ent_xsamp(x, y, dim, radius, 1) else ent_samp(x, dim, radius)
    },
    approximate = {
      extra$radius <- radius
      extra$estimator <- if (isCross) "cross approximate" else "approximate"
      if (isCross) ent_xap(x, y, dim, radius, 1) else ent_ap(x, dim, radius)
    },
    permutation = {
      rejectCross("permutation")
      p <- ent_permu(x, dim, delay)
      extra$delay <- delay
      extra$patternCounts <- p$hist
      extra$estimator <- "permutation"
      p$permEnt
    },
    symbolic = {
      rejectCross("symbolic")
      u <- unique(x)
      if (!all(u %in% c(0, 1))) {
        stop(sprintf(paste0('algorithm "symbolic" needs a binary series; x has %d distinct ',
                            "values in [%g, %g]. Threshold it first, for example ",
                            "x = as.numeric(x > median(x))."), length(u), min(u), max(u)),
             call. = FALSE)
      }
      extra$estimator <- "symbolic"
      ent_symbolic(x, dim)
    },
    multiscale = {
      rejectCross("multiscale")
      ms <- ent_ms_plus(x, scale, dim, radius)
      extra$radius <- radius
      extra$scale <- scale
      extra$rcmse <- ms$RCMSE
      extra$cmse <- ms$CMSE
      extra$mse <- ms$MSE
      extra$msfe <- ms$MSFE
      extra$gmse <- ms$GMSE
      extra$estimator <- "refined composite multiscale sample"
      ms$RCMSE
    },
    stop(sprintf(paste0('Unknown algorithm "%s". Supported: "sample", "approximate", ',
                        '"permutation", "symbolic", "multiscale".'), algorithm), call. = FALSE))
  list(value = value, extra = extra)
}
