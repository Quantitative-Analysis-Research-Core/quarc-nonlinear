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
