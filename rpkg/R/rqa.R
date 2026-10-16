# Recurrence quantification analysis: rqa (auto), crqa (cross), jrqa
# (joint) and mdrqa (multidimensional), with line_hist and set_radius.

#' Line histograms and percent recurrence of a thresholded distance matrix
#'
#' @param x Unused; kept for parity with MATLAB's signature.
#' @param a Negated distance matrix, or a list of them (joint recurrence:
#'   the thresholded matrices are multiplied together).
#' @param radius Recurrence radius.
#' @param type `"rqa"`, `"crqa"`, `"jrqa"` or `"mdrqa"`. Cross recurrence
#'   counts every diagonal and ignores `theiler`.
#' @param theiler Diagonals within `theiler` of the main diagonal are
#'   excluded from the diagonal histogram and from percent recurrence; the
#'   main diagonal (line of identity) always is, except in cross recurrence.
#' @return A list with `recurrence` (percent), `diag_hist` and
#'   `vertical_hist` (line lengths) and `a`, the 0/1 recurrence matrix.
#' @export
line_hist <- function(x, a, radius, type, theiler = 0) {
  isCross <- grepl("crqa", type, fixed = TRUE)
  if (isCross) theiler <- 0
  if (!is.list(a)) a <- list(a)
  a <- lapply(a, function(m) (m + radius >= 0) * 1)
  if (length(a) > 1) {
    s <- min(vapply(a, ncol, integer(1)))
    A <- a[[1]][1:s, 1:s]
    for (i in 2:length(a)) A <- A * a[[i]][1:s, 1:s]
  } else {
    A <- a[[1]]
  }
  n <- nrow(A)
  runs <- function(v) {
    r <- rle(v)
    r$lengths[r$values == 1]
  }
  diag_hist <- integer(0)
  for (k in -(n - 1):(n - 1)) {
    if (!isCross && abs(k) <= theiler) next
    diag_hist <- c(diag_hist, runs(diag_k(A, k)))
  }
  if (!isCross && length(diag_hist) == 0) diag_hist <- 0
  vertical_hist <- integer(0)
  for (j in seq_len(ncol(A))) {
    v <- A[, j]
    if (theiler > 0) v[max(1, j - theiler):min(n, j + theiler)] <- 0
    vertical_hist <- c(vertical_hist, runs(v))
  }
  if (!isCross) {
    band <- sum(diag_k(A, 0))
    for (k in seq_len(theiler)) band <- band + sum(diag_k(A, k)) + sum(diag_k(A, -k))
    excluded <- n + 2 * sum(n - seq_len(theiler))
    recurrence <- 100 * (sum(A) - band) / (n^2 - excluded)
  } else {
    recurrence <- 100 * sum(A) / n^2
  }
  list(recurrence = recurrence, diag_hist = diag_hist, vertical_hist = vertical_hist, a = A)
}

# MATLAB's diag(A, k): k > 0 above the main diagonal.
diag_k <- function(A, k) {
  n <- nrow(A)
  if (k >= 0) {
    i <- seq_len(n - k)
    A[cbind(i, i + k)]
  } else {
    i <- seq_len(n + k)
    A[cbind(i - k, i)]
  }
}

#' Radius giving a target percent recurrence
#'
#' Brackets the target from both sides, then bisects `iter` times.
#'
#' @inheritParams line_hist
#' @param radius_start,radius_end Initial lower and upper radii.
#' @param threshold Target percent recurrence.
#' @param iter Bisection steps.
#' @return A list with `rec`, `diag_hist`, `vertical_hist`, `rad_final` and
#'   `A`, from the last bisection step.
#' @export
set_radius <- function(x, a, radius_start, radius_end, threshold, type, iter, theiler = 0) {
  maxAdjust <- 200
  rec <- line_hist(x, a, radius_start, type, theiler)$recurrence
  k <- 0
  while (rec == 0 || rec > threshold) {
    k <- k + 1
    if (k > maxAdjust) {
      stop(sprintf(paste0("Could not bracket %%REC = %g from below after %d adjustments ",
                          "(radius %g gives %%REC = %g). The distance matrix takes too few ",
                          "distinct values for this target to be reachable, which is usual ",
                          'for binary or categorical data. Set the radius directly with ',
                          'param = "rad", or pick a target that lies on an achievable step.'),
                   threshold, maxAdjust, radius_start, rec), call. = FALSE)
    }
    radius_start <- if (rec == 0) radius_start * 2 else radius_start / 1.5
    rec <- line_hist(x, a, radius_start, type, theiler)$recurrence
  }
  rec <- line_hist(x, a, radius_end, type, theiler)$recurrence
  k <- 0
  while (rec < threshold) {
    k <- k + 1
    if (k > maxAdjust) {
      stop(sprintf("Could not bracket %%REC = %g from above after %d adjustments (radius %g gives %%REC = %g).",
                   threshold, maxAdjust, radius_end, rec), call. = FALSE)
    }
    radius_end <- radius_end * 2
    rec <- line_hist(x, a, radius_end, type, theiler)$recurrence
  }
  lv <- radius_start
  hv <- radius_end
  for (i in seq_len(iter)) {
    mid <- (lv + hv) / 2
    lh <- line_hist(x, a, mid, type, theiler)
    if (lh$recurrence < threshold) lv <- mid else hv <- mid
  }
  message(sprintf("%% recurrence = %s, radius = %s", format(lh$recurrence, digits = 5),
                  format(mid, digits = 5)))
  list(rec = lh$recurrence, diag_hist = lh$diag_hist, vertical_hist = lh$vertical_hist,
       rad_final = mid, A = lh$a)
}

#' Recurrence quantification analysis
#'
#' Auto-recurrence of one series (`rqa`), cross recurrence of two (`crqa`),
#' joint recurrence of several (`jrqa`) and multidimensional recurrence of
#' several channels in one phase space (`mdrqa`).
#'
#' @param x For `rqa`, a numeric vector (or a phase space with
#'   `phasespace = TRUE`). For `crqa`, a two-column matrix. For `jrqa`, a
#'   matrix with one channel per column, or with `phasespace = TRUE` a list
#'   of phase spaces. For `mdrqa`, a matrix with one channel per column.
#' @param delay,dim Embedding delay and dimension. For `jrqa`, one per
#'   channel. A channel with `dim = 1` is used unembedded.
#' @param param `"rec"` sets the radius to reach `threshold` percent
#'   recurrence; `"rad"` uses `threshold` as the radius.
#' @param threshold Target percent recurrence, or the radius.
#' @param zscore 1 (default) standardises the data first.
#' @param norm Distance normalisation: `"none"` (default), `"euc"`, `"min"`
#'   or `"max"`.
#' @param dmin,vmin Shortest diagonal and vertical line counted. Default 2.
#' @param iter Bisection steps for `param = "rec"`. Default 20.
#' @param phasespace `TRUE` if `x` is already a phase space.
#' @param theiler Theiler window, `rqa` only. Default 0.
#' @return A list with `RP`, the recurrence plot as displayed (rotated as
#'   MATLAB's `imrotate(1 - A, 90)`, so 0 marks a recurrence), and
#'   `RESULTS`: `DIM`, `EMB`, `DEL`, `RADIUS`, `NORM`, `ZSCORE`, `Size`, `REC`,
#'   `DET`, `MeanL`, `MaxL`, `EntrL`, `LAM`, `MeanV`, `MaxV`, `EntrV` and
#'   (except `jrqa`) `EntrW`.
#' @references Marwan, N., Romano, M. C., Thiel, M. and Kurths, J. (2007).
#'   Recurrence plots for the analysis of complex systems. Physics Reports,
#'   438(5-6), 237-329.
#' @export
rqa <- function(x, delay = 1, dim = 1, param = "rec", threshold = 2.5, zscore = 1,
                norm = "none", dmin = 2, vmin = 2, iter = 20, phasespace = FALSE,
                theiler = 0) {
  x <- as.matrix(x)
  if (!phasespace && ncol(x) != 1) stop("Data must be single column vector.", call. = FALSE)
  if (zscore) x <- zscore_cols(x)
  if (phasespace) {
    if (norm == "euc") {
      warning(paste0("PhaseSpace is true, so DATA was not rebuilt from TAU and DIM and DIM ",
                     'cannot be confirmed to match its width. Norm="euc" scales by DIM ',
                     'regardless; pass DIM equal to ncol(data) or use Norm="none", "min" or "max".'),
              call. = FALSE)
    }
  } else if (dim > 1) {
    x <- psr(x, delay, dim)
  }
  a <- rqa_normalise(-pairwise_dist(x, x), norm, 1 * dim)
  rqa_finish(x, a, a, param, threshold, "rqa", iter, theiler, dmin, vmin,
             list(DIM = 1, EMB = dim, DEL = delay, NORM = norm, ZSCORE = zscore), TRUE)
}

#' @rdname rqa
#' @export
crqa <- function(x, delay = 1, dim = 1, param = "rec", threshold = 2.5, zscore = 1,
                 norm = "none", dmin = 2, vmin = 2, iter = 20, phasespace = FALSE) {
  x <- as.matrix(x)
  if (!phasespace && ncol(x) != 2) stop("Data must be two column vectors.", call. = FALSE)
  if (zscore) x <- zscore_cols(x)
  if (phasespace) {
    if (norm == "euc") {
      warning(paste0("PhaseSpace is true, so DATA was not rebuilt from TAU and DIM and DIM ",
                     'cannot be confirmed to match its width. Norm="euc" scales by DIM ',
                     'regardless; pass DIM equal to ncol(data)/2 or use Norm="none", "min" or "max".'),
              call. = FALSE)
    }
  } else if (dim > 1) {
    x <- psr(x, delay, dim)
  }
  odd <- seq(1, ncol(x), by = 2)
  a <- rqa_normalise(-pairwise_dist(x[, odd, drop = FALSE], x[, odd + 1, drop = FALSE]),
                     norm, 2 * dim)
  rqa_finish(x, a, a, param, threshold, "crqa", iter, 0, dmin, vmin,
             list(DIM = 1, EMB = dim, DEL = delay, NORM = norm, ZSCORE = zscore), TRUE)
}

#' @rdname rqa
#' @export
jrqa <- function(x, delay = c(1, 1), dim = c(1, 1), param = "rec", threshold = 2.5,
                 zscore = 1, norm = "none", dmin = 2, vmin = 2, iter = 20,
                 phasespace = FALSE) {
  if (phasespace) {
    if (!is.list(x)) {
      stop(paste0("PhaseSpace is true, so DATA must be a list with one already-reconstructed ",
                  "phase space per channel, e.g. list(Y1, Y2), rather than a raw multi-column series."),
           call. = FALSE)
    }
    if (length(x) < 2) stop("DATA must contain at least two channels.", call. = FALSE)
  } else {
    if (is.list(x) && !is.data.frame(x)) {
      stop("DATA must be a numeric matrix unless PhaseSpace is true.", call. = FALSE)
    }
    x <- as.matrix(x)
    if (ncol(x) < 2) stop("Data must be at least two column vectors.", call. = FALSE)
  }
  if (zscore) x <- if (phasespace) lapply(x, zscore_cols) else zscore_cols(x)
  if (phasespace) {
    data2 <- lapply(x, as.matrix)
    if (norm == "euc") {
      warning(paste0("PhaseSpace is true, so DATA was not rebuilt from TAU and DIM and DIM ",
                     'cannot be confirmed to match its width. Norm="euc" scales by DIM ',
                     "regardless; pass DIM to match each phase space's width or use ",
                     'Norm="none", "min" or "max".'), call. = FALSE)
    }
  } else {
    data2 <- lapply(seq_len(ncol(x)), function(i) {
      if (dim[i] > 1) psr(x[, i], delay[i], dim[i]) else x[, i, drop = FALSE]
    })
  }
  DIM <- length(data2)
  if (norm == "euc" && length(unique(dim)) > 1) {
    stop(paste0('Norm="euc" scales by DIM, which is per channel here and differs between ',
                'channels. Use Norm="none", "min" or "max".'), call. = FALSE)
  }
  a <- lapply(data2, function(d) rqa_normalise(-pairwise_dist(d, d), norm, DIM * dim[1]))
  rqa_finish(x, a, a[[1]], param, threshold, "jrqa", iter, 0, dmin, vmin,
             list(DIM = 1, EMB = dim, DEL = delay, NORM = norm, ZSCORE = zscore), FALSE)
}

#' @rdname rqa
#' @export
mdrqa <- function(x, delay = 1, dim = 1, param = "rec", threshold = 2.5, zscore = 1,
                  norm = "none", dmin = 2, vmin = 2, iter = 20, phasespace = FALSE) {
  x <- as.matrix(x)
  if (zscore) x <- zscore_cols(x)
  if (phasespace) {
    if (norm == "euc") {
      warning(paste0("PhaseSpace is true, so DATA was not rebuilt from TAU and DIM and DIM ",
                     'cannot be confirmed to match its width. Norm="euc" scales by DIM ',
                     'regardless; pass DIM equal to ncol(data) or use Norm="none", "min" or "max".'),
              call. = FALSE)
    }
  } else if (dim > 1) {
    x <- psr(x, delay, dim)
  }
  a <- rqa_normalise(-pairwise_dist(x, x), norm, 1 * dim)
  rqa_finish(x, a, a, param, threshold, "mdrqa", iter, 0, dmin, vmin,
             list(DIM = 1, EMB = dim, DEL = delay, NORM = norm, ZSCORE = zscore), TRUE)
}

# ------------------------------------------------------------------ helpers

zscore_cols <- function(x) {
  x <- as.matrix(x)
  sweep(sweep(x, 2, colMeans(x)), 2, apply(x, 2, stats::sd), "/")
}

# Euclidean distances between the rows of X and the rows of Y.
pairwise_dist <- function(X, Y) {
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  s <- matrix(0, nrow(X), nrow(Y))
  for (k in seq_len(ncol(X))) s <- s + outer(X[, k], Y[, k], `-`)^2
  sqrt(s)
}

rqa_normalise <- function(a, norm, scaleDim) {
  if (norm == "euc") {
    b <- mean(a[a < 0])
    b <- -sqrt(abs(b^2 + 2 * scaleDim))
    a / abs(b)
  } else if (norm == "min") {
    a / abs(max(a[a < 0]))
  } else if (norm == "max") {
    a / abs(min(a[a < 0]))
  } else {
    a
  }
}

rqa_finish <- function(x, a, wrp, param, threshold, type, iter, theiler, dmin, vmin,
                       RESULTS, withEntrW) {
  if (param == "rad") {
    radius <- threshold
    lh <- line_hist(x, a, threshold, type, theiler)
    recurrence <- lh$recurrence
    dh <- lh$diag_hist
    vh <- lh$vertical_hist
    A <- lh$a
  } else if (param == "rec") {
    sr <- set_radius(x, a, 0.01, 0.5, threshold, type, iter, theiler)
    recurrence <- sr$rec
    dh <- sr$diag_hist
    vh <- sr$vertical_hist
    radius <- sr$rad_final
    A <- sr$A
  } else {
    stop('param must be "rad" or "rec".', call. = FALSE)
  }
  RESULTS$RADIUS <- radius
  RESULTS$Size <- nrow(A)
  RESULTS$REC <- recurrence
  line_entropy <- function(h) {
    p <- tabulate(h)
    p <- p[p > 0] / sum(p)
    -sum(p * log2(p))
  }
  if (recurrence > 0) {
    d <- dh[dh >= dmin]
    v <- vh[vh >= vmin]
    RESULTS$DET <- 100 * sum(d) / sum(dh)
    RESULTS$MeanL <- mean(d)
    RESULTS$MaxL <- max(d)
    RESULTS$EntrL <- line_entropy(d)
    RESULTS$LAM <- 100 * sum(v) / sum(vh)
    RESULTS$MeanV <- mean(v)
    RESULTS$MaxV <- max(v)
    RESULTS$EntrV <- line_entropy(v)
    if (withEntrW) RESULTS$EntrW <- ent_weighted(wrp)
  } else {
    for (f in c("DET", "MeanL", "MaxL", "EntrL", "LAM", "MeanV", "MaxV", "EntrV")) RESULTS[[f]] <- NaN
    if (withEntrW) RESULTS$EntrW <- NaN
  }
  M <- 1 - A
  list(RP = t(M)[ncol(M):1, , drop = FALSE], RESULTS = RESULTS)
}
