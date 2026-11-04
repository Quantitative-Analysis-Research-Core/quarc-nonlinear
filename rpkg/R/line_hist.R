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
