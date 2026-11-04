#' Sample entropy
#'
#' Richman and Moorman (2000). Both the matches of length `dim` (B) and of
#' length `dim + 1` (A) are counted over the same `N - dim` templates, and
#' self-matches are excluded.
#'
#' @param x Numeric vector.
#' @param dim Template length (`m` in the paper).
#' @param radius Tolerance (`r` in the paper). Default 0.2.
#' @param flag `"prop"` (default) scales `radius` by `sd(x)`; anything else
#'   uses it as an absolute distance.
#' @return Sample entropy in nats. `Inf` if no template of length `dim`
#'   matches.
#' @references Richman, J. S. and Moorman, J. R. (2000). Physiological
#'   time-series analysis using approximate entropy and sample entropy.
#'   American Journal of Physiology, 278(6), H2039-H2049.
#' @export
ent_samp <- function(x, dim, radius = 0.2, flag = "prop") {
  x <- as.numeric(x)
  N <- length(x)
  if (N < dim + 1) stop("The time series must have at least m+1 data points.", call. = FALSE)
  if (tolower(flag) == "prop") radius <- radius * stats::sd(x)

  .ent_samp_core(x, dim, radius)
}

# The R code for ent_samp, kept as the readable reference for the compiled
# kernel (cpp/src/ent_samp.cpp); the tests check that the two agree.
ent_samp_r <- function(x, dim, radius) {
  N <- length(x)
  n <- N - dim
  # Template j of length L starts at sample j; column k holds x[j + k - 1].
  B <- count_matches(x, n, dim, radius)
  A <- count_matches(x, n, dim + 1, radius)
  if (B == 0) Inf else -log(A / B)
}

# Number of pairs i < j among the first n templates of length L whose
# Chebyshev distance is within radius.
count_matches <- function(x, n, L, radius) {
  total <- 0
  for (i in seq_len(n - 1)) {
    j <- (i + 1):n
    d <- abs(x[j] - x[i])
    for (k in seq_len(L - 1)) {
      d <- pmax(d, abs(x[j + k] - x[i + k]))
    }
    total <- total + sum(d <= radius)
  }
  total
}
