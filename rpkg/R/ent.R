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

#' Approximate entropy
#'
#' Pincus (1991). Self-matches are counted, as in the original definition.
#'
#' @param x Numeric vector.
#' @param dim Template length (`m` in the paper).
#' @param radius Tolerance as a proportion of `sd(x)` (`r` in the paper).
#' @return Approximate entropy in nats.
#' @references Pincus, S. M. (1991). Approximate entropy as a measure of
#'   system complexity. PNAS, 88(6), 2297-2301.
#' @export
ent_ap <- function(x, dim, radius) {
  x <- as.numeric(x)
  radius <- radius * stats::sd(x)
  N <- length(x)
  phim <- numeric(2)
  for (j in 1:2) {
    m <- dim + j - 1
    n <- N - m + 1
    phi <- numeric(n)
    for (i in seq_len(n)) {
      d <- abs(x[seq_len(n)] - x[i])
      for (k in seq_len(m - 1)) {
        d <- pmax(d, abs(x[seq_len(n) + k] - x[i + k]))
      }
      phi[i] <- sum(d <= radius) / n
    }
    phim[j] <- sum(log(phi)) / n
  }
  phim[1] - phim[2]
}

#' Permutation entropy
#'
#' Bandt and Pompe (2002), in bits.
#'
#' @param x Numeric vector.
#' @param dim Order of the ordinal patterns.
#' @param delay Delay between the samples of a pattern.
#' @return A list with `permEnt`, the entropy in bits, and `hist`, the count
#'   of each ordinal pattern. `hist` is in the order of MATLAB's
#'   `perms(1:dim)` (reverse lexicographic), so it can be compared entry by
#'   entry with the MATLAB output.
#' @references Bandt, C. and Pompe, B. (2002). Permutation entropy: a
#'   natural complexity measure for time series. Physical Review Letters,
#'   88(17), 174102.
#' @export
ent_permu <- function(x, dim, delay) {
  x <- as.numeric(x)
  N <- length(x)
  perm <- matlab_perms(dim)
  nWin <- N - delay * (dim - 1)
  # The ordinal pattern of each window, as a string key; order() is stable,
  # like MATLAB's sort, so ties resolve the same way.
  keys <- character(nWin)
  for (s in seq_len(nWin)) {
    keys[s] <- paste(order(x[seq(s, s + delay * (dim - 1), by = delay)]), collapse = ",")
  }
  permKeys <- apply(perm, 1, paste, collapse = ",")
  hist <- as.numeric(table(factor(keys, levels = permKeys)))
  p <- hist[hist != 0] / sum(hist)
  list(permEnt = -sum(p * log2(p)), hist = hist)
}

# MATLAB's perms(1:n): all permutations, rows in reverse lexicographic order.
matlab_perms <- function(n) {
  if (n == 1) return(matrix(1L, 1, 1))
  sub <- matlab_perms(n - 1)
  out <- NULL
  for (first in n:1) {
    rest <- setdiff(seq_len(n), first)   # ascending: sub value k is the k-th smallest
    out <- rbind(out, cbind(first, matrix(rest[sub], nrow(sub))))
  }
  unname(out)
}
