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
