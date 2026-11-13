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
