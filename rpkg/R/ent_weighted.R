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
