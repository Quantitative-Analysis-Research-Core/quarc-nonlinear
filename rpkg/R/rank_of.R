# Ranks 1..n, ties broken by position (order() is stable, like MATLAB sort).
rank_of <- function(v) {
  r <- integer(length(v))
  r[order(v)] <- seq_along(v)
  r
}
