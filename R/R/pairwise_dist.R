# Euclidean distances between the rows of X and the rows of Y.
pairwise_dist <- function(X, Y) {
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  s <- matrix(0, nrow(X), nrow(Y))
  for (k in seq_len(ncol(X))) s <- s + outer(X[, k], Y[, k], `-`)^2
  sqrt(s)
}
