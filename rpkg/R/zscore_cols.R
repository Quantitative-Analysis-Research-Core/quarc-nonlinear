zscore_cols <- function(x) {
  x <- as.matrix(x)
  sweep(sweep(x, 2, colMeans(x)), 2, apply(x, 2, stats::sd), "/")
}
