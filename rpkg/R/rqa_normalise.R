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
