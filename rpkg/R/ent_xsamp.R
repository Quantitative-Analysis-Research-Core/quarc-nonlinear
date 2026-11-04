#' Cross sample entropy
#'
#' Richman and Moorman (2000): how often a template of `y` is matched by
#' templates of `x`, for lengths `dim` and `dim + 1`.
#'
#' @param x,y Numeric vectors of equal length.
#' @param dim Template length.
#' @param radius Tolerance.
#' @param norm Normalisation: `1` rescales both series to the range 0-1 and
#'   sets the tolerance to `radius` times the mean of their SDs; `2`
#'   standardises both and uses `radius` as is; anything else leaves the data
#'   alone and uses `radius` times the mean SD of the raw series.
#' @return Cross sample entropy in nats.
#' @export
ent_xsamp <- function(x, y, dim, radius, norm) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  if (length(x) != length(y)) stop("The data series need to be the same length!", call. = FALSE)
  N <- length(x)
  if (norm == 1) {
    xn <- (x - min(x)) / (max(x) - min(x))
    yn <- (y - min(y)) / (max(y) - min(y))
    r <- radius * ((stats::sd(xn) + stats::sd(yn)) / 2)
  } else if (norm == 2) {
    xn <- (x - mean(x)) / stats::sd(x)
    yn <- (y - mean(y)) / stats::sd(y)
    r <- radius
  } else {
    xn <- x
    yn <- y
    r <- radius * ((stats::sd(x) + stats::sd(y)) / 2)
  }
  n <- N - dim
  Bm <- numeric(n)
  Am <- numeric(n)
  for (i in seq_len(n)) {
    dj <- abs(xn[seq_len(n)] - yn[i])
    for (k in seq_len(dim - 1)) dj <- pmax(dj, abs(xn[seq_len(n) + k] - yn[i + k]))
    dj1 <- pmax(dj, abs(xn[seq_len(n) + dim] - yn[i + dim]))
    Bm[i] <- sum(dj <= r) / n
    Am[i] <- sum(dj1 <= r) / n
  }
  -log((sum(Am) / n) / (sum(Bm) / n))
}
