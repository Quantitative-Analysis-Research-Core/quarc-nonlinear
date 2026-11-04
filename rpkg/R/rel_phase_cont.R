#' Continuous relative phase
#'
#' Phase angle between two oscillators, from position and velocity and from
#' the analytic signal (Hilbert transform). Each series is centred on the
#' midpoint of its range first.
#'
#' @param x,y Numeric vectors of equal length.
#' @param fs Sampling frequency, Hz.
#' @return A list with `crp`, from position and velocity (one shorter than
#'   the input), and `crpH`, from the Hilbert transform. Both in degrees.
#' @export
rel_phase_cont <- function(x, y, fs) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  x <- x - min(x) - (max(x) - min(x)) / 2
  y <- y - min(y) - (max(y) - min(y)) / 2
  n <- length(x)
  x0 <- x[-n]
  y0 <- y[-n]
  dx <- diff(x) * fs
  dy <- diff(y) * fs
  crp <- atan((dx * y0 - dy * x0) / (x0 * y0 - dy * dx)) * 180 / pi
  hx <- Im(hilbert(x))
  hy <- Im(hilbert(y))
  crpH <- atan((hx * y - hy * x) / (x * y + hx * hy)) * 180 / pi
  list(crp = crp, crpH = crpH)
}
