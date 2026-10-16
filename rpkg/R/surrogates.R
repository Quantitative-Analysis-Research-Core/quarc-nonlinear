#' Surrogate data, Theiler's algorithms 0, 1 and 2
#'
#' Theiler et al. (1992). Algorithm 0 shuffles the series (destroys all
#' temporal structure, keeps the distribution). Algorithm 1 randomises the
#' Fourier phases (keeps the amplitude spectrum exactly, so the
#' autocorrelation; not the distribution). Algorithm 2 is the amplitude
#' adjusted Fourier transform (keeps the distribution exactly and the
#' spectrum approximately).
#'
#' Draws come from R's generator, so a surrogate differs from MATLAB's even
#' at the same seed; what each algorithm preserves is the same, and that is
#' what the tests check.
#'
#' @param x Numeric vector, or a matrix with one series per column.
#' @param algorithm 0, 1 or 2.
#' @return A surrogate of the same shape as `x`.
#' @references Theiler, J., Eubank, S., Longtin, A., Galdrikian, B. and
#'   Farmer, J. D. (1992). Testing for nonlinearity in time series: the
#'   method of surrogate data. Physica D, 58(1-4), 77-94.
#' @export
surr_theiler <- function(x, algorithm) {
  isVec <- is.null(dim(x))
  x <- as.matrix(x)
  z <- switch(as.character(algorithm),
    "0" = apply(x, 2, function(v) v[order(stats::rnorm(length(v)))]),
    "1" = surr1(x, 1),
    "2" = surr1(x, 2),
    stop(sprintf("Unknown algorithm %s. Supported: 0, 1, 2.", format(algorithm)), call. = FALSE))
  z <- matrix(z, nrow(x))
  if (isVec) as.numeric(z) else z
}

surr1 <- function(x, algorithm) {
  r <- nrow(x)
  c <- ncol(x)
  if (algorithm == 2) {
    # Gaussianise: replace each value by the same-rank value of a sorted
    # normal sample.
    sr <- apply(matrix(stats::rnorm(r * c), r, c), 2, sort)
    sx <- apply(x, 2, sort)
    y <- sapply(seq_len(c), function(k) sr[rank_of(x[, k]), k])
    y <- matrix(y, r, c)
  } else {
    sx <- NULL
    y <- x
  }
  m <- colMeans(y)
  y <- sweep(y, 2, m)
  fy <- stats::mvfft(y)
  nHalf <- floor((r - 1) / 2)
  rot <- rep(1 + 0i, r)
  if (nHalf > 0) {
    ph <- 2 * pi * stats::runif(nHalf)
    rot[2:(nHalf + 1)] <- exp(1i * ph)
    rot[r:(r - nHalf + 1)] <- Conj(exp(1i * ph))
  }
  yy <- Re(stats::mvfft(fy * rot, inverse = TRUE)) / r
  yy <- sweep(yy, 2, m, "+")
  if (algorithm == 2) {
    yy <- sapply(seq_len(c), function(k) sx[rank_of(yy[, k]), k])
  }
  matrix(yy, r, c)
}

#' Pseudo-periodic surrogate
#'
#' Small, Yu and Harrison (2001). A new series is built by hopping through
#' the reconstructed attractor: from the current point, the next is chosen
#' with probability falling off as `exp(-distance / rho)`, and its successor
#' is emitted.
#'
#' @param x Numeric vector.
#' @param delay,dim Embedding delay and dimension.
#' @param rho Noise radius. See [surr_find_rho()].
#' @return A list with `ys`, the surrogate, and `yi`, the index of `x` that
#'   each sample was taken from.
#' @references Small, M., Yu, D. and Harrison, R. G. (2001). Surrogate test
#'   for pseudoperiodic time series data. Physical Review Letters, 87(18),
#'   188101.
#' @export
surr_pseudo_periodic <- function(x, delay, dim, rho) {
  x <- as.numeric(x)
  Y <- psr(x, delay, dim)
  nY <- nrow(Y)
  xi <- floor(stats::runif(1) * nY) + 1
  ys <- numeric(nY)
  yi <- numeric(nY)
  ys[1] <- x[xi]
  yi[1] <- xi
  M <- nY - 2
  for (i in 2:nY) {
    prob <- exp(-sqrt(rowSums((Y[1:M, , drop = FALSE] -
                               matrix(Y[xi, ], M, ncol(Y), byrow = TRUE))^2)) / rho)
    sum3 <- cumsum(prob)
    xi_n <- NA
    ind <- 1
    while (is.na(xi_n) || xi_n > nY || xi_n == xi) {
      hit <- which(sum3 < stats::runif(1) * sum3[M])
      xi_n <- if (length(hit)) max(hit) + 2 else NA
      ind <- ind + 1
      if (ind == 100) {
        stop("a new value of xi could not be found, check that rho is not too low", call. = FALSE)
      }
    }
    xi <- xi_n
    ys[i] <- x[xi]
    yi[i] <- xi
  }
  list(ys = ys, yi = yi)
}

#' Noise radius for pseudo-periodic surrogates
#'
#' Bisects for the `rho` that maximises the number of short segments copied
#' from the original series, the criterion of Small et al. (2001).
#'
#' @inheritParams surr_pseudo_periodic
#' @param rhoL,rhoH Search bounds. Defaults 0.1 and 1.
#' @return A list with `rho` and `out`, one row per trial: trial, rho and
#'   the segment count.
#' @export
surr_find_rho <- function(x, delay, dim, rhoL = 0.1, rhoH = 1) {
  if (length(rhoL) != 1 || length(rhoH) != 1 || !(rhoL > 0) || !(rhoH > rhoL)) {
    stop("Search bounds must satisfy 0 < rhoL < rhoH.", call. = FALSE)
  }
  di_of <- function(yi, n = 2) sum(diff(which(diff(yi) != 1)) > n)
  diH <- di_of(surr_pseudo_periodic(x, delay, dim, rhoH)$yi)
  diL <- di_of(surr_pseudo_periodic(x, delay, dim, rhoL)$yi)
  out <- rbind(c(1, rhoH, diH), c(2, rhoL, diL))
  if (diH > diL) { dmax <- diH; rho <- rhoH } else { dmax <- diL; rho <- rhoL }
  ind <- 3
  while (abs(rhoH - rhoL) / rhoL > 0.02) {
    rhoi <- (rhoH + rhoL) / 2
    di <- di_of(surr_pseudo_periodic(x, delay, dim, rhoi)$yi)
    out <- rbind(out, c(ind, rhoi, di))
    ind <- ind + 1
    if (di > dmax) { dmax <- di; rho <- rhoi }
    if (diL < diH) { diL <- di; rhoL <- rhoi } else { diH <- di; rhoH <- rhoi }
  }
  lo <- out[2, 2]
  hi <- out[1, 2]
  margin <- 0.15 * log(hi / lo)
  if (log(rho / lo) < margin || log(hi / rho) < margin) {
    warning(sprintf(paste0("The best rho found (%g) lies against a search bound [%g %g], so the ",
                           "optimum is probably not bracketed and this value may not be a ",
                           "maximum of di(rho). Re-run with wider bounds, e.g. ",
                           "surr_find_rho(y, delay, dim, %g, %g)."), rho, lo, hi, lo / 10, hi * 2),
            call. = FALSE)
  }
  colnames(out) <- c("trial", "rho", "di")
  list(rho = rho, out = out)
}

#' Fractional Gaussian noise
#'
#' Davies-Harte circulant embedding: exact fGn of length `n` with Hurst
#' exponent `H`, zero mean and unit variance.
#'
#' @param n Length.
#' @param H Hurst exponent, in (0, 1).
#' @return Numeric vector of length `n`.
#' @references Davies, R. B. and Harte, D. S. (1987). Tests for Hurst
#'   effect. Biometrika, 74(1), 95-101.
#' @export
fgn_sim <- function(n, H) {
  z <- stats::rnorm(2 * n)
  zr <- z[1:n]
  zi <- z[(n + 1):(2 * n)]
  zic <- -zi
  zi[1] <- 0
  zr[1] <- zr[1] * sqrt(2)
  zi[n] <- 0
  zr[n] <- zr[n] * sqrt(2)
  back <- if (n > 2) (n - 1):2 else integer(0)
  zr <- c(zr[1:n], zr[back])
  zi <- c(zi[1:n], zic[back])
  z <- complex(real = zr, imaginary = zi)
  k <- 0:(n - 1)
  gammak <- (abs(k - 1)^(2 * H) - 2 * abs(k)^(2 * H) + abs(k + 1)^(2 * H)) / 2
  ind <- c(0:(n - 2), n - 1, if (n > 2) (n - 2):1 else integer(0))
  # MATLAB's ifft(v)*length(z) is an unnormalised inverse transform.
  gk <- Re(stats::fft(gammak[ind + 1], inverse = TRUE)) * length(z) / length(ind)
  if (!all(gk > 0)) stop("Re(gk)-vector not positive", call. = FALSE)
  z <- z * sqrt(gk)
  z <- stats::fft(z, inverse = TRUE) / length(z) * length(z)
  z <- 0.5 * (n - 1)^(-0.5) * z
  Re(z[1:n])
}

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

#' Analytic signal by the FFT method
#'
#' The same construction as MATLAB's `hilbert`: zero the negative
#' frequencies, double the positive ones, invert.
#'
#' @param x Numeric vector.
#' @return Complex vector; `Re` is `x`, `Im` its Hilbert transform.
#' @export
hilbert <- function(x) {
  n <- length(x)
  h <- numeric(n)
  if (n %% 2 == 0) {
    h[c(1, n / 2 + 1)] <- 1
    h[2:(n / 2)] <- 2
  } else {
    h[1] <- 1
    h[2:((n + 1) / 2)] <- 2
  }
  stats::fft(stats::fft(x) * h, inverse = TRUE) / n
}
