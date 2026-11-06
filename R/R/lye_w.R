#' Largest Lyapunov exponent, Wolf's method
#'
#' Wolf, Swift, Swinney and Vastano (1985): follow a reference trajectory
#' and a neighbour for `evolve` samples, record the log stretch, then replace
#' the neighbour with one near the evolved point and close in direction.
#'
#' @param x Numeric vector, or a phase space.
#' @param fs Sampling frequency, Hz.
#' @param delay,dim Embedding delay and dimension, used when `x` is a vector.
#' @param evolve Propagation length in samples. No default: use 1 for a map.
#' @param SCALEMX,ANGLMX,ZMULT Maximum replacement distance (default a tenth
#'   of the data range), maximum replacement angle in radians (default 30
#'   degrees) and Wolf's scale multiplier. Give all three or none.
#' @return A list with `out`, one row per replacement step (iteration,
#'   point, neighbour, start and end distance, running exponent, current
#'   scale, best angle and maximum angle in degrees), and `LyE`, the exponent
#'   in bits per unit time.
#' @references Wolf, A., Swift, J. B., Swinney, H. L. and Vastano, J. A.
#'   (1985). Determining Lyapunov exponents from a time series. Physica D,
#'   16(3), 285-317.
#' @export
lye_w <- function(x, fs, delay, dim, evolve, SCALEMX = NULL, ANGLMX = NULL, ZMULT = NULL) {
  given <- !c(is.null(SCALEMX), is.null(ANGLMX), is.null(ZMULT))
  if (!any(given)) {
    SCALEMX <- (max(x) - min(x)) / 10
    ANGLMX <- 30 * pi / 180
    ZMULT <- 1
  } else if (!all(given)) {
    stop("not enough input arguements", call. = FALSE)
  }
  DT <- 1 / fs
  x <- as.matrix(x)
  if (ncol(x) == 1) {
    Y <- psr(x, delay, dim)
    NPT <- nrow(x) - (dim - 1) * delay - evolve
    Y <- Y[1:(NPT + evolve), , drop = FALSE]
  } else {
    Y <- x
    NPT <- nrow(Y) - evolve
  }
  out <- matrix(0, floor(NPT / evolve), 9)
  ITS <- 0
  distSUM <- 0
  thbest <- 0
  OUTMX <- SCALEMX
  pair <- NA_integer_
  dist_to <- function(p) sqrt(rowSums((Y[1:NPT, , drop = FALSE] -
                                       matrix(Y[p, ], NPT, ncol(Y), byrow = TRUE))^2))
  minIndex <- function(v) { w <- which.min(v); if (length(w)) w else 1L }
  LyE <- NA_real_
  for (i in seq(1, NPT, by = evolve)) {
    cp <- i
    if (cp == 1) {
      d <- dist_to(cp)
      ex <- (cp - 10):(cp + 10)
      d[d <= 0] <- NA
      d[ex[ex >= 1 & ex <= NPT]] <- NA
      pair <- minIndex(d)
    }
    start_dist <- sqrt(sum((Y[cp, ] - Y[pair, ])^2))
    end_dist <- sqrt(sum((Y[cp + evolve, ] - Y[pair + evolve, ])^2))
    distSUM <- distSUM + log2(end_dist / start_dist) / (evolve * DT)
    ITS <- ITS + 1
    LyE <- distSUM / ITS
    row <- floor(i / evolve) + 1
    if (row > nrow(out)) out <- rbind(out, matrix(0, row - nrow(out), 9))
    out[row, ] <- c(ITS, cp, pair, start_dist, end_dist, LyE, OUTMX, thbest * 180 / pi,
                    ANGLMX * 180 / pi)
    ZMULT <- 1
    if (end_dist < SCALEMX) {
      pair <- pair + evolve
      if (pair > NPT) {
        pair <- pair - evolve
        nx <- wolf_next_point(1, Y, cp, pair, NPT, evolve, SCALEMX, ANGLMX)
        pair <- nx$next_point; ANGLMX <- nx$ANGLMX; thbest <- nx$thbest; OUTMX <- nx$SCALEMX
      }
      next
    }
    nx <- wolf_next_point(0, Y, cp, pair, NPT, evolve, SCALEMX, ANGLMX)
    pair <- nx$next_point; ANGLMX <- nx$ANGLMX; thbest <- nx$thbest; OUTMX <- nx$SCALEMX
  }
  list(out = out, LyE = LyE)
}

# Replacement neighbour near the evolved point, preferring a small angle to
# the current separation (GetNextPoint and find_next_point in MATLAB).
wolf_next_point <- function(flag, Y, cp, pair, NPT, evolve, SCALEMX, ANGLMX) {
  p1 <- Y[cp + evolve, ]
  Vnew <- matrix(p1, NPT, ncol(Y), byrow = TRUE) - Y[1:NPT, , drop = FALSE]
  d <- sqrt(rowSums(Vnew^2))
  ex <- (cp + evolve - 10):(cp + evolve + 10)
  d[ex[ex >= 1 & ex <= NPT]] <- NA
  Vcurr <- p1 - Y[pair + evolve, ]
  end_dist <- sqrt(sum(Vcurr^2))
  cosTheta <- abs(as.numeric(Vnew %*% Vcurr) / (d * end_dist))
  # MATLAB's acos of a cosine rounded just above 1 is complex: comparisons
  # use its real part (0) and abs() its modulus (acosh).
  above <- !is.na(cosTheta) & cosTheta > 1
  thetaRe <- ifelse(above, 0, acos(pmin(cosTheta, 1)))
  thetaAbs <- ifelse(above, acosh(pmax(cosTheta, 1)), thetaRe)
  minIndex <- function(v) { w <- which.min(v); if (length(w)) w else 1L }
  if (flag == 0) {
    pot <- d
    pot[!is.na(d) & (d <= 0 | thetaRe >= ANGLMX)] <- NA
    np <- minIndex(pot)
    if (!is.na(pot[np]) && pot[np] <= SCALEMX) {
      return(list(next_point = np, ANGLMX = 30 * pi / 180, thbest = abs(thetaAbs[np]),
                  SCALEMX = SCALEMX))
    }
  }
  pot <- d
  pot[!is.na(d) & d <= 0] <- NA
  list(next_point = minIndex(pot), ANGLMX = ANGLMX, thbest = ANGLMX, SCALEMX = SCALEMX)
}
