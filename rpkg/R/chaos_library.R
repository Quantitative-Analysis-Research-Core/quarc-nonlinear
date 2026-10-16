#' Trajectories of standard chaotic systems
#'
#' Flows are integrated with [ode45()] or [ode23()], which follow MATLAB's
#' solvers of the same names and default tolerances (relative 1e-3,
#' absolute 1e-6); maps are iterated directly.
#'
#' A flow agrees with MATLAB closely over a short horizon, then drifts:
#' rounding differences grow exponentially on a chaotic attractor, which no
#' two implementations avoid. A map is plain arithmetic and agrees exactly.
#'
#' @param S System name: `"Lorenz"`, `"Rossler"`, `"Aizawa"` (ode23),
#'   `"DequanLi"` (ode23), `"NoseHoover"`, `"QiChen"`, `"YuWang"`, `"LuChen"`,
#'   `"Arneodo"`, `"TSUCSI"`, or the maps `"Hennon"` and `"Logistic"`.
#' @param t Output times for a flow. For a map, the number of iterates, or
#'   a vector whose length sets it.
#' @param IC Initial condition.
#' @param p Parameters, in the order of the MATLAB function.
#' @return A list with `t` and `y`, one row per time and one column per
#'   variable.
#' @export
chaos_library <- function(S, t, IC, p) {
  flow <- function(f, solver = ode45) {
    out <- solver(function(tt, y) f(tt, y, p), t, IC)
    list(t = out$t, y = out$y)
  }
  switch(S,
    Rossler = flow(rossler_rhs),
    Lorenz = flow(lorenz_rhs),
    Hennon = {
      n <- map_length(t)
      x <- z <- numeric(n)
      x[1] <- IC[1]
      z[1] <- IC[2]
      for (c in seq_len(n - 1)) {
        x[c + 1] <- 1 - p[1] * x[c]^2 + z[c]
        z[c + 1] <- p[2] * x[c]
      }
      list(t = map_time(t, n), y = cbind(x, z, deparse.level = 0))
    },
    Logistic = {
      n <- map_length(t)
      x <- numeric(n)
      x[1] <- IC[1]
      for (c in seq_len(n - 1)) x[c + 1] <- p[1] * x[c] * (1 - x[c])
      list(t = map_time(t, n), y = matrix(x))
    },
    Aizawa = flow(aizawa_rhs, ode23),
    DequanLi = flow(dequanli_rhs, ode23),
    NoseHoover = flow(nosehoover_rhs),
    QiChen = flow(qichen_rhs),
    YuWang = flow(yuwang_rhs),
    LuChen = flow(luchen_rhs),
    Arneodo = flow(arneodo_rhs),
    TSUCSI = flow(tsucsi_rhs),
    stop(sprintf("Unknown system '%s'.", S), call. = FALSE))
}

map_length <- function(t) {
  n <- if (length(t) == 1) round(t) else length(t)
  if (n < 1) stop("t must ask for at least one data point.", call. = FALSE)
  n
}
map_time <- function(t, n) if (length(t) == 1) 0:(n - 1) else as.numeric(t)

rossler_rhs <- function(t, y, p) c(-y[2] - y[3], y[1] + p[1] * y[2], p[2] + y[3] * (y[1] - p[3]))
lorenz_rhs <- function(t, y, p) c(p[1] * (y[2] - y[1]), -y[1] * y[3] + p[2] * y[1] - y[2],
                                  y[1] * y[2] - p[3] * y[3])
dequanli_rhs <- function(t, y, p) c(p[1] * (y[2] - y[1]) + p[2] * (y[1] * y[3]),
                                    p[3] * y[1] + p[4] * y[2] - y[1] * y[3],
                                    p[5] * y[3] + y[1] * y[2] - p[6] * (y[1] * y[1]))
nosehoover_rhs <- function(t, y, p) c(y[2], -y[1] + y[2] * y[3], p[1] - y[2] * y[2])
qichen_rhs <- function(t, y, p) c(p[1] * (y[2] - y[1]) + y[2] * y[3], p[3] * y[1] + y[2] - y[1] * y[3],
                                  y[1] * y[2] - p[2] * y[3])
yuwang_rhs <- function(t, y, p) c(p[1] * (y[2] - y[1]), p[2] * y[1] - p[3] * y[1] * y[3],
                                  exp(y[1] * y[2]) - p[4] * y[3])
luchen_rhs <- function(t, y, p) c(-1 * (p[1] * p[2] * y[1] / (p[1] + p[2])) - y[2] * y[3] + p[3],
                                  p[1] * y[2] + y[1] * y[3], p[2] * y[3] + y[1] * y[2])
arneodo_rhs <- function(t, y, p) c(y[2], y[3], -p[1] * y[1] - p[2] * y[2] - y[3] + p[3] * y[1]^3)
tsucsi_rhs <- function(t, y, p) c(p[1] * (y[2] - y[1]) + p[3] * y[1] * y[3], p[5] * y[2] - y[1] * y[3],
                                  p[2] * y[3] + y[1] * y[2] - p[4] * y[1] * y[1])
aizawa_rhs <- function(t, y, p) c((y[3] - p[2]) * y[1] - p[4] * y[2],
                                  p[4] * y[1] + (y[3] - p[2]) * y[2],
                                  p[3] + p[1] * y[3] - y[3]^3 / 3 - y[1]^2 + p[6] * y[3] * y[1]^3)

#' Explicit Runge-Kutta solvers with MATLAB's step control
#'
#' `ode45` is the Dormand-Prince 5(4) pair and `ode23` the Bogacki-Shampine
#' 3(2) pair, with the error control, initial step, step-size update and
#' continuous extension described by Shampine and Reichelt (1997) for
#' MATLAB's solvers of the same names. Solutions are reported at the times
#' in `tspan` by interpolation; with two times, at every step and three
#' points between (MATLAB's `Refine`, 4 for ode45, 1 for ode23).
#'
#' @param f Function of `(t, y)` returning the derivative vector.
#' @param tspan Output times, increasing or decreasing.
#' @param y0 Initial state.
#' @param rtol,atol Relative and absolute tolerance. Defaults 1e-3, 1e-6.
#' @return A list with `t` and `y` (one row per time).
#' @references Shampine, L. F. and Reichelt, M. W. (1997). The MATLAB ODE
#'   suite. SIAM Journal on Scientific Computing, 18(1), 1-22.
#'
#'   Dormand, J. R. and Prince, P. J. (1980). A family of embedded
#'   Runge-Kutta formulae. Journal of Computational and Applied
#'   Mathematics, 6(1), 19-26.
#'
#'   Bogacki, P. and Shampine, L. F. (1989). A 3(2) pair of Runge-Kutta
#'   formulas. Applied Mathematics Letters, 2(4), 321-325.
#' @export
ode45 <- function(f, tspan, y0, rtol = 1e-3, atol = 1e-6) {
  tab <- list(
    pow = 1 / 5, refine = 4,
    A = c(1 / 5, 3 / 10, 4 / 5, 8 / 9, 1, 1),
    B = matrix(c(1 / 5, 3 / 40, 44 / 45, 19372 / 6561, 9017 / 3168, 35 / 384,
                 0, 9 / 40, -56 / 15, -25360 / 2187, -355 / 33, 0,
                 0, 0, 32 / 9, 64448 / 6561, 46732 / 5247, 500 / 1113,
                 0, 0, 0, -212 / 729, 49 / 176, 125 / 192,
                 0, 0, 0, 0, -5103 / 18656, -2187 / 6784,
                 0, 0, 0, 0, 0, 11 / 84,
                 0, 0, 0, 0, 0, 0), 7, 6, byrow = TRUE),
    E = c(71 / 57600, 0, -71 / 16695, 71 / 1920, -17253 / 339200, 22 / 525, -1 / 40),
    BI = matrix(c(1, -183 / 64, 37 / 12, -145 / 128,
                  0, 0, 0, 0,
                  0, 1500 / 371, -1000 / 159, 1000 / 371,
                  0, -125 / 32, 125 / 12, -375 / 64,
                  0, 9477 / 3392, -729 / 106, 25515 / 6784,
                  0, -11 / 7, 11 / 3, -55 / 28,
                  0, 3 / 2, -4, 5 / 2), 7, 4, byrow = TRUE))
  rk_solve(f, tspan, y0, rtol, atol, tab)
}

#' @rdname ode45
#' @export
ode23 <- function(f, tspan, y0, rtol = 1e-3, atol = 1e-6) {
  tab <- list(
    pow = 1 / 3, refine = 1,
    A = c(1 / 2, 3 / 4, 1),
    B = matrix(c(1 / 2, 0, 2 / 9,
                 0, 3 / 4, 1 / 3,
                 0, 0, 4 / 9,
                 0, 0, 0), 4, 3, byrow = TRUE),
    E = c(-5 / 72, 1 / 12, 1 / 9, -1 / 8),
    BI = matrix(c(1, -4 / 3, 5 / 9,
                  0, 1, -2 / 3,
                  0, 4 / 3, -8 / 9,
                  0, -1, 1), 4, 3, byrow = TRUE))
  rk_solve(f, tspan, y0, rtol, atol, tab)
}

rk_solve <- function(f, tspan, y0, rtol, atol, tab) {
  tspan <- as.numeric(tspan)
  y <- as.numeric(y0)
  t0 <- tspan[1]
  tfinal <- tspan[length(tspan)]
  tdir <- sign(tfinal - t0)
  rtol <- max(rtol, 100 * .Machine$double.eps)
  threshold <- atol / rtol
  pow <- tab$pow
  nst <- length(tab$A) + 1
  refineMode <- length(tspan) == 2
  hmax <- abs(0.1 * (tfinal - t0))
  t <- t0
  f1 <- f(t, y)
  K <- matrix(0, length(y), nst)
  K[, 1] <- f1
  htspan <- abs(tspan[2] - tspan[1])
  absh <- min(hmax, htspan)
  rh <- max(abs(f1 / pmax(abs(y), threshold))) / (0.8 * rtol^pow)
  if (absh * rh > 1) absh <- 1 / rh
  eps_at <- function(v) { v <- abs(v); if (v == 0) 2^-1074 else 2^(floor(log2(v)) - 52) }
  absh <- max(absh, 16 * eps_at(t))

  tout <- t0
  yout <- list(y)
  nextIdx <- 2
  done <- FALSE
  while (!done) {
    hmin <- 16 * eps_at(t)
    absh <- min(hmax, max(hmin, absh))
    h <- tdir * absh
    if (1.1 * absh >= abs(tfinal - t)) {
      h <- tfinal - t
      absh <- abs(h)
      done <- TRUE
    }
    nofailed <- TRUE
    repeat {
      hB <- h * tab$B
      for (j in 2:(nst - 1)) {
        K[, j] <- f(t + h * tab$A[j - 1], y + K[, 1:(j - 1), drop = FALSE] %*% hB[1:(j - 1), j - 1])
      }
      tnew <- t + h * tab$A[nst - 1]
      if (done) tnew <- tfinal
      h <- tnew - t
      ynew <- as.numeric(y + K[, 1:(nst - 1), drop = FALSE] %*% hB[1:(nst - 1), nst - 1])
      K[, nst] <- f(tnew, ynew)
      err <- absh * max(abs(as.numeric(K %*% tab$E) / pmax(pmax(abs(y), abs(ynew)), threshold)))
      if (err > rtol) {
        if (absh <= hmin) stop(sprintf("Unable to meet integration tolerances at t = %g.", t),
                               call. = FALSE)
        if (nofailed) {
          nofailed <- FALSE
          absh <- max(hmin, absh * max(if (pow == 1 / 5) 0.1 else 0.5, 0.8 * (rtol / err)^pow))
        } else {
          absh <- max(hmin, 0.5 * absh)
        }
        h <- tdir * absh
        done <- FALSE
      } else {
        break
      }
    }
    interp <- function(ti) {
      s <- (ti - t) / h
      y + as.numeric(K %*% (h * tab$BI) %*% cumprod(rep(s, ncol(tab$BI))))
    }
    if (refineMode) {
      for (ti in t + h * (seq_len(tab$refine - 1) / tab$refine)) {
        tout <- c(tout, ti)
        yout[[length(yout) + 1]] <- interp(ti)
      }
      tout <- c(tout, tnew)
      yout[[length(yout) + 1]] <- ynew
    } else {
      while (nextIdx <= length(tspan) && tdir * (tnew - tspan[nextIdx]) >= 0) {
        ti <- tspan[nextIdx]
        tout <- c(tout, ti)
        yout[[length(yout) + 1]] <- if (ti == tnew) ynew else interp(ti)
        nextIdx <- nextIdx + 1
      }
    }
    if (nofailed) {
      temp <- 1.25 * (err / rtol)^pow
      absh <- if (temp > 0.2) absh / temp else 5 * absh
    }
    t <- tnew
    y <- ynew
    K[, 1] <- K[, nst]
  }
  list(t = tout, y = do.call(rbind, yout))
}
