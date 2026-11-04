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
