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
