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
