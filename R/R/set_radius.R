#' Radius giving a target percent recurrence
#'
#' Brackets the target from both sides, then bisects `iter` times.
#'
#' @inheritParams line_hist
#' @param radius_start,radius_end Initial lower and upper radii.
#' @param threshold Target percent recurrence.
#' @param iter Bisection steps.
#' @return A list with `rec`, `diag_hist`, `vertical_hist`, `rad_final` and
#'   `A`, from the last bisection step.
#' @export
set_radius <- function(x, a, radius_start, radius_end, threshold, type, iter, theiler = 0) {
  maxAdjust <- 200
  rec <- line_hist(x, a, radius_start, type, theiler)$recurrence
  k <- 0
  while (rec == 0 || rec > threshold) {
    k <- k + 1
    if (k > maxAdjust) {
      stop(sprintf(paste0("Could not bracket %%REC = %g from below after %d adjustments ",
                          "(radius %g gives %%REC = %g). The distance matrix takes too few ",
                          "distinct values for this target to be reachable, which is usual ",
                          'for binary or categorical data. Set the radius directly with ',
                          'param = "rad", or pick a target that lies on an achievable step.'),
                   threshold, maxAdjust, radius_start, rec), call. = FALSE)
    }
    radius_start <- if (rec == 0) radius_start * 2 else radius_start / 1.5
    rec <- line_hist(x, a, radius_start, type, theiler)$recurrence
  }
  rec <- line_hist(x, a, radius_end, type, theiler)$recurrence
  k <- 0
  while (rec < threshold) {
    k <- k + 1
    if (k > maxAdjust) {
      stop(sprintf("Could not bracket %%REC = %g from above after %d adjustments (radius %g gives %%REC = %g).",
                   threshold, maxAdjust, radius_end, rec), call. = FALSE)
    }
    radius_end <- radius_end * 2
    rec <- line_hist(x, a, radius_end, type, theiler)$recurrence
  }
  lv <- radius_start
  hv <- radius_end
  for (i in seq_len(iter)) {
    mid <- (lv + hv) / 2
    lh <- line_hist(x, a, mid, type, theiler)
    if (lh$recurrence < threshold) lv <- mid else hv <- mid
  }
  message(sprintf("%% recurrence = %s, radius = %s", format(lh$recurrence, digits = 5),
                  format(mid, digits = 5)))
  list(rec = lh$recurrence, diag_hist = lh$diag_hist, vertical_hist = lh$vertical_hist,
       rad_final = mid, A = lh$a)
}
