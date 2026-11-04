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
