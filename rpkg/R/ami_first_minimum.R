# Lag of the first strict local minimum, falling back to the first lag at
# which AMI drops to `fraction` of AMI(0) (Abarbanel et al. 1993).
ami_first_minimum <- function(curve, fraction, caller) {
  v <- curve[, 2]
  k <- length(v)
  isMin <- c(FALSE, v[2:(k - 1)] < v[1:(k - 2)] & v[2:(k - 1)] <= v[3:k], FALSE)
  if (k < 3) isMin <- rep(FALSE, k)
  idx <- which(isMin)
  info <- list(allMinima = curve[idx, , drop = FALSE], usedFallback = FALSE,
               fractionLag = NA_real_)
  below <- which(v <= fraction * v[1])
  if (length(below)) info$fractionLag <- curve[below[1], 1]
  if (length(idx)) {
    tau <- curve[idx[1], 1]
  } else {
    tau <- info$fractionLag
    info$usedFallback <- TRUE
    if (is.na(tau)) {
      warning(sprintf(paste0("%s: no local minimum in AMI over lags 0:%d and the curve ",
                             "never fell to %.0f%% of its lag-0 value."),
                      caller, as.integer(curve[k, 1]), 100 * fraction), call. = FALSE)
    }
  }
  list(tau = unname(tau), curve = curve, info = info)
}
