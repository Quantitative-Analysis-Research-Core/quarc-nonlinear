#' Detrended fluctuation analysis
#'
#' Peng et al. (1995). The profile is split into non-overlapping windows of
#' each scale, a polynomial is removed from each window, and the root mean
#' square of the residuals is averaged across windows.
#'
#' @param x Numeric vector.
#' @param scales Window sizes in samples.
#' @param order Order of the detrending polynomial. Default 1.
#' @return A list with `scales`, `fluctuation` (one per scale) and `alpha`,
#'   the slope of log fluctuation on log scale.
#' @references Peng, C.-K., Havlin, S., Stanley, H. E. and Goldberger, A. L.
#'   (1995). Quantification of scaling exponents and crossover phenomena in
#'   nonstationary heartbeat time series. Chaos, 5(1), 82-87.
#' @export
dfa <- function(x, scales, order = 1) {
  x <- as.numeric(x)
  profile <- cumsum(x - mean(x))
  fluctuation <- numeric(length(scales))
  for (s in seq_along(scales)) {
    scale <- scales[s]
    chunks <- floor(length(x) / scale)
    idx <- seq_len(scale)
    V <- outer(idx, order:0, `^`)          # Vandermonde, as MATLAB's polyfit
    qrV <- qr(V)
    ms <- 0
    for (i in seq_len(chunks)) {
      seg <- profile[((i - 1) * scale + 1):(i * scale)]
      ms <- ms + mean(qr.resid(qrV, seg)^2)
    }
    fluctuation[s] <- sqrt(ms / chunks)
  }
  ls <- log(scales)
  lf <- log(fluctuation)
  keep <- !is.nan(lf)
  alpha <- unname(stats::coef(stats::lm.fit(cbind(ls[keep], 1), lf[keep]))[1])
  list(scales = scales, fluctuation = fluctuation, alpha = alpha)
}
