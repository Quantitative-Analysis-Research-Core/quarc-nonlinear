#' @rdname rqa
#' @export
jrqa <- function(x, delay = c(1, 1), dim = c(1, 1), param = "rec", threshold = 2.5,
                 zscore = 1, norm = "none", dmin = 2, vmin = 2, iter = 20,
                 phasespace = FALSE) {
  if (phasespace) {
    if (!is.list(x)) {
      stop(paste0("PhaseSpace is true, so DATA must be a list with one already-reconstructed ",
                  "phase space per channel, e.g. list(Y1, Y2), rather than a raw multi-column series."),
           call. = FALSE)
    }
    if (length(x) < 2) stop("DATA must contain at least two channels.", call. = FALSE)
  } else {
    if (is.list(x) && !is.data.frame(x)) {
      stop("DATA must be a numeric matrix unless PhaseSpace is true.", call. = FALSE)
    }
    x <- as.matrix(x)
    if (ncol(x) < 2) stop("Data must be at least two column vectors.", call. = FALSE)
  }
  if (zscore) x <- if (phasespace) lapply(x, zscore_cols) else zscore_cols(x)
  if (phasespace) {
    data2 <- lapply(x, as.matrix)
    if (norm == "euc") {
      warning(paste0("PhaseSpace is true, so DATA was not rebuilt from TAU and DIM and DIM ",
                     'cannot be confirmed to match its width. Norm="euc" scales by DIM ',
                     "regardless; pass DIM to match each phase space's width or use ",
                     'Norm="none", "min" or "max".'), call. = FALSE)
    }
  } else {
    data2 <- lapply(seq_len(ncol(x)), function(i) {
      if (dim[i] > 1) psr(x[, i], delay[i], dim[i]) else x[, i, drop = FALSE]
    })
  }
  DIM <- length(data2)
  if (norm == "euc" && length(unique(dim)) > 1) {
    stop(paste0('Norm="euc" scales by DIM, which is per channel here and differs between ',
                'channels. Use Norm="none", "min" or "max".'), call. = FALSE)
  }
  a <- lapply(data2, function(d) rqa_normalise(-pairwise_dist(d, d), norm, DIM * dim[1]))
  rqa_finish(x, a, a[[1]], param, threshold, "jrqa", iter, 0, dmin, vmin,
             list(DIM = 1, EMB = dim, DEL = delay, NORM = norm, ZSCORE = zscore), FALSE)
}
