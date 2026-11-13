#' @rdname rqa
#' @export
mdrqa <- function(x, delay = 1, dim = 1, param = "rec", threshold = 2.5, zscore = 1,
                  norm = "none", dmin = 2, vmin = 2, iter = 20, phasespace = FALSE) {
  x <- as.matrix(x)
  if (zscore) x <- zscore_cols(x)
  if (phasespace) {
    if (norm == "euc") {
      warning(paste0("PhaseSpace is true, so DATA was not rebuilt from TAU and DIM and DIM ",
                     'cannot be confirmed to match its width. Norm="euc" scales by DIM ',
                     'regardless; pass DIM equal to ncol(data) or use Norm="none", "min" or "max".'),
              call. = FALSE)
    }
  } else if (dim > 1) {
    x <- psr(x, delay, dim)
  }
  a <- rqa_normalise(-pairwise_dist(x, x), norm, 1 * dim)
  rqa_finish(x, a, a, param, threshold, "mdrqa", iter, 0, dmin, vmin,
             list(DIM = 1, EMB = dim, DEL = delay, NORM = norm, ZSCORE = zscore), TRUE)
}
