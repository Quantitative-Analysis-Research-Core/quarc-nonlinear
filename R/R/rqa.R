#' Recurrence quantification analysis
#'
#' Auto-recurrence of one series (`rqa`), cross recurrence of two (`crqa`),
#' joint recurrence of several (`jrqa`) and multidimensional recurrence of
#' several channels in one phase space (`mdrqa`).
#'
#' @param x For `rqa`, a numeric vector (or a phase space with
#'   `phasespace = TRUE`). For `crqa`, a two-column matrix. For `jrqa`, a
#'   matrix with one channel per column, or with `phasespace = TRUE` a list
#'   of phase spaces. For `mdrqa`, a matrix with one channel per column.
#' @param delay,dim Embedding delay and dimension. For `jrqa`, one per
#'   channel. A channel with `dim = 1` is used unembedded.
#' @param param `"rec"` sets the radius to reach `threshold` percent
#'   recurrence; `"rad"` uses `threshold` as the radius.
#' @param threshold Target percent recurrence, or the radius.
#' @param zscore 1 (default) standardises the data first.
#' @param norm Distance normalisation: `"none"` (default), `"euc"`, `"min"`
#'   or `"max"`.
#' @param dmin,vmin Shortest diagonal and vertical line counted. Default 2.
#' @param iter Bisection steps for `param = "rec"`. Default 20.
#' @param phasespace `TRUE` if `x` is already a phase space.
#' @param theiler Theiler window, `rqa` only. Default 0.
#' @return A list with `RP`, the recurrence plot as displayed (rotated as
#'   MATLAB's `imrotate(1 - A, 90)`, so 0 marks a recurrence), and
#'   `RESULTS`: `DIM`, `EMB`, `DEL`, `RADIUS`, `NORM`, `ZSCORE`, `Size`, `REC`,
#'   `DET`, `MeanL`, `MaxL`, `EntrL`, `LAM`, `MeanV`, `MaxV`, `EntrV` and
#'   (except `jrqa`) `EntrW`.
#' @references Marwan, N., Romano, M. C., Thiel, M. and Kurths, J. (2007).
#'   Recurrence plots for the analysis of complex systems. Physics Reports,
#'   438(5-6), 237-329.
#' @export
rqa <- function(x, delay = 1, dim = 1, param = "rec", threshold = 2.5, zscore = 1,
                norm = "none", dmin = 2, vmin = 2, iter = 20, phasespace = FALSE,
                theiler = 0) {
  x <- as.matrix(x)
  if (!phasespace && ncol(x) != 1) stop("Data must be single column vector.", call. = FALSE)
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
  rqa_finish(x, a, a, param, threshold, "rqa", iter, theiler, dmin, vmin,
             list(DIM = 1, EMB = dim, DEL = delay, NORM = norm, ZSCORE = zscore), TRUE)
}
