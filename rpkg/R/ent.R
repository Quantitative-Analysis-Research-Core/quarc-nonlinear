#' Entropy of a time series, by algorithm
#'
#' One entry point for the entropy family.
#'
#' @param x Numeric vector, no `NA`.
#' @param algorithm `"sample"` (default), `"approximate"`, `"permutation"`,
#'   `"symbolic"` (binary series only) or `"multiscale"` (refined composite
#'   multiscale sample entropy).
#' @param y Optional second series: gives the cross entropy. Only with
#'   `"sample"` or `"approximate"`.
#' @param dim Template length or pattern order. Default 2.
#' @param radius Tolerance, as a proportion of the SD. Default 0.2.
#' @param delay Delay for `"permutation"`. Default 1.
#' @param scale Largest scale for `"multiscale"`. Default 20.
#' @return A list with `value` and `extra`, a list of diagnostics.
#' @export
ent <- function(x, algorithm = "sample", y = NULL, dim = 2, radius = 0.2, delay = 1,
                scale = 20) {
  x <- as.numeric(x)
  if (anyNA(x)) stop("x contains NA. Entropy is undefined for missing data; remove or impute first.",
                     call. = FALSE)
  if (length(x) <= dim + 1) {
    stop(sprintf("Series of %d samples is too short for dim = %d.", length(x), as.integer(dim)),
         call. = FALSE)
  }
  isCross <- !is.null(y)
  if (isCross) {
    y <- as.numeric(y)
    if (anyNA(y)) stop("y contains NA.", call. = FALSE)
    if (length(y) != length(x)) {
      stop(sprintf("x has %d samples and y has %d; a cross entropy needs equal lengths.",
                   length(x), length(y)), call. = FALSE)
    }
  }
  algorithm <- tolower(algorithm)
  rejectCross <- function(name) {
    if (isCross) {
      stop(sprintf(paste0('algorithm "%s" has no cross-entropy form. Supply y only with ',
                          '"sample" or "approximate".'), name), call. = FALSE)
    }
  }
  extra <- list(algorithm = algorithm, dim = dim, cross = isCross)
  value <- switch(algorithm,
    sample = {
      extra$radius <- radius
      extra$estimator <- if (isCross) "cross sample" else "sample"
      if (isCross) ent_xsamp(x, y, dim, radius, 1) else ent_samp(x, dim, radius)
    },
    approximate = {
      extra$radius <- radius
      extra$estimator <- if (isCross) "cross approximate" else "approximate"
      if (isCross) ent_xap(x, y, dim, radius, 1) else ent_ap(x, dim, radius)
    },
    permutation = {
      rejectCross("permutation")
      p <- ent_permu(x, dim, delay)
      extra$delay <- delay
      extra$patternCounts <- p$hist
      extra$estimator <- "permutation"
      p$permEnt
    },
    symbolic = {
      rejectCross("symbolic")
      u <- unique(x)
      if (!all(u %in% c(0, 1))) {
        stop(sprintf(paste0('algorithm "symbolic" needs a binary series; x has %d distinct ',
                            "values in [%g, %g]. Threshold it first, for example ",
                            "x = as.numeric(x > median(x))."), length(u), min(u), max(u)),
             call. = FALSE)
      }
      extra$estimator <- "symbolic"
      ent_symbolic(x, dim)
    },
    multiscale = {
      rejectCross("multiscale")
      ms <- ent_ms_plus(x, scale, dim, radius)
      extra$radius <- radius
      extra$scale <- scale
      extra$rcmse <- ms$RCMSE
      extra$cmse <- ms$CMSE
      extra$mse <- ms$MSE
      extra$msfe <- ms$MSFE
      extra$gmse <- ms$GMSE
      extra$estimator <- "refined composite multiscale sample"
      ms$RCMSE
    },
    stop(sprintf(paste0('Unknown algorithm "%s". Supported: "sample", "approximate", ',
                        '"permutation", "symbolic", "multiscale".'), algorithm), call. = FALSE))
  list(value = value, extra = extra)
}
