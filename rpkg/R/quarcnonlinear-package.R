#' quarcnonlinear: nonlinear time series analysis
#'
#' R port of the QUARC Nonlinear MATLAB library. Function names, argument
#' names and defaults follow the MATLAB functions of the same name, which are
#' the specification (see `PORTING.md` in the repository).
#'
#' Where a MATLAB function returns several outputs, the R function returns a
#' named list whose element names are the MATLAB output names.
#'
#' @keywords internal
#' @useDynLib quarcnonlinear, .registration = TRUE
#' @importFrom Rcpp sourceCpp
"_PACKAGE"
