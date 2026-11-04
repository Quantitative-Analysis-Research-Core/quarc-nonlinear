check_ami_input <- function(x, maxlag, caller) {
  x <- as.numeric(x)
  if (anyNA(x)) {
    stop(sprintf(paste0("%s: x contains NA. Mutual information is undefined for ",
                        "missing data; remove or impute first."), caller), call. = FALSE)
  }
  if (maxlag >= length(x)) {
    stop(sprintf("%s: maxlag (%d) must be smaller than the number of samples (%d).",
                 caller, as.integer(maxlag), length(x)), call. = FALSE)
  }
  x
}
