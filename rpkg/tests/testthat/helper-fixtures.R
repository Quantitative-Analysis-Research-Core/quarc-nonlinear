# The shared fixtures live in the repository's tests/fixtures/, outside the
# package. Walk up from the test directory to find them. When the package is
# tested away from the repository (R CMD check on a built tarball) they are
# not there, and the agreement tests skip rather than fail.
fixture_dir <- function() {
  d <- normalizePath(testthat::test_path("."), mustWork = FALSE)
  for (i in 1:6) {
    cand <- file.path(d, "tests", "fixtures")
    if (file.exists(file.path(cand, "matlab_reference.json"))) return(cand)
    d <- dirname(d)
  }
  NULL
}

skip_without_fixtures <- function() {
  testthat::skip_if_not_installed("jsonlite")
  testthat::skip_if(is.null(fixture_dir()), "shared fixtures not found (not in the repository)")
}

reference <- function() {
  jsonlite::fromJSON(file.path(fixture_dir(), "matlab_reference.json"), simplifyVector = TRUE)
}

fixture <- function(name) {
  as.numeric(utils::read.csv(file.path(fixture_dir(), paste0(name, ".csv")),
                             header = FALSE)[[1]])
}

# PORTING.md: differences at 1e-15 are summation order; 1e-4 is a defect.
# 1e-10 relative leaves room for the first and none for the second.
TOL <- 1e-10
