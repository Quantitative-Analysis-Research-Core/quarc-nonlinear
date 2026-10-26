# The compiled kernels (cpp/) against the R code they replace, and the
# package's copy of the core against cpp/ itself.

test_that("compiled ent_samp agrees with the R reference code", {
  set.seed(1)
  for (n in c(10, 200, 1000)) {
    x <- cumsum(rnorm(n))
    for (dim in 1:3) {
      r <- 0.2 * stats::sd(x)
      expect_equal(ent_samp(x, dim, 0.2), quarcnonlinear:::ent_samp_r(x, dim, r),
                   tolerance = 1e-12, label = sprintf("n=%d dim=%d", n, dim))
    }
  }
})

test_that("rpkg/src carries an up-to-date copy of cpp/", {
  skip_without_fixtures()
  core <- file.path(dirname(dirname(fixture_dir())), "cpp")
  skip_if(!dir.exists(core), "cpp/ not found (not in the repository)")
  src <- file.path(dirname(dirname(fixture_dir())), "rpkg", "src")
  for (f in list.files(file.path(core, "src"), pattern = "\\.cpp$")) {
    expect_identical(readLines(file.path(src, paste0("core_", f))),
                     readLines(file.path(core, "src", f)), label = f)
  }
  for (f in list.files(file.path(core, "include", "quarc"))) {
    expect_identical(readLines(file.path(src, "quarc", f)),
                     readLines(file.path(core, "include", "quarc", f)), label = f)
  }
})
