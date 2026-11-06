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

test_that("compiled line_hist agrees with the R reference code", {
  set.seed(2)
  for (n in c(1, 2, 7, 60)) {
    A <- matrix(rbinom(n * n, 1, 0.4), n)
    A <- A * t(A)
    for (theiler in c(0, 1, 3)) {
      for (cross in c(FALSE, TRUE)) {
        th <- if (cross) 0 else theiler
        lab <- sprintf("n=%d theiler=%d cross=%s", n, theiler, cross)
        got <- quarcnonlinear:::.line_hist_core(A, th, cross)
        want <- quarcnonlinear:::line_hist_lines_r(A, th, cross)
        expect_equal(got$diag, as.numeric(want$diag), label = lab)
        expect_equal(got$vert, as.numeric(want$vert), label = lab)
      }
    }
  }
})

test_that("compiled corr_dim pair loops agree with the R reference code", {
  set.seed(3)
  x <- cumsum(rnorm(400))
  for (p in list(c(2, 1), c(3, 4), c(5, 2))) {
    dim <- p[1]; delay <- p[2]
    r <- quarcnonlinear:::.corr_dim_range_core(x, dim, delay)
    expect_equal(r, quarcnonlinear:::corr_dim_pairs_r(x, dim, delay), tolerance = 1e-12)
    e <- exp(seq(log(r[1]), log(r[2]), length.out = 50))
    e[c(1, 50)] <- r
    expect_equal(quarcnonlinear:::.corr_dim_counts_core(x, dim, delay, e),
                 quarcnonlinear:::corr_dim_pairs_r(x, dim, delay, e))
  }
})

test_that("compiled fnn agrees with the R reference code, k-d tree included", {
  set.seed(4)
  x <- cumsum(rnorm(1300))      # n > 512, so the tree has several leaves
  n <- length(x) - 2 * 4
  for (j in 1:4) {
    expect_equal(quarcnonlinear:::.fnn_core(x, 2, j, n, 15, 2, stats::sd(x)),
                 quarcnonlinear:::fnn_false_fraction_r(x, 2, j, n, 15, 2, stats::sd(x)),
                 label = sprintf("j=%d", j))
  }
})

test_that("compiled lye_r agrees with the R reference code", {
  set.seed(5)
  for (n in c(50, 300)) {
    Y <- psr(cumsum(rnorm(n)), 3, 3)
    for (theiler in c(0, 2, 10)) {
      core <- quarcnonlinear:::.lye_r_core(Y, theiler)
      ref <- quarcnonlinear:::lye_r_divergence_r(Y, theiler)
      lab <- sprintf("n=%d theiler=%d", n, theiler)
      expect_equal(core$nn, as.numeric(ref$out[, 2]), label = lab)
      expect_equal(core$avg, ref$out[, 3], tolerance = 1e-12, label = lab)
      expect_equal(core$count, length(ref$AveLnDiv), label = lab)
    }
  }
})

test_that("R/src carries an up-to-date copy of cpp/", {
  skip_without_fixtures()
  core <- file.path(dirname(dirname(fixture_dir())), "cpp")
  skip_if(!dir.exists(core), "cpp/ not found (not in the repository)")
  src <- file.path(dirname(dirname(fixture_dir())), "R", "src")
  for (f in list.files(file.path(core, "src"), pattern = "\\.cpp$")) {
    expect_identical(readLines(file.path(src, paste0("core_", f))),
                     readLines(file.path(core, "src", f)), label = f)
  }
  for (f in list.files(file.path(core, "include", "quarc"))) {
    expect_identical(readLines(file.path(src, "quarc", f)),
                     readLines(file.path(core, "include", "quarc", f)), label = f)
  }
})
