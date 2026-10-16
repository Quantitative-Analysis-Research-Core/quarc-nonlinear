# Every number here was produced by MATLAB on the same samples
# (tests/matlab/make_reference.m). Agreement is the definition of a correct
# port: see PORTING.md.

test_that("psr matches MATLAB", {
  skip_without_fixtures()
  r <- reference()$psr$lorenz_2048_d10_m3
  Y <- psr(fixture("lorenz_2048"), 10, 3)
  expect_equal(dim(Y), as.integer(r$size))
  expect_equal(Y[1, ], r$first, tolerance = TOL)
  expect_equal(Y[nrow(Y), ], r$last, tolerance = TOL)
  expect_equal(colSums(Y), r$colsum, tolerance = TOL)
})

test_that("ent_samp matches MATLAB", {
  skip_without_fixtures()
  ref <- reference()
  expect_equal(ent_samp(fixture("ar1_phi70_512"), 2, 0.2), ref$ent_samp$ar1_phi70_512,
               tolerance = TOL)
  for (s in names(ref$ent_samp_dims)) {
    x <- fixture(s)
    got <- c(ent_samp(x, 1, 0.2), ent_samp(x, 2, 0.2), ent_samp(x, 3, 0.15))
    expect_equal(got, ref$ent_samp_dims[[s]], tolerance = TOL, label = s)
  }
})

test_that("ent_ap matches MATLAB", {
  skip_without_fixtures()
  ref <- reference()$ent_ap
  for (s in names(ref)) {
    x <- fixture(s)
    got <- c(ent_ap(x, 1, 0.2), ent_ap(x, 2, 0.2), ent_ap(x, 3, 0.15))
    expect_equal(got, ref[[s]], tolerance = TOL, label = s)
  }
})

test_that("ent_permu matches MATLAB, value and pattern histogram", {
  skip_without_fixtures()
  ref <- reference()$ent_permu
  x <- fixture("lorenz_2048")
  a <- ent_permu(x, 3, 1)
  b <- ent_permu(x, 4, 2)
  expect_equal(a$permEnt, ref$lorenz_2048$dim3_delay1, tolerance = TOL)
  expect_equal(a$hist, as.numeric(ref$lorenz_2048$hist3))
  expect_equal(b$permEnt, ref$lorenz_2048$dim4_delay2, tolerance = TOL)
  expect_equal(b$hist, as.numeric(ref$lorenz_2048$hist4))
  w <- ent_permu(fixture("white_4096"), 4, 1)
  expect_equal(w$permEnt, ref$white_4096$dim4_delay1, tolerance = TOL)
  expect_equal(w$hist, as.numeric(ref$white_4096$hist4))
})

test_that("dfa matches MATLAB at every scale", {
  skip_without_fixtures()
  ref <- reference()$dfa
  for (s in names(ref)) {
    out <- dfa(fixture(s), ref[[s]]$scales, 1)
    expect_equal(out$fluctuation, ref[[s]]$fluctuation, tolerance = TOL, label = s)
    expect_equal(out$alpha, ref[[s]]$alpha, tolerance = TOL, label = s)
  }
})

test_that("ami_histogram matches MATLAB: curve, tau and bin count", {
  skip_without_fixtures()
  ref <- reference()
  for (s in names(ref$ami_histogram)) {
    x <- fixture(s)
    out <- ami_histogram(x, 20)
    expect_equal(out$info$bins, ref$ami_histogram[[s]]$bins, label = s)
    expect_equal(out$curve[, 2], ref$ami_histogram[[s]]$curve, tolerance = TOL, label = s)
    expect_equal(out$tau, ref$ami_histogram[[s]]$tau, label = s)
    out16 <- ami_histogram(x, 20, bins = 16)
    expect_equal(out16$curve[, 2], ref$ami_histogram_bins16[[s]]$curve, tolerance = TOL,
                 label = paste(s, "bins=16"))
    expect_equal(out16$tau, ref$ami_histogram_bins16[[s]]$tau, label = paste(s, "bins=16"))
  }
})

test_that("ami_kde matches MATLAB: curve and tau", {
  skip_without_fixtures()
  ref <- reference()
  for (s in names(ref$ami_kde)) {
    maxlag <- length(ref$ami_kde[[s]]$curve) - 1
    out <- ami_kde(fixture(s), maxlag)
    expect_equal(out$curve[, 2], ref$ami_kde[[s]]$curve, tolerance = TOL, label = s)
    expect_equal(out$tau, ref$ami_kde[[s]]$tau, label = s)
  }
  # The older single-value reference from the original AMI_Thomas.
  expect_equal(ami(fixture("lorenz_2048"), 20, algorithm = "kde")$tau,
               ref$ami_thomas$lorenz_2048)
})

test_that("ami_partition, ami_knn and ami_copula match MATLAB", {
  skip_without_fixtures()
  ref <- reference()
  cases <- list(
    ami_partition     = function(x) ami_partition(x, 15),
    ami_partition_a01 = function(x) ami_partition(x, 15, alpha = 0.01),
    ami_knn           = function(x) ami_knn(x, 10),
    ami_knn_k3        = function(x) ami_knn(x, 6, k = 3),
    ami_copula        = function(x) ami_copula(x, 30))
  for (key in names(cases)) {
    for (s in names(ref[[key]])) {
      out <- cases[[key]](fixture(s))
      lab <- paste(key, s)
      expect_equal(out$curve[, 2], ref[[key]][[s]]$curve, tolerance = TOL, label = lab)
      expect_equal(out$tau, ref[[key]][[s]]$tau, label = lab)
    }
  }
})

test_that("ami_kde with adaptive bandwidth matches MATLAB", {
  skip_without_fixtures()
  ref <- reference()$ami_kde_adaptive
  for (s in names(ref)) {
    maxlag <- length(ref[[s]]$curve) - 1
    out <- ami_kde(fixture(s), maxlag, bandwidth = "adaptive")
    expect_equal(out$curve[, 2], ref[[s]]$curve, tolerance = TOL, label = s)
    expect_equal(out$tau, ref[[s]]$tau, label = s)
  }
})
