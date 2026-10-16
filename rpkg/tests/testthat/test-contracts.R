# Behaviour that does not need MATLAB to check. PORTING.md: fail loudly
# rather than return a plausible wrong answer, and return the same shape
# whatever the data.

test_that("AMI refuses missing data and impossible lags", {
  x <- sin(1:200 / 5)
  x[17] <- NA
  expect_error(ami_histogram(x, 10), "contains NA")
  expect_error(ami_kde(x, 10), "contains NA")
  expect_error(ami_histogram(rnorm(50), 60), "must be smaller")
})

test_that("ami dispatches, aliases and rejects unknown algorithms", {
  x <- sin(1:600 / 7) + 0.1 * cos(1:600 / 3)
  expect_equal(ami(x, 10, algorithm = "stergiou")$info$estimator, "histogram")
  expect_equal(ami(x, 10, algorithm = "thomas")$info$estimator, "kde")
  expect_error(ami(x, 10, algorithm = "wavelet"), "Unknown algorithm")
  expect_warning(ami(x, 10, algorithm = "kde", bins = 8), "ignored")
})

test_that("ami returns a scalar tau and a fixed-shape curve", {
  x <- sin(1:600 / 7)
  out <- ami(x, 15)
  expect_length(out$tau, 1)
  expect_equal(dim(out$curve), c(16L, 2L))
  expect_true(out$tau %in% out$curve[, 1])
})

test_that("psr refuses a series too short to embed", {
  expect_error(psr(1:5, 3, 3), "too short")
  expect_equal(dim(psr(1:10, 2, 3)), c(6L, 3L))
  expect_equal(dim(psr(cbind(1:10, 11:20), 2, 3)), c(6L, 6L))
})

test_that("ent_samp is Inf when no template matches", {
  expect_equal(ent_samp(c(1, 100, 1e4, 1e6, 1e8), 2, 1e-9, flag = "abs"), Inf)
})

test_that("new AMI estimators dispatch through ami, with aliases", {
  x <- sin(1:600 / 7) + 0.1 * cos(1:600 / 3)
  pairs <- c("kde-adaptive" = "kde-adaptive", adaptive = "adaptive",
             partition = "adaptive", knn = "knn", ksg = "knn",
             copula = "copula", gcmi = "copula")
  for (a in names(pairs)) {
    expect_equal(ami(x, 8, algorithm = a)$info$estimator, unname(pairs[a]), label = a)
  }
  expect_equal(ami(x, 8, algorithm = "knn", k = 5)$info$k, 5)
  expect_warning(ami(x, 8, algorithm = "histogram", k = 5), "ignored")
})

test_that("rank-based estimators ignore monotone transforms exactly", {
  x <- sin(1:800 / 9) + 0.3 * sin(1:800 / 2.3)
  for (f in list(ami_partition, ami_copula)) {
    expect_equal(f(exp(x), 10)$curve, f(x, 10)$curve, tolerance = 1e-12)
  }
})

test_that("ami_knn leaves the global random stream alone", {
  set.seed(42)
  before <- .Random.seed
  ami_knn(sin(1:300 / 5) + 0.2 * cos(1:300 / 1.7), 15)
  expect_identical(.Random.seed, before)
})
