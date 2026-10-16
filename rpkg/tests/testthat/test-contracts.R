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

test_that("ent dispatches to each estimator and guards its inputs", {
  x <- sin(1:300 / 4) + 0.3 * cos(1:300 / 1.3)
  y <- cos(1:300 / 4)
  expect_equal(ent(x)$value, ent_samp(x, 2, 0.2))
  expect_equal(ent(x, "approximate")$value, ent_ap(x, 2, 0.2))
  expect_equal(ent(x, "permutation", dim = 3)$value, ent_permu(x, 3, 1)$permEnt)
  expect_equal(ent(x, y = y)$value, ent_xsamp(x, y, 2, 0.2, 1))
  expect_equal(ent(x, "multiscale", scale = 3)$value, ent_ms_plus(x, 3, 2, 0.2)$RCMSE)
  expect_error(ent(x, "permutation", y = y), "no cross-entropy form")
  expect_error(ent(x, "symbolic"), "binary")
  expect_error(ent(x, y = y[-1]), "equal lengths")
  expect_error(ent(c(1, NA, 3, 4, 5)), "contains NA")
  expect_error(ent(x, "fuzzy"), "Unknown algorithm")
})

test_that("rqa hits its recurrence target and guards its inputs", {
  x <- sin(1:400 / 6) + 0.2 * cos(1:400 / 1.9)
  out <- suppressMessages(rqa(x, 6, 3, "rec", 2.5))
  expect_equal(out$RESULTS$REC, 2.5, tolerance = 0.02)
  expect_equal(dim(out$RP), rep(400 - 12, 2))
  expect_error(rqa(cbind(x, x)), "single column")
  expect_error(crqa(x), "two column")
  expect_error(jrqa(x), "at least two")
  expect_error(jrqa(cbind(x, x), c(1, 1), c(1, 2), "rad", 0.5, norm = "euc"), "differs between")
})

test_that("lyapunov requires evolve for wolf and guards its inputs", {
  x <- sin(1:800 / 7) + 0.3 * sin(1:800 / 2.1)
  expect_error(lyapunov(x, 1, algorithm = "wolf", delay = 5, dim = 3), "needs evolve")
  expect_error(lyapunov(c(x[1:10], NA), 1), "contains NA")
  expect_error(lyapunov(x, 1, algorithm = "ftle"), "Unknown algorithm")
  expect_warning(lyapunov(psr(x, 5, 3), 1), "without theiler")
  expect_equal(lyapunov(x, 1, algorithm = "w", delay = 5, dim = 3, evolve = 3)$extra$estimator,
               "wolf")
})

test_that("chaos_library names its failures and the solvers reach tfinal", {
  expect_error(chaos_library("Chua", 1:10, 1, 1), "Unknown system")
  expect_error(chaos_library("Logistic", 0, 0.2, 3.9), "at least one")
  # Exponential decay has an exact answer.
  s <- ode45(function(t, y) -y, c(0, 2), 1, rtol = 1e-8, atol = 1e-10)
  expect_equal(tail(s$t, 1), 2)
  expect_equal(tail(s$y[, 1], 1), exp(-2), tolerance = 1e-7)
  s3 <- ode23(function(t, y) -y, seq(0, 1, by = 0.1), 1, rtol = 1e-6, atol = 1e-9)
  expect_equal(s3$y[, 1], exp(-seq(0, 1, by = 0.1)), tolerance = 1e-5)
})
