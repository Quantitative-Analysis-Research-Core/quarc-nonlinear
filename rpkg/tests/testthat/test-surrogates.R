# Surrogates and fGn are random, and R cannot draw MATLAB's numbers, so they
# are tested by what each algorithm promises to preserve -- the same
# contracts as the MATLAB suite (quarctest.surrogateContract).

amp <- function(v) Mod(stats::fft(v - mean(v)))

test_that("algorithm 1 keeps the amplitude spectrum exactly", {
  set.seed(1)
  x <- cumsum(rnorm(512)) + sin(1:512 / 5)
  for (i in 1:5) {
    z <- surr_theiler(x, 1)
    expect_equal(amp(z), amp(x), tolerance = 1e-9)
    expect_equal(mean(z), mean(x), tolerance = 1e-12)
    expect_false(isTRUE(all.equal(z, x)))
  }
})

test_that("algorithms 0 and 2 keep the distribution exactly", {
  set.seed(2)
  x <- rexp(400)^2
  for (alg in c(0, 2)) {
    z <- surr_theiler(x, alg)
    expect_identical(sort(z), sort(x), label = paste("algorithm", alg))
    expect_false(identical(z, x))
  }
})

test_that("algorithm 2's spectrum is close but not exact", {
  set.seed(3)
  x <- sin(1:1024 / 8) + 0.3 * rnorm(1024)
  z <- surr_theiler(x, 2)
  r <- stats::cor(amp(z)[2:200], amp(x)[2:200])
  expect_gt(r, 0.9)
})

test_that("surrogates keep the shape of their input and reject unknown algorithms", {
  set.seed(4)
  m <- cbind(rnorm(100), rnorm(100))
  for (alg in 0:2) expect_equal(dim(surr_theiler(m, alg)), c(100L, 2L))
  expect_length(surr_theiler(rnorm(99), 1), 99)
  expect_error(surr_theiler(rnorm(50), 3), "Unknown algorithm")
})

test_that("pseudo-periodic surrogates copy values from the series", {
  set.seed(5)
  x <- sin(2 * pi * (1:600) / 40) + 0.05 * rnorm(600)
  s <- surr_pseudo_periodic(x, 10, 3, 0.2)
  expect_length(s$ys, 600 - 20)
  expect_equal(s$ys, x[s$yi])
  expect_true(all(s$yi >= 1 & s$yi <= 580))
})

test_that("surr_find_rho returns a rho inside its bounds", {
  set.seed(6)
  x <- sin(2 * pi * (1:400) / 40) + 0.05 * rnorm(400)
  r <- suppressWarnings(surr_find_rho(x, 10, 3, 0.05, 1))
  expect_gte(r$rho, 0.05)
  expect_lte(r$rho, 1)
  expect_error(surr_find_rho(x, 10, 3, 1, 0.5), "0 < rhoL < rhoH")
})

test_that("fgn_sim gives unit-variance noise whose DFA slope recovers H", {
  set.seed(7)
  for (H in c(0.3, 0.5, 0.85)) {
    x <- fgn_sim(8192, H)
    expect_length(x, 8192)
    expect_equal(stats::sd(x), 1, tolerance = 0.1)
    sc <- unique(round(2^seq(4, 10, by = 0.5)))
    expect_equal(dfa(x, sc)$alpha, H, tolerance = 0.08 / H, label = paste("H =", H))
  }
})
