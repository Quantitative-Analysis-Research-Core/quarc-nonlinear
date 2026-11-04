# fGn is random, so it is tested by what it promises: unit variance and the
# DFA slope recovering H.

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
