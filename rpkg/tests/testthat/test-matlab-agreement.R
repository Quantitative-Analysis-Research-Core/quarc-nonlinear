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

test_that("cross, symbolic and weighted entropies match MATLAB", {
  skip_without_fixtures()
  ref <- reference()
  a <- fixture("ar1_phi70_512")
  h <- fixture("henon_512")
  expect_equal(c(ent_xsamp(a, h, 2, 0.2, 0), ent_xsamp(a, h, 2, 0.2, 1),
                 ent_xsamp(a, h, 2, 0.2, 2), ent_xsamp(a, h, 1, 0.3, 2)),
               ref$ent_xsamp$ar1_henon, tolerance = TOL)
  expect_equal(c(ent_xap(a, h, 2, 0.2, 1), ent_xap(a, h, 1, 0.3, 1), ent_xap(a, h, 2, 0.25, 2)),
               ref$ent_xap$ar1_henon, tolerance = TOL)
  b <- as.numeric(a > stats::median(a))
  expect_equal(c(ent_symbolic(b, 2), ent_symbolic(b, 3), ent_symbolic(b, 4)),
               ref$ent_symbolic$ar1_median, tolerance = TOL)
  y <- fixture("lorenz_2048")
  Y <- psr(y[1:400], 10, 3)
  D <- as.matrix(stats::dist(Y))
  expect_equal(ent_weighted(exp(-D / stats::sd(y))), ref$ent_weighted$lorenz_400,
               tolerance = TOL)
})

test_that("ent_ms_plus matches MATLAB on all five outputs", {
  skip_without_fixtures()
  ref <- reference()$ent_ms_plus
  for (s in names(ref)) {
    r <- ref[[s]]
    out <- ent_ms_plus(fixture(s), r$scale, r$dim, r$radius)
    for (f in c("rcmse", "cmse", "mse", "msfe", "gmse")) {
      expect_equal(out[[toupper(f)]], r[[f]], tolerance = TOL, label = paste(s, f))
    }
  }
})

test_that("fnn matches MATLAB, k-d tree included", {
  skip_without_fixtures()
  ref <- reference()$fnn
  for (s in names(ref)) {
    r <- ref[[s]]
    x <- fixture(s)
    a <- suppressMessages(fnn(x, r$delay, r$maxdim, 15, 2, 0))
    b <- suppressMessages(fnn(x, r$delay, r$maxdim, 10, 2, 1))
    expect_equal(a$dE, r$dE_speed0, tolerance = TOL, label = paste(s, "speed 0"))
    expect_equal(a$dim, r$dim_speed0, label = paste(s, "speed 0"))
    expect_equal(b$dE, r$dE_speed1, tolerance = TOL, label = paste(s, "speed 1"))
    expect_equal(b$dim, r$dim_speed1, label = paste(s, "speed 1"))
  }
})

test_that("corr_dim matches MATLAB: curve, fitted region and slope", {
  skip_without_fixtures()
  ref <- reference()$corr_dim
  for (s in names(ref)) {
    r <- ref[[s]]
    out <- corr_dim(fixture(s), r$delay, r$dim)
    expect_equal(out$extra$logEps, r$logEps, tolerance = TOL, label = s)
    expected <- r$logC
    expected[is.na(expected)] <- -Inf       # MATLAB's jsonencode writes -Inf as null
    expect_equal(out$extra$logC, expected, tolerance = TOL, label = s)
    expect_equal(out$extra$idx, r$idx, label = s)
    expect_equal(out$CoD, r$CoD, tolerance = TOL, label = s)
    expect_equal(out$extra$r2, r$r2, tolerance = TOL, label = s)
  }
})

test_that("the RQA family matches MATLAB: plot and every measure", {
  skip_without_fixtures()
  ref <- reference()
  lz <- fixture("lorenz_2048")
  hn <- fixture("henon_512")
  ar <- fixture("ar1_phi70_512")
  q <- function(f) suppressMessages(f)
  cases <- list(
    rqa = list(
      lorenz_rec = q(rqa(lz[1:600], 10, 3, "rec", 2.5)),
      lorenz_rad_euc_theiler = rqa(lz[1:600], 10, 3, "rad", 0.3, norm = "euc", theiler = 3),
      henon_rec_max = q(rqa(hn, 1, 2, "rec", 5, norm = "max")),
      ar1_rad_min = rqa(ar, 1, 1, "rad", 0.5, norm = "min", dmin = 3, vmin = 3, zscore = 0)),
    crqa = list(
      ar1_henon_rec = q(crqa(cbind(ar[1:400], hn[1:400]), 1, 2, "rec", 3)),
      lorenz_rad_euc = crqa(cbind(lz[1:400], lz[6:405]), 4, 3, "rad", 0.4, norm = "euc")),
    jrqa = list(
      lorenz_rec = q(jrqa(cbind(lz[1:300], lz[11:310]), c(5, 5), c(2, 2), "rec", 5)),
      ar1_henon_rad = jrqa(cbind(ar[1:300], hn[1:300]), c(1, 1), c(1, 1), "rad", 0.5, norm = "max")),
    mdrqa = list(
      lorenz_henon_rec = q(mdrqa(cbind(lz[1:400], hn[1:400]), 1, 1, "rec", 4)),
      lorenz_henon_rad = mdrqa(cbind(lz[1:400], hn[1:400]), 2, 2, "rad", 0.5, norm = "max")))
  fields <- c("REC", "RADIUS", "Size", "DET", "MeanL", "MaxL", "EntrL", "LAM", "MeanV",
              "MaxV", "EntrV", "EntrW")
  na_is_nan <- function(v) { v[is.null(v)] <- NaN; if (is.null(v) || is.na(v)) NaN else v }
  for (fam in names(cases)) for (s in names(cases[[fam]])) {
    out <- cases[[fam]][[s]]
    r <- ref[[fam]][[s]]
    lab <- paste(fam, s)
    expect_equal(rowSums(out$RP), as.numeric(r$rowsum), label = paste(lab, "RP rows"))
    expect_equal(colSums(out$RP), as.numeric(r$colsum), label = paste(lab, "RP cols"))
    for (f in fields) {
      got <- out$RESULTS[[f]]
      want <- if (is.null(r[[f]])) NaN else r[[f]]
      if (is.null(got)) got <- NaN                    # jrqa has no EntrW
      expect_equal(got, want, tolerance = TOL, label = paste(lab, f))
    }
  }
})

test_that("lye_r matches MATLAB: neighbours and divergence curve", {
  skip_without_fixtures()
  ref <- reference()$lye_r
  lz <- fixture("lorenz_2048")
  cases <- list(lorenz = lye_r(lz, 100, 10, 3),
                lorenz_theiler25 = lye_r(lz, 100, 10, 3, theiler = 25),
                henon = lye_r(fixture("henon_512"), 1, 1, 2))
  for (s in names(cases)) {
    expect_equal(cases[[s]][, 2], as.numeric(ref[[s]]$pairs), label = s)
    expect_equal(cases[[s]][, 3], ref[[s]]$divergence, tolerance = TOL, label = s)
  }
  sl <- lye_r(lz, 100, 10, 3, slope = c(0, 0.4, 1, 3), MeanPeriod = 0.75)
  expect_equal(sl$LyES, ref$lorenz_slopes$LyES, tolerance = TOL)
  expect_equal(sl$LyEL, ref$lorenz_slopes$LyEL, tolerance = TOL)
})

test_that("lye_w matches MATLAB: every replacement step", {
  skip_without_fixtures()
  ref <- reference()$lye_w
  lz <- fixture("lorenz_2048")
  cases <- list(lorenz_e10 = lye_w(lz, 100, 10, 3, 10),
                henon_e1 = lye_w(fixture("henon_512"), 1, 1, 2, 1),
                lorenz_custom = lye_w(lz[1:1500], 100, 8, 4, 5, 2, 20 * pi / 180, 1))
  for (s in names(cases)) {
    want <- matrix(ref[[s]]$out, ncol = ref[[s]]$ncol, byrow = TRUE)
    expect_equal(cases[[s]]$out, want, tolerance = TOL, label = s)
    expect_equal(cases[[s]]$LyE, ref[[s]]$LyE, tolerance = TOL, label = s)
  }
})

test_that("lyapunov matches MATLAB for both methods", {
  skip_without_fixtures()
  ref <- reference()$lyapunov
  lz <- fixture("lorenz_2048")
  r <- lyapunov(lz, 100, delay = 10, dim = 3)
  expect_equal(r$lambda, ref$lorenz_rosenstein$lambda, tolerance = TOL)
  expect_equal(r$extra$scalingRegion, as.integer(ref$lorenz_rosenstein$scalingRegion))
  expect_equal(r$extra$fitR2, ref$lorenz_rosenstein$fitR2, tolerance = TOL)
  h <- lyapunov(fixture("henon_512"), 1, delay = 1, dim = 2)
  expect_equal(h$lambda, ref$henon_rosenstein$lambda, tolerance = TOL)
  expect_equal(h$extra$scalingRegion, as.integer(ref$henon_rosenstein$scalingRegion))
  w <- lyapunov(lz, 100, algorithm = "wolf", delay = 10, dim = 3, evolve = 10)
  expect_equal(w$lambda, ref$lorenz_wolf$lambda, tolerance = TOL)
})
