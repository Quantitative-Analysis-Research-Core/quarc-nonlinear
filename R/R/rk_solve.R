rk_solve <- function(f, tspan, y0, rtol, atol, tab) {
  tspan <- as.numeric(tspan)
  y <- as.numeric(y0)
  t0 <- tspan[1]
  tfinal <- tspan[length(tspan)]
  tdir <- sign(tfinal - t0)
  rtol <- max(rtol, 100 * .Machine$double.eps)
  threshold <- atol / rtol
  pow <- tab$pow
  nst <- length(tab$A) + 1
  refineMode <- length(tspan) == 2
  hmax <- abs(0.1 * (tfinal - t0))
  t <- t0
  f1 <- f(t, y)
  K <- matrix(0, length(y), nst)
  K[, 1] <- f1
  htspan <- abs(tspan[2] - tspan[1])
  absh <- min(hmax, htspan)
  rh <- max(abs(f1 / pmax(abs(y), threshold))) / (0.8 * rtol^pow)
  if (absh * rh > 1) absh <- 1 / rh
  eps_at <- function(v) { v <- abs(v); if (v == 0) 2^-1074 else 2^(floor(log2(v)) - 52) }
  absh <- max(absh, 16 * eps_at(t))

  tout <- t0
  yout <- list(y)
  nextIdx <- 2
  done <- FALSE
  while (!done) {
    hmin <- 16 * eps_at(t)
    absh <- min(hmax, max(hmin, absh))
    h <- tdir * absh
    if (1.1 * absh >= abs(tfinal - t)) {
      h <- tfinal - t
      absh <- abs(h)
      done <- TRUE
    }
    nofailed <- TRUE
    repeat {
      hB <- h * tab$B
      for (j in 2:(nst - 1)) {
        K[, j] <- f(t + h * tab$A[j - 1], y + K[, 1:(j - 1), drop = FALSE] %*% hB[1:(j - 1), j - 1])
      }
      tnew <- t + h * tab$A[nst - 1]
      if (done) tnew <- tfinal
      h <- tnew - t
      ynew <- as.numeric(y + K[, 1:(nst - 1), drop = FALSE] %*% hB[1:(nst - 1), nst - 1])
      K[, nst] <- f(tnew, ynew)
      err <- absh * max(abs(as.numeric(K %*% tab$E) / pmax(pmax(abs(y), abs(ynew)), threshold)))
      if (err > rtol) {
        if (absh <= hmin) stop(sprintf("Unable to meet integration tolerances at t = %g.", t),
                               call. = FALSE)
        if (nofailed) {
          nofailed <- FALSE
          absh <- max(hmin, absh * max(if (pow == 1 / 5) 0.1 else 0.5, 0.8 * (rtol / err)^pow))
        } else {
          absh <- max(hmin, 0.5 * absh)
        }
        h <- tdir * absh
        done <- FALSE
      } else {
        break
      }
    }
    interp <- function(ti) {
      s <- (ti - t) / h
      y + as.numeric(K %*% (h * tab$BI) %*% cumprod(rep(s, ncol(tab$BI))))
    }
    if (refineMode) {
      for (ti in t + h * (seq_len(tab$refine - 1) / tab$refine)) {
        tout <- c(tout, ti)
        yout[[length(yout) + 1]] <- interp(ti)
      }
      tout <- c(tout, tnew)
      yout[[length(yout) + 1]] <- ynew
    } else {
      while (nextIdx <= length(tspan) && tdir * (tnew - tspan[nextIdx]) >= 0) {
        ti <- tspan[nextIdx]
        tout <- c(tout, ti)
        yout[[length(yout) + 1]] <- if (ti == tnew) ynew else interp(ti)
        nextIdx <- nextIdx + 1
      }
    }
    if (nofailed) {
      temp <- 1.25 * (err / rtol)^pow
      absh <- if (temp > 0.2) absh / temp else 5 * absh
    }
    t <- tnew
    y <- ynew
    K[, 1] <- K[, nst]
  }
  list(t = tout, y = do.call(rbind, yout))
}
