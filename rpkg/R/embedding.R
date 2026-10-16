#' Delay embedding with arbitrary lags
#'
#' @param z Numeric vector.
#' @param v Lags, or, with `w`, the number of lags.
#' @param w Spacing between lags: when given, the lags are `0:(v-1) * w`.
#' @return A matrix with one row per lag, in increasing lag order, and one
#'   column per embedded time. Row `i` is `z` delayed by the `i`-th smallest
#'   lag, so the lag-0 row holds the latest sample of each vector. This is the
#'   transpose of [psr()]'s layout, as in MATLAB.
#' @export
embed <- function(z, v = c(0, 1, 2), w = NULL) {
  if (!is.null(w)) v <- seq(0, w * (v - 1), by = w)
  z <- as.numeric(z)
  lags <- sort(v)
  dim <- length(lags)
  n <- length(z)
  if (n < lags[dim]) stop("Vector is too small to be embedded with the given lags", call. = FALSE)
  m <- n - (lags[dim] - lags[1])
  t <- (1:m) + lags[dim]
  x <- matrix(0, dim, m)
  for (i in seq_len(dim)) x[i, ] <- z[t - lags[i]]
  x
}

#' False nearest neighbours
#'
#' Kennel, Brown and Abarbanel (1992). For each dimension up to `maxdim`, the
#' fraction of nearest neighbours that separate when one more delay
#' coordinate is added. Neighbours are found with a k-d tree of 512-point
#' leaves, reproduced from the MATLAB implementation step for step.
#'
#' @param x Numeric vector.
#' @param delay Embedding delay in samples.
#' @param maxdim Largest dimension tested.
#' @param Rtol Relative distance threshold. Typically 10-15.
#' @param Atol Attractor-size threshold. Typically 2.
#' @param speed `1` stops at the first dimension meeting the criteria; `0`
#'   computes every dimension first, then chooses.
#' @return A list with `dim`, the selected embedding dimension, and `dE`,
#'   the fraction of false neighbours at each dimension `1:maxdim` (zero
#'   beyond the stopping dimension when `speed = 1`).
#' @references Kennel, M. B., Brown, R. and Abarbanel, H. D. I. (1992).
#'   Determining embedding dimension for phase-space reconstruction using a
#'   geometrical construction. Physical Review A, 45(6), 3403-3411.
#' @export
fnn <- function(x, delay, maxdim, Rtol, Atol, speed) {
  x <- as.numeric(x)
  n <- length(x) - delay * maxdim
  RA <- stats::sd(x)
  z <- x[1:n]
  y <- NULL
  dim <- NULL
  dE <- numeric(maxdim)
  for (j in seq_len(maxdim)) {
    y <- cbind(y, z)
    z <- x[(1 + delay * j):(n + delay * j)]
    L <- logical(n)
    tree <- kd_part(unname(y), z, 512)
    for (i in seq_len(n)) {
      nb <- kd_nearest(tree, y[i, ], 2)
      distance <- nb$pqz[1] - nb$pqz[2]
      if (abs(distance) > nb$pqd[2] * Rtol) L[i] <- TRUE
      if (sqrt(nb$pqd[2]^2 + distance^2) / RA > Atol) L[i] <- TRUE
    }
    dE[j] <- sum(L) / n
    if (speed == 1) {
      if (dE[j] <= 0) { dim <- j; break }
      if (j >= 3 && (dE[j - 2] > dE[j - 1] && dE[j - 1] < dE[j])) { dim <- j - 1; break }
      if (j >= 2 && abs(dE[j] - dE[j - 1]) <= 0.001) { dim <- j; break }
    }
  }
  if (speed == 0 && length(dE) >= 3) {
    for (i in 2:(length(dE) - 1)) {
      if (dE[i] <= 0 || (dE[i - 1] > dE[i] && dE[i] < dE[i + 1])) { dim <- i; break }
      if (dE[i - 1] - dE[i] <= 0.001) { dim <- i - 1; break }
    }
  }
  if (is.null(dim)) {
    dim <- maxdim
    message("no dimension found, dim set to maxdim")
  }
  list(dim = dim, dE = dE)
}

# Build the k-d tree: points are reordered in place so each node owns a
# contiguous block; a node is split at the median of its widest coordinate
# while it holds at least bin_size points.
kd_part <- function(y, z, bin_size) {
  y <- as.matrix(y)
  d <- ncol(y)
  node_list <- matrix(c(1, nrow(y), 0, 0), 1)
  sort_list <- matrix(c(0, 0), 1)
  node <- 1
  last <- 1
  while (node <= last) {
    segment <- node_list[node, 1]:node_list[node, 2]
    yt <- y[segment, , drop = FALSE]
    range <- apply(yt, 2, max) - apply(yt, 2, min)
    if (max(range) > 0 && length(segment) >= bin_size) {
      widest <- order(range)[d]                 # stable: last of any ties
      ord <- order(yt[, widest])
      y_sort <- yt[ord, widest]
      tlen <- nrow(yt)
      cut <- if (tlen %% 2) y_sort[(tlen + 1) / 2] else (y_sort[tlen / 2] + y_sort[tlen / 2 + 1]) / 2
      L <- y_sort <= cut
      if (sum(L) == tlen) {
        L <- y_sort < cut
        cut <- (cut + max(y_sort[L])) / 2
      }
      y[segment, ] <- yt[ord, ]
      z[segment] <- z[segment][ord]
      sort_list[node, ] <- c(widest, cut)
      node_list[node, 3:4] <- c(last + 1, last + 2)
      last <- last + 2
      node_list <- rbind(node_list, c(segment[1], segment[1] + sum(L) - 1, 0, 0),
                         c(segment[1] + sum(L), segment[tlen], 0, 0))
      sort_list <- rbind(sort_list, c(0, 0), c(0, 0))
    }
    node <- node + 1
  }
  list(y = y, z = z, node_list = node_list, sort_list = sort_list)
}

# The m_search nearest points to yq, as the MATLAB kdsearch finds them,
# including its stopping rule: a leaf ends the search as soon as any
# coordinate of yq lies within the current m-th distance of a cell bound.
kd_nearest <- function(tree, yq, m_search) {
  st <- new.env()
  st$pqd <- rep(Inf, m_search)
  st$pqz <- numeric(0)
  st$done <- FALSE
  st$b_upper <- rep(Inf, length(yq))
  st$b_lower <- rep(-Inf, length(yq))
  kd_search(1, tree, yq, m_search, st)
  list(pqd = st$pqd, pqz = st$pqz)
}

kd_search <- function(node, tree, yq, m_search, st) {
  if (st$done) return(invisible())
  nl <- tree$node_list
  if (nl[node, 3] == 0) {
    rows <- nl[node, 1]:nl[node, 2]
    yt <- tree$y[rows, , drop = FALSE]
    dist <- sqrt(rowSums((yt - matrix(yq, nrow(yt), length(yq), byrow = TRUE))^2))
    pqd <- c(dist, st$pqd)
    pqz <- c(tree$z[rows], st$pqz)
    index <- order(pqd)
    pqd <- pqd[index]
    len <- length(pqz)
    pqz <- if (length(index) > len) pqz[index[seq_len(len)]] else pqz[index]
    if (length(pqd) > m_search) pqd <- pqd[1:m_search]
    if (length(pqz) > m_search) pqz <- pqz[1:m_search]
    st$pqd <- pqd
    st$pqz <- pqz
    if (any(abs(yq - st$b_lower) <= pqd[m_search] | abs(yq - st$b_upper) <= pqd[m_search])) {
      st$done <- TRUE
    }
    return(invisible())
  }
  disc <- tree$sort_list[node, 1]
  part <- tree$sort_list[node, 2]
  if (yq[disc] <= part) {
    temp <- st$b_upper[disc]; st$b_upper[disc] <- part
    kd_search(nl[node, 3], tree, yq, m_search, st)
    st$b_upper[disc] <- temp
  } else {
    temp <- st$b_lower[disc]; st$b_lower[disc] <- part
    kd_search(nl[node, 4], tree, yq, m_search, st)
    st$b_lower[disc] <- temp
  }
  if (st$done) return(invisible())
  if (yq[disc] <= part) {
    temp <- st$b_lower[disc]; st$b_lower[disc] <- part
    if (kd_overlap(yq, m_search, st)) kd_search(nl[node, 4], tree, yq, m_search, st)
    st$b_lower[disc] <- temp
  } else {
    temp <- st$b_upper[disc]; st$b_upper[disc] <- part
    if (kd_overlap(yq, m_search, st)) kd_search(nl[node, 3], tree, yq, m_search, st)
    st$b_upper[disc] <- temp
  }
  invisible()
}

# Does the ball of radius pqd[m_search] around yq reach the current cell?
kd_overlap <- function(yq, m_search, st) {
  dist <- st$pqd[m_search]^2
  s <- 0
  for (i in seq_along(yq)) {
    if (yq[i] < st$b_lower[i]) {
      s <- s + (yq[i] - st$b_lower[i])^2
      if (s > dist) return(FALSE)
    } else if (yq[i] > st$b_upper[i]) {
      s <- s + (yq[i] - st$b_upper[i])^2
      if (s > dist) return(FALSE)
    }
  }
  TRUE
}

#' Correlation dimension
#'
#' Grassberger and Procaccia (1983). The correlation sum is computed over
#' 200 log-spaced radii, excluding pairs closer than `dim * delay` samples in
#' time, and the slope is fitted over the band of `ln C` from the midpoint of
#' its range to a quarter-range above it.
#'
#' The region is selected by height, not linearity: check `extra$r2`. A low
#' value means the number is a slope through a bend, not a scaling exponent.
#'
#' @param x Numeric vector.
#' @param delay Embedding delay in samples.
#' @param dim Embedding dimension.
#' @return A list with `CoD`, the correlation dimension, and `extra`: `logEps`,
#'   `logC`, `idx` (the points fitted), `slope`, `intercept`, `r2` and
#'   `bounds`.
#' @references Grassberger, P. and Procaccia, I. (1983). Measuring the
#'   strangeness of strange attractors. Physica D, 9(1-2), 189-208.
#' @export
corr_dim <- function(x, delay, dim) {
  x <- as.numeric(x)
  bins <- 200
  n <- length(x) - (dim - 1) * delay
  y <- embed(x, dim, delay)
  k <- dim * delay
  pair_dist <- function(i) {
    j <- (i + k + 1):n
    sqrt(colSums((y[, j, drop = FALSE] - y[, i])^2))
  }
  eps1 <- 0
  eps2 <- Inf
  for (i in seq_len(n - k - 1)) {
    d <- pair_dist(i)
    eps1 <- max(max(d), eps1)
    eps2 <- min(min(d), eps2)
  }
  if (eps2 == 0) eps2 <- .Machine$double.eps
  lo <- log(eps2)
  hi <- log(eps1)
  e <- lo + (0:(bins - 1)) * (hi - lo) / (bins - 1)
  e[bins] <- hi
  epsilon <- exp(1)^e
  epsilon[c(1, bins)] <- c(eps2, eps1)   # exact end radii; see the MATLAB note
  CI <- numeric(bins)
  for (i in seq_len(n - k - 1)) {
    d <- pair_dist(i)
    # histc: edges[b] <= d < edges[b+1]; the last bin counts d == edges[end].
    b <- findInterval(d, epsilon)
    b[d == epsilon[bins]] <- bins
    b <- b[b >= 1 & (b < bins | d == epsilon[bins])]
    CI <- CI + tabulate(b, bins)
  }
  CI <- log(cumsum(CI) / ((n - k)^2))
  CI <- CI + log(n - k - 1)
  epsilon <- log(epsilon)
  if (is.infinite(min(CI))) {
    i <- max(which(CI == -Inf))
    minCI <- min(CI[(i + 1):bins])
  } else {
    minCI <- min(CI)
  }
  mid <- (max(CI) - minCI) / 2 + minCI
  quarter <- (max(CI) - minCI) / 4
  bounds <- c(mid, mid + quarter)
  intervals <- which(CI > bounds[1] & CI < bounds[2])
  idx <- min(intervals):max(intervals)
  fit <- stats::lm.fit(cbind(epsilon[idx], 1), CI[idx])
  cf <- unname(fit$coefficients)
  yHat <- cf[1] * epsilon[idx] + cf[2]
  ssTot <- sum((CI[idx] - mean(CI[idx]))^2)
  r2 <- if (ssTot > 0 && length(idx) > 2) 1 - sum((CI[idx] - yHat)^2) / ssTot else NA_real_
  list(CoD = cf[1],
       extra = list(logEps = epsilon, logC = CI, idx = idx, slope = cf[1], intercept = cf[2],
                    r2 = r2, bounds = bounds))
}
