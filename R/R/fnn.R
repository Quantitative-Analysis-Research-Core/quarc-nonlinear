#' False nearest neighbours
#'
#' Kennel, Brown and Abarbanel (1992). For each dimension up to `maxdim`, the
#' fraction of nearest neighbours that separate when one more delay
#' coordinate is added. Neighbours are found exactly with a k-d tree of
#' 512-point leaves, reproduced from the MATLAB implementation step for step.
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
  dim <- NULL
  dE <- numeric(maxdim)
  for (j in seq_len(maxdim)) {
    dE[j] <- if (n >= 2) .fnn_core(x, delay, j, n, Rtol, Atol, RA)
             else fnn_false_fraction_r(x, delay, j, n, Rtol, Atol, RA)
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

# The R code for one dimension's false-neighbour fraction, kept as the
# readable reference for the compiled kernel (cpp/src/fnn.cpp), which ports
# the same k-d tree step for step; the tests check that the two agree.
fnn_false_fraction_r <- function(x, delay, j, n, Rtol, Atol, RA) {
  y <- sapply(0:(j - 1), function(c) x[(1:n) + delay * c])
  y <- matrix(y, n, j)
  z <- x[(1:n) + delay * j]
  L <- logical(n)
  tree <- kd_part(y, z, 512)
  for (i in seq_len(n)) {
    nb <- kd_nearest(tree, y[i, ], 2)
    distance <- nb$pqz[1] - nb$pqz[2]
    if (abs(distance) > nb$pqd[2] * Rtol) L[i] <- TRUE
    if (sqrt(nb$pqd[2]^2 + distance^2) / RA > Atol) L[i] <- TRUE
  }
  sum(L) / n
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

# The m_search nearest points to yq: descend to yq's leaf, then visit
# neighbouring cells only while the search ball can reach them.
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
    # Ball-within-bounds: stop only when the ball lies wholly inside the cell.
    if (all(abs(yq - st$b_lower) > pqd[m_search] & abs(yq - st$b_upper) > pqd[m_search])) {
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
