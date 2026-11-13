rqa_finish <- function(x, a, wrp, param, threshold, type, iter, theiler, dmin, vmin,
                       RESULTS, withEntrW) {
  if (param == "rad") {
    radius <- threshold
    lh <- line_hist(x, a, threshold, type, theiler)
    recurrence <- lh$recurrence
    dh <- lh$diag_hist
    vh <- lh$vertical_hist
    A <- lh$a
  } else if (param == "rec") {
    sr <- set_radius(x, a, 0.01, 0.5, threshold, type, iter, theiler)
    recurrence <- sr$rec
    dh <- sr$diag_hist
    vh <- sr$vertical_hist
    radius <- sr$rad_final
    A <- sr$A
  } else {
    stop('param must be "rad" or "rec".', call. = FALSE)
  }
  RESULTS$RADIUS <- radius
  RESULTS$Size <- nrow(A)
  RESULTS$REC <- recurrence
  line_entropy <- function(h) {
    p <- tabulate(h)
    p <- p[p > 0] / sum(p)
    -sum(p * log2(p))
  }
  if (recurrence > 0) {
    d <- dh[dh >= dmin]
    v <- vh[vh >= vmin]
    RESULTS$DET <- 100 * sum(d) / sum(dh)
    RESULTS$MeanL <- mean(d)
    RESULTS$MaxL <- max(d)
    RESULTS$EntrL <- line_entropy(d)
    RESULTS$LAM <- 100 * sum(v) / sum(vh)
    RESULTS$MeanV <- mean(v)
    RESULTS$MaxV <- max(v)
    RESULTS$EntrV <- line_entropy(v)
    if (withEntrW) RESULTS$EntrW <- ent_weighted(wrp)
  } else {
    for (f in c("DET", "MeanL", "MaxL", "EntrL", "LAM", "MeanV", "MaxV", "EntrV")) RESULTS[[f]] <- NaN
    if (withEntrW) RESULTS$EntrW <- NaN
  }
  M <- 1 - A
  list(RP = t(M)[ncol(M):1, , drop = FALSE], RESULTS = RESULTS)
}
