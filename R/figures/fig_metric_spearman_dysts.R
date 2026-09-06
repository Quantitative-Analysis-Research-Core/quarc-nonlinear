# fig_metric_spearman_dysts.R
# Spearman correlation between the battery's headline metrics, twice: within systems
# (each system's metrics correlated across its ~100 realizations, then the elementwise
# median of those matrices) and between systems (each system collapsed to its median
# metric vector, correlated across the catalogue). The two panels answer different
# questions -- do metrics co-fluctuate under IC noise, and do they rank systems the
# same way -- and the contrast is the point: between-system structure is strong
# (Wolf and Rosenstein agree at rho 0.84) while within-system correlations are
# mostly noise, except the estimator pair (0.47) and a small embedding-into-lambda
# coupling that the per-realization embedding protocol makes visible.
#
# INPUT: tests/reports/characterization_dysts_per_realization.csv, which is
# deliberately untracked (26 MB of per-realization draws). Regenerate it with the
# characterize_series run documented in tests/README.md; the committed summary csv
# cannot substitute, because the within-system panel needs the raw draws.
if (!exists(".COMMON_LOADED")) {
  .a <- commandArgs(FALSE); .ff <- sub("--file=", "", grep("--file=", .a, value = TRUE)[1])
  source(file.path(dirname(normalizePath(.ff)), "_common.R"))
}

per_csv <- res("characterization_dysts_per_realization.csv")
if (!file.exists(per_csv)) {
  cat("  fig_metric_spearman_dysts.png SKIPPED:", basename(per_csv), "not present\n")
} else {

d <- utils::read.csv(per_csv, stringsAsFactors = FALSE)

METRICS <- c(lyap_ros = "Rosenstein λ", lyap_wolf = "Wolf λ",
             corr_dim = "D2", ami_delay = "AMI delay", fnn_dim = "FNN dim",
             dfa_alpha = "DFA α", ent_samp = "SampEn", ent_permu = "PermEn",
             rqa_det = "RQA %DET", rqa_meanL = "RQA mean L")
k <- names(METRICS)

d <- d[d$metric %in% k & d$seriesUsable == 1, ]

wide_of <- function(s) {
  w <- reshape(s[, c("realization", "metric", "value")],
               idvar = "realization", timevar = "metric", direction = "wide")
  names(w) <- sub("^value\\.", "", names(w))
  w[, intersect(k, names(w)), drop = FALSE]
}

systems <- split(d, d$system)
mats <- Filter(Negate(is.null), lapply(systems, function(s) {
  w <- wide_of(s)
  if (nrow(w) < 10) return(NULL)
  suppressWarnings(cor(w, method = "spearman", use = "pairwise.complete.obs"))
}))

aligned <- lapply(mats, function(m) {
  out <- matrix(NA_real_, length(k), length(k), dimnames = list(k, k))
  out[rownames(m), colnames(m)] <- m
  out
})
within_med <- apply(simplify2array(aligned), c(1, 2), median, na.rm = TRUE)

medvecs <- t(sapply(systems, function(s)
  sapply(k, function(m) median(s$value[s$metric == m], na.rm = TRUE))))
between <- suppressWarnings(cor(medvecs, method = "spearman",
                                use = "pairwise.complete.obs"))

melt_mat <- function(m, panel) {
  g <- expand.grid(x = k, y = k, stringsAsFactors = FALSE)
  g$rho <- mapply(function(a, b) m[a, b], g$x, g$y)
  g$rho[g$x == g$y] <- NA               # the diagonal carries no information
  g$panel <- panel
  g
}
long <- rbind(
  melt_mat(within_med, sprintf("Within systems (median of %d matrices)", length(aligned))),
  melt_mat(between,    sprintf("Between systems (%d system medians)", nrow(medvecs))))
long$x <- factor(long$x, levels = k, labels = METRICS[k])
long$y <- factor(long$y, levels = rev(k), labels = METRICS[rev(k)])
long$panel <- factor(long$panel, levels = unique(long$panel))
long$lbl <- ifelse(is.na(long$rho), "", sprintf("%.2f", long$rho))

# Correlation is a polarity, so the fill is diverging: the two method colours the
# other figures already use as poles, a neutral near-white at zero.
p <- ggplot(long, aes(x, y, fill = rho)) +
  geom_tile(colour = "white", linewidth = 1.2) +
  geom_text(aes(label = lbl), size = 2.9, colour = "grey15") +
  facet_wrap(~panel) +
  scale_fill_gradient2(low = "#0072B2", mid = "grey96", high = "#D55E00",
                       limits = c(-1, 1), na.value = "grey88",
                       name = "Spearman ρ") +
  coord_fixed() +
  labs(x = NULL, y = NULL) +
  theme_pub +
  theme(axis.text = element_text(size = 11),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.line = element_blank(), axis.ticks = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 12, face = "bold"),
        legend.position = "right", legend.title = element_text(size = 12))

save_fig(p, "fig_metric_spearman_dysts.png", w = 12.5, h = 6.2)
}
