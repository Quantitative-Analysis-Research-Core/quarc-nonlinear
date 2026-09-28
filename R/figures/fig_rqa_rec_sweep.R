# fig_rqa_rec_sweep.R
# What the RQA recurrence target does to the battery, across 2.5 / 5 / 10 percent on the
# same stored dysts ensembles. Two panels, two findings.
#
# Left: the determinism ceiling, drawn as distance from 100 on a log axis because the
# interesting variation lives in the last fraction of a percent. Raising the target does
# not unsaturate DET -- it saturates it harder, because a looser radius admits recurrences
# that overwhelmingly join existing diagonals. The knob a user would reach for goes the
# wrong way; the ceiling belongs to dmin (every diagonal of length 2 counts, and at 40
# samples per period that is 1/20th of an orbit), which the min-line-length experiments
# take up separately.
#
# Right: mean diagonal line length is the threshold-robust member of the family -- the
# system ordering barely moves between the tightest and loosest target, so conclusions
# built on meanL do not hinge on the recurrence-rate choice.
if (!exists(".COMMON_LOADED")) {
  .a <- commandArgs(FALSE); .ff <- sub("--file=", "", grep("--file=", .a, value = TRUE)[1])
  source(file.path(dirname(normalizePath(.ff)), "_common.R"))
}

d <- utils::read.csv(res("rqa_rec_sweep.csv"), stringsAsFactors = FALSE)
d$categoryLabel <- factor(nice_cat(d$category), levels = names(category_cols))

det <- d[d$metric == "rqa_det" & is.finite(d$median), ]
det$gap <- pmax(100 - det$median, 1e-3)
det$rt <- factor(det$recTarget, levels = c(2.5, 5, 10),
                 labels = c("2.5%", "5%", "10%"))

pA <- ggplot(det, aes(rt, gap, group = system)) +
  geom_line(colour = "grey70", linewidth = 0.3, alpha = 0.6) +
  geom_point(aes(colour = categoryLabel), size = 1.3) +
  stat_summary(aes(group = 1), fun = median, geom = "line",
               colour = "black", linewidth = 1.1) +
  scale_y_continuous(trans = scales::compose_trans("log10", "reverse"),
                     breaks = c(0.01, 0.1, 1),
                     labels = c("99.99", "99.9", "99")) +
  scale_colour_manual(values = category_cols, drop = TRUE) +
  guides(colour = guide_legend(override.aes = list(size = 2))) +
  labs(x = "Recurrence target", y = "Determinism (%)") +
  theme_pub

ml <- d[d$metric == "rqa_meanL" & is.finite(d$median), ]
w <- merge(ml[ml$recTarget == 2.5, c("system", "categoryLabel", "median")],
           ml[ml$recTarget == 10,  c("system", "median")], by = "system",
           suffixes = c("_lo", "_hi"))
rho <- cor(w$median_lo, w$median_hi, method = "spearman")

pB <- ggplot(w, aes(median_lo, median_hi, colour = categoryLabel)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey55") +
  geom_point(size = 1.6) +
  scale_x_log10() + scale_y_log10() +
  scale_colour_manual(values = category_cols, drop = TRUE) +
  guides(colour = guide_legend(override.aes = list(size = 2))) +
  annotate("text", x = min(w$median_lo) * 1.15, y = max(w$median_hi) * 0.85,
           hjust = 0, size = 4.2,
           label = sprintf("Spearman ρ = %.3f", rho)) +
  labs(x = "Mean line length at 2.5%", y = "Mean line length at 10%") +
  theme_pub

if (requireNamespace("patchwork", quietly = TRUE)) {
  library(patchwork)
  p <- pA + pB + plot_layout(guides = "collect") & theme(legend.position = "top")
  save_fig(p, "fig_rqa_rec_sweep.png", w = 11.5, h = 5.4)
} else {
  save_fig(pA, "fig_rqa_rec_sweep_det.png",   w = 6.4, h = 5.4)
  save_fig(pB, "fig_rqa_rec_sweep_meanL.png", w = 6.4, h = 5.4)
}
