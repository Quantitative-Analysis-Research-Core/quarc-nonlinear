# fig_rqa_dmin_sweep.R
# What the minimum line length does to determinism, on line histograms counted exactly.
#
# Left: every system's median DET against dmin. At the default of 2 -- 1/20th of an orbit
# at this sampling -- the catalogue is pressed into the low-90s-and-up band; by dmin=40,
# one full orbital period, it spans nearly the whole axis. The recurrence-target sweep
# showed the radius saturates DET harder; this is the knob that actually spreads it.
#
# Right: what the spread buys. Spearman correlation between DET and the reference largest
# Lyapunov exponent PER SAMPLE (lambda/fs; line lengths are measured in samples, and the
# export fixes sampling at 40 per period, so per-unit-time lambda would smuggle every
# system's natural timescale into the comparison), per dmin. MeanL is drawn alongside.
if (!exists(".COMMON_LOADED")) {
  .a <- commandArgs(FALSE); .ff <- sub("--file=", "", grep("--file=", .a, value = TRUE)[1])
  source(file.path(dirname(normalizePath(.ff)), "_common.R"))
}

d <- utils::read.csv(res("rqa_dmin_sweep.csv"), stringsAsFactors = FALSE)
d <- d[d$n > 0, ]
d$categoryLabel <- factor(nice_cat(d$category), levels = names(category_cols))

ref <- utils::read.csv(res("characterization_dysts.csv"), stringsAsFactors = FALSE)
ref <- ref[ref$metric == "lyap_wolf", c("system", "reference", "fs")]
ref$reference <- ref$reference / ref$fs      # nats per sample

pA <- ggplot(d, aes(dmin, det_median, group = system)) +
  geom_line(colour = "grey70", linewidth = 0.3, alpha = 0.6) +
  geom_point(aes(colour = categoryLabel), size = 1.2) +
  stat_summary(aes(group = 1), fun = median, geom = "line",
               colour = "black", linewidth = 1.1) +
  scale_x_log10(breaks = c(2, 5, 10, 20, 40)) +
  scale_colour_manual(values = category_cols, drop = TRUE) +
  guides(colour = guide_legend(override.aes = list(size = 2))) +
  labs(x = "Minimum line length (samples)", y = "Determinism (%)") +
  theme_pub

m <- merge(d, ref, by = "system")
m <- m[is.finite(m$reference) & m$reference > 0, ]
rows <- do.call(rbind, lapply(split(m, m$dmin), function(g) {
  data.frame(dmin = g$dmin[1],
             DET     = cor(g$det_median,   g$reference, method = "spearman"),
             `MeanL` = cor(g$meanL_median, g$reference, method = "spearman"),
             check.names = FALSE)
}))
long <- reshape(rows, direction = "long", varying = c("DET", "MeanL"),
                v.names = "rho", timevar = "metric", times = c("DET", "MeanL"))
long$metric <- factor(long$metric, levels = c("DET", "MeanL"))

metric_cols <- c(DET = "#D55E00", MeanL = "#0072B2")

pB <- ggplot(long, aes(dmin, rho, colour = metric)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey55") +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.2) +
  scale_x_log10(breaks = c(2, 5, 10, 20, 40)) +
  scale_colour_manual(values = metric_cols) +
  labs(x = "Minimum line length (samples)",
       y = "Spearman ρ with reference λ per sample") +
  theme_pub

if (requireNamespace("patchwork", quietly = TRUE)) {
  library(patchwork)
  p <- pA + pB & theme(legend.position = "top")
  save_fig(p, "fig_rqa_dmin_sweep.png", w = 11.5, h = 5.4)
} else {
  save_fig(pA, "fig_rqa_dmin_sweep_det.png", w = 6.4, h = 5.4)
  save_fig(pB, "fig_rqa_dmin_sweep_rho.png", w = 6.4, h = 5.4)
}
