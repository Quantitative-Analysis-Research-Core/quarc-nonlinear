# fig_rqa_vmin_sweep.R
# Laminarity against its minimum line length. Where the dmin sweep found a knob, this
# finds a cliff: median LAM is 95.3 at vmin=2, 41.2 at 5, and effectively zero from 10
# on. Laminar episodes on these flows do not outlive a quarter period at 40 samples per
# period, and no member of the vertical family correlates with the exponent at any vmin
# (|rho| <= 0.24) -- laminarity measures sojourn, not divergence. One panel, because the
# finding is one sentence; the table behind it is rqa_vmin_sweep.csv.
if (!exists(".COMMON_LOADED")) {
  .a <- commandArgs(FALSE); .ff <- sub("--file=", "", grep("--file=", .a, value = TRUE)[1])
  source(file.path(dirname(normalizePath(.ff)), "_common.R"))
}

d <- utils::read.csv(res("rqa_vmin_sweep.csv"), stringsAsFactors = FALSE)
d <- d[is.finite(d$median_lam), ]
d$categoryLabel <- factor(nice_cat(d$category), levels = names(category_cols))

p <- ggplot(d, aes(vmin, median_lam, group = system)) +
  geom_line(colour = "grey70", linewidth = 0.3, alpha = 0.6) +
  geom_point(aes(colour = categoryLabel), size = 1.2) +
  stat_summary(aes(group = 1), fun = median, geom = "line",
               colour = "black", linewidth = 1.1) +
  scale_x_log10(breaks = c(2, 5, 10, 20, 40)) +
  scale_colour_manual(values = category_cols, drop = TRUE) +
  guides(colour = guide_legend(override.aes = list(size = 2))) +
  labs(x = "Minimum vertical line length (samples)", y = "Laminarity (%)") +
  theme_pub

save_fig(p, "fig_rqa_vmin_sweep.png", w = 6.8, h = 5.4)
