# fig_corr_dim_reference.R
# Correlation dimension per system, as a ratio to Sprott's published D2.
# One row per system, ordered by the ratio, coloured by system category.
if (!exists(".COMMON_LOADED")) {
  .a <- commandArgs(FALSE); .ff <- sub("--file=", "", grep("--file=", .a, value = TRUE)[1])
  source(file.path(dirname(normalizePath(.ff)), "_common.R"))
}

d <- read_characterization("full")
d <- d[d$metric == "corr_dim" & is.finite(d$reference) & d$reference != 0, ]

# The point is the ensemble mean and the whisker its 95% CI, both divided by the published
# value so every system lands on one comparable axis. The CI is computed from the stored
# moments rather than raw draws, which is why ci_half is used in place of agg_ci.
d$ratio <- d$mean / d$reference
d$err   <- ci_half(d$sd, d$n) / d$reference

d$systemLabel <- factor(d$systemLabel, levels = d$systemLabel[order(d$ratio)])
d$categoryLabel <- factor(d$categoryLabel, levels = names(category_cols))

# An incomplete ensemble means the surviving realizations are conditioned on not degenerating,
# so those systems are drawn hollow rather than dropped or silently mixed in.
d$complete <- d$n >= d$R

d$ensemble <- factor(ifelse(d$complete, "Complete Ensemble", "Incomplete Ensemble"),
                     levels = c("Complete Ensemble", "Incomplete Ensemble"))

p <- ggplot(d, aes(x = ratio, y = systemLabel, colour = categoryLabel)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey55") +
  geom_errorbar(aes(xmin = ratio - err, xmax = ratio + err),
                orientation = "y", width = 0, linewidth = 0.6, na.rm = TRUE) +
  geom_point(aes(shape = ensemble), size = 2.1, fill = "white", stroke = 0.7) +
  scale_shape_manual(values = c(`Complete Ensemble` = 16, `Incomplete Ensemble` = 21)) +
  scale_colour_manual(values = category_cols, drop = FALSE) +
  scale_x_continuous(breaks = seq(0.4, 2.0, 0.2)) +
  labs(x = "Correlation Dimension ÷ Published Value", y = "System") +
  theme_pub +
  theme(axis.text.y = element_text(size = 9, margin = margin(r = 5)),
        legend.box = "vertical", legend.box.just = "left",
        legend.margin = margin(b = 0), legend.spacing.y = unit(1, "pt")) +
  guides(colour = guide_legend(nrow = 2, byrow = TRUE, order = 1),
         shape  = guide_legend(nrow = 1, order = 2,
                               override.aes = list(colour = "grey30")))

save_fig(p, "fig_corr_dim_reference.png", w = 7.6, h = 10.2)
