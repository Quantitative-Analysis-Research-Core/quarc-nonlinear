# fig_corr_dim_reference_dysts.R
# Correlation dimension per dysts system, as a ratio to the published D2.
# One row per system, ordered by the ratio, coloured by system category — the dysts analogue
# of fig_corr_dim_reference.R. Unlike the Sprott version the axis is logarithmic: with 125
# systems the ratios run from 0.09 to 9.8, and a linear axis stops being readable.
if (!exists(".COMMON_LOADED")) {
  .a <- commandArgs(FALSE); .ff <- sub("--file=", "", grep("--file=", .a, value = TRUE)[1])
  source(file.path(dirname(normalizePath(.ff)), "_common.R"))
}

d <- read_characterization("dysts")
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

FLOOR <- 0.05

p <- ggplot(d, aes(x = ratio, y = systemLabel, colour = categoryLabel)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey55") +
  geom_errorbar(aes(xmin = pmax(ratio - err, FLOOR), xmax = ratio + err),
                orientation = "y", width = 0, linewidth = 0.5, na.rm = TRUE) +
  geom_point(aes(shape = ensemble), size = 1.8, fill = "white", stroke = 0.7) +
  scale_shape_manual(values = c(`Complete Ensemble` = 16, `Incomplete Ensemble` = 21)) +
  scale_colour_manual(values = category_cols, drop = TRUE) +
  scale_x_log10(breaks = c(0.1, 0.3, 1, 3, 10),
                labels = c("0.1", "0.3", "1", "3", "10")) +
  labs(x = "Correlation Dimension ÷ dysts Reference", y = "System") +
  theme_pub +
  theme(axis.text.y = element_text(size = 6.8, margin = margin(r = 5)),
        legend.box = "vertical", legend.box.just = "left",
        legend.margin = margin(b = 0), legend.spacing.y = unit(1, "pt")) +
  guides(colour = guide_legend(nrow = 1, order = 1),
         shape  = guide_legend(nrow = 1, order = 2,
                               override.aes = list(colour = "grey30")))

save_fig(p, "fig_corr_dim_reference_dysts.png", w = 7.6, h = 17.5)
