# fig_noise_sweep_rqa.R
# RQA metrics against observational noise, under two embedding arms.
#
# The companion to fig_noise_sweep.R for the family that sweep left out. The RQA
# metrics have no published reference, so the quantity drawn is each system's
# value at a noise level divided by its own clean value: how far noise moves the
# measurement from what the same series gave without it. One panel per metric;
# the line is the median of that ratio across systems, the band the interquartile
# range. Colour is the arm, as in the Lyapunov figure: "fixed" holds the clean
# series' delay and dimension, "reembed" re-estimates both on the noisy series.
#
# WHY A RATIO TO THE CLEAN VALUE AND NOT THE VALUE. Determinism runs 60 to 100
# across systems and maximal line length from tens to thousands of samples; drawn
# raw, the between-system spread would swamp the within-system effect of noise,
# which is the thing being asked. Dividing by the clean value puts every system
# at 1 on the left and makes the panels comparable.
#
# The radius is included because it is the one RQA input the sweep re-chooses: it
# is set to hit a target recurrence rate on each series, so noise that inflates
# the radius is a mechanism for the other metrics' changes, not just another
# outcome.
if (!exists(".COMMON_LOADED")) {
  .a <- commandArgs(FALSE); .ff <- sub("--file=", "", grep("--file=", .a, value = TRUE)[1])
  source(file.path(dirname(normalizePath(.ff)), "_common.R"))
}

f <- res("noise_sweep_rqa.csv")
stopifnot(file.exists(f))
d <- utils::read.csv(f, stringsAsFactors = FALSE)

SHOW <- c("rqa_radius", "rqa_det", "rqa_lam", "rqa_maxL_tw", "rqa_entL")
METLBL <- c(rqa_radius = "Radius", rqa_det = "Determinism", rqa_lam = "Laminarity",
            rqa_maxL_tw = "Max Line (Theiler)", rqa_entL = "Line Entropy")

d <- d[d$metric %in% SHOW & is.finite(d$median), ]

# Each system's clean value, taken from the fixed arm at zero noise (the two arms
# are identical there by construction).
clean <- d[d$noise == 0 & d$arm == "fixed", c("system", "metric", "median")]
names(clean)[3] <- "clean"
d <- merge(d, clean, by = c("system", "metric"))
d <- d[is.finite(d$clean) & d$clean != 0, ]
d$ratio <- d$median / d$clean

d$metricLabel <- factor(METLBL[d$metric], levels = unname(METLBL))
d$armLabel <- factor(ifelse(d$arm == "fixed", "Fixed embedding", "Re-embedded"),
                     levels = c("Fixed embedding", "Re-embedded"))
d$noisePct <- 100 * d$noise

qs <- function(x, p) as.numeric(stats::quantile(x, p, na.rm = TRUE))
agg <- do.call(rbind, lapply(
  split(d, list(d$noisePct, d$armLabel, d$metricLabel), drop = TRUE),
  function(g) data.frame(
    noisePct = g$noisePct[1], armLabel = g$armLabel[1], metricLabel = g$metricLabel[1],
    med = median(g$ratio, na.rm = TRUE),
    lo  = qs(g$ratio, 0.25),
    hi  = qs(g$ratio, 0.75),
    n   = sum(is.finite(g$ratio)))))

nsys <- length(unique(d$system))
arm_cols <- c("Fixed embedding" = "#0072B2", "Re-embedded" = "#D55E00")

p <- ggplot(agg, aes(x = noisePct, y = med, colour = armLabel, fill = armLabel)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey55") +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.16, colour = NA) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 1.9) +
  facet_wrap(~ metricLabel, nrow = 1, scales = "free_y") +
  scale_colour_manual(values = arm_cols) +
  scale_fill_manual(values = arm_cols) +
  labs(x = "Observational Noise (% of series SD)", y = "Value ÷ Clean Value",
       caption = sprintf("Line: median across %d systems. Band: interquartile range.", nsys)) +
  theme_pub +
  theme(panel.spacing = unit(1.1, "lines"),
        strip.background = element_blank(),
        strip.text = element_text(size = 14),
        plot.caption = element_text(size = 11, colour = "grey35"))

save_fig(p, "fig_noise_sweep_rqa.png", w = 14, h = 4.8)
