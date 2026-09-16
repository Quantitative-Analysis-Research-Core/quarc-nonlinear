# fig_noise_sweep.R
# Estimator accuracy against observational noise, under two embedding arms.
#
# One panel per referenced metric. The line is the median across systems, the band the
# interquartile range across systems. Colour is the arm: "fixed" holds the clean series'
# delay and dimension, "reembed" re-estimates both on the noisy series the way a
# practitioner would. The gap between the two lines is the part of the damage that belongs
# to the embedding step rather than to the estimator.
#
# WHY A BAND AND NOT ONE LINE PER SYSTEM. The first draft drew all 126 systems faintly
# behind the medians. At 10% noise a few systems reach ratios of 100 and -150 -- Wolf
# overshooting on near-zero exponents, Rosenstein returning the wrong sign -- and they set
# the axis, flattening every median into an unreadable smear at the bottom. The quartile
# band shows the same spread where the reader can see it, and the tails are reported in
# the sweep's table rather than drawn.
#
# The y axis is linear, not logarithmic as in fig_lyapunov_reference_dysts.R: that figure
# spans orders of magnitude across systems, while this one asks how a ratio near 1 decays
# with noise, and a log axis would flatten exactly the part being read.
if (!exists(".COMMON_LOADED")) {
  .a <- commandArgs(FALSE); .ff <- sub("--file=", "", grep("--file=", .a, value = TRUE)[1])
  source(file.path(dirname(normalizePath(.ff)), "_common.R"))
}

f <- res("noise_sweep.csv")
stopifnot(file.exists(f))
d <- utils::read.csv(f, stringsAsFactors = FALSE)

REFERENCED <- c("lyap_wolf", "lyap_ros", "corr_dim")
METLBL <- c(lyap_wolf = "Wolf", lyap_ros = "Rosenstein", corr_dim = "Correlation Dimension")

d <- d[d$metric %in% REFERENCED & is.finite(d$medianOverReference), ]
d$metricLabel <- factor(METLBL[d$metric], levels = unname(METLBL))
d$armLabel <- factor(ifelse(d$arm == "fixed", "Fixed embedding", "Re-embedded"),
                     levels = c("Fixed embedding", "Re-embedded"))
d$noisePct <- 100 * d$noise

qs <- function(x, p) as.numeric(stats::quantile(x, p, na.rm = TRUE))
agg <- do.call(rbind, lapply(
  split(d, list(d$noisePct, d$armLabel, d$metricLabel), drop = TRUE),
  function(g) data.frame(
    noisePct = g$noisePct[1], armLabel = g$armLabel[1], metricLabel = g$metricLabel[1],
    med = median(g$medianOverReference, na.rm = TRUE),
    lo  = qs(g$medianOverReference, 0.25),
    hi  = qs(g$medianOverReference, 0.75),
    n   = sum(is.finite(g$medianOverReference)))))

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
  labs(x = "Observational Noise (% of series SD)", y = "Estimate ÷ dysts Reference",
       caption = sprintf("Line: median across %d systems. Band: interquartile range.", nsys)) +
  theme_pub +
  theme(panel.spacing = unit(1.1, "lines"),
        strip.background = element_blank(),
        strip.text = element_text(size = 14),
        plot.caption = element_text(size = 11, colour = "grey35"))

save_fig(p, "fig_noise_sweep.png", w = 11, h = 4.8)
