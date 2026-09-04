# fig_evolve_sweep.R
# Wolf's exponent against its renormalisation interval, one line per system.
#
# EVOLVE is how long lye_w lets a pair separate before replacing the neighbour.
# Wolf's method assumes the pair stays in the linear regime over that interval,
# and lyapunov.m applies one default of 10 samples to every system. On a flow
# decimated to ~40 samples per orbit that is a quarter turn; on a map it is 10
# iterations, over which a pair on the logistic map separates by 2^10.
if (!exists(".COMMON_LOADED")) {
  .a <- commandArgs(FALSE); .ff <- sub("--file=", "", grep("--file=", .a, value = TRUE)[1])
  source(file.path(dirname(normalizePath(.ff)), "_common.R"))
}

d <- utils::read.csv(res("evolve_sweep.csv"), stringsAsFactors = FALSE)
d <- d[is.finite(d$ratio) & d$ratio > 0, ]
d$systemLabel <- nice(d$system)
d$kindLabel <- factor(ifelse(d$kind == "map", "Maps", "Flows"),
                      levels = c("Maps", "Flows"))

kind_cols <- c(Maps = "#D55E00", Flows = "#0072B2")

p <- ggplot(d, aes(evolve, ratio, group = system, colour = kindLabel)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey55") +
  geom_line(linewidth = 0.45, alpha = 0.75) +
  geom_point(size = 1.1, alpha = 0.85) +
  scale_colour_manual(values = kind_cols) +
  scale_x_log10(breaks = c(1, 2, 3, 5, 10, 20, 30)) +
  scale_y_log10(breaks = c(0.05, 0.1, 0.25, 0.5, 1, 2),
                labels = c("0.05", "0.1", "0.25", "0.5", "1", "2")) +
  facet_wrap(~ kindLabel, ncol = 2) +
  labs(x = "Renormalisation Interval (Samples)",
       y = "Wolf Exponent ÷ Published Value") +
  theme_pub +
  theme(strip.background = element_blank(), strip.text = element_text(size = 14),
        legend.position = "none")

save_fig(p, "fig_evolve_sweep.png", w = 8.4, h = 5.0)

# ---- the same thing as one number per system: where the ratio is closest to 1

best <- do.call(rbind, lapply(split(d, d$system), function(s) {
  i <- which.min(abs(log(s$ratio)))
  data.frame(system = s$system[1], kindLabel = s$kindLabel[1],
             bestEvolve = s$evolve[i], bestRatio = s$ratio[i],
             defaultRatio = s$ratio[s$evolve == 10][1],
             stringsAsFactors = FALSE)
}))
best$systemLabel <- factor(nice(best$system),
                           levels = nice(best$system[order(best$bestEvolve,
                                                           best$system)]))

p2 <- ggplot(best, aes(bestEvolve, systemLabel, colour = kindLabel)) +
  geom_point(size = 2.2) +
  geom_vline(xintercept = 10, linetype = "dashed", colour = "grey55") +
  scale_colour_manual(values = kind_cols) +
  scale_x_log10(breaks = c(1, 2, 3, 5, 10, 20, 30)) +
  labs(x = "Renormalisation Interval Closest to the Published Value",
       y = "System") +
  theme_pub +
  theme(axis.text.y = element_text(size = 9, margin = margin(r = 5)))

save_fig(p2, "fig_evolve_best.png", w = 7.4, h = 8.2)
