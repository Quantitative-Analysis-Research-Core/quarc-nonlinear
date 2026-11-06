# fig_lyapunov_reference.R
# Largest Lyapunov exponent per system, as a ratio to Sprott's published value, for both
# estimators. One row per system, ordered by Rosenstein's ratio, coloured by method.
#
# The x axis is logarithmic because the ratios span three orders of magnitude: the linear
# congruential generator, whose exponent is ln 7141, defeats both estimators by a factor of
# roughly a thousand. A linear axis would compress every other system onto the reference line
# and hide the structure that matters.
if (!exists(".COMMON_LOADED")) {
  .a <- commandArgs(FALSE); .ff <- sub("--file=", "", grep("--file=", .a, value = TRUE)[1])
  source(file.path(dirname(normalizePath(.ff)), "_common.R"))
}

d <- read_characterization("full")
d <- d[d$metric %in% c("lyap_wolf", "lyap_ros") & is.finite(d$reference) & d$reference != 0, ]

d$ratio <- d$mean / d$reference
d$err   <- ci_half(d$sd, d$n) / d$reference
d$method <- factor(ifelse(d$metric == "lyap_wolf", "Wolf", "Rosenstein"),
                   levels = c("Wolf", "Rosenstein"))

ord <- d[d$method == "Rosenstein", ]
ord <- ord$systemLabel[order(ord$ratio)]
d$systemLabel <- factor(d$systemLabel, levels = ord)

d <- d[is.finite(d$ratio), ]

# A log axis cannot show a non-positive ratio, and dropping one silently would let an
# estimator that returned the WRONG SIGN on a chaotic system vanish from the figure. Those
# rows are pinned to the left edge with a distinct symbol and labelled instead.
FLOOR <- 0.005
neg <- d[d$ratio <= 0, ]
pos <- d[d$ratio  > 0, ]
if (nrow(neg)) neg$ratio <- FLOOR

method_cols <- c(Wolf = "#D55E00", Rosenstein = "#0072B2")

p <- ggplot(pos, aes(x = ratio, y = systemLabel, colour = method)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey55") +
  geom_errorbar(aes(xmin = pmax(ratio - err, FLOOR), xmax = ratio + err),
                orientation = "y", width = 0, linewidth = 0.6,
                position = position_dodge(width = 0.55), na.rm = TRUE) +
  geom_point(size = 1.9, position = position_dodge(width = 0.55))

if (nrow(neg)) {
  p <- p +
    geom_point(data = neg, shape = 4, size = 2.4, stroke = 1, show.legend = FALSE) +
    geom_text(data = neg, aes(label = "Mean Negative"),
              hjust = -0.14, size = 3.1, show.legend = FALSE)
}

p <- p +
  scale_colour_manual(values = method_cols) +
  scale_x_log10(breaks = c(0.01, 0.1, 1, 10),
                labels = c("0.01", "0.1", "1", "10"),
                limits = c(FLOOR * 0.85, 12)) +
  labs(x = "Lyapunov Exponent ÷ Published Value", y = "System") +
  theme_pub +
  theme(axis.text.y = element_text(size = 9, margin = margin(r = 5)))

save_fig(p, "fig_lyapunov_reference.png", w = 7.6, h = 10.2)
