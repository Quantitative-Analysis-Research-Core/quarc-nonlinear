# fig_scaling_regions.R
# Evidence of linearity: the curves each slope was actually fitted to, with the
# fitted region highlighted and its R^2 stated. Two panels of small multiples,
# one per estimator, saved separately.
#
# A correlation dimension or a Lyapunov exponent is a slope, and a slope means
# nothing unless the stretch it was taken from is straight. The report gives the
# number and its R^2; this shows the shape behind both.
if (!exists(".COMMON_LOADED")) {
  .a <- commandArgs(FALSE); .ff <- sub("--file=", "", grep("--file=", .a, value = TRUE)[1])
  source(file.path(dirname(normalizePath(.ff)), "_common.R"))
}

fit_col  <- "#D55E00"   # the fitted region
rest_col <- "grey65"    # the rest of the curve

# ---- correlation sum curves -------------------------------------------------

cd <- utils::read.csv(res("curves_corr_dim.csv"), stringsAsFactors = FALSE)

# At the smallest radii no pair of points falls inside epsilon, so C(epsilon) is
# zero and ln(C) is -Inf. Those bins carry no information and cannot be drawn;
# they are dropped here rather than left to poison the panel ranges. They are
# never inside a fitted region, so nothing about the fit is affected.
dropped <- sum(!is.finite(cd$logC) | !is.finite(cd$logEps))
if (dropped) cat("  dropped", dropped, "non-finite ln(C) bins\n")
cd <- cd[is.finite(cd$logC) & is.finite(cd$logEps), ]
cd$systemLabel <- nice(cd$system)

# Curvature is computed here from the fitted points themselves, scaled by the
# window's own extent, because R^2 does not measure linearity: a smooth bend
# explains nearly all the variance. On synthetic data a genuine bend scores
# R^2 = 0.99993 and straight-plus-noise scores 0.99993 as well, while their
# curvatures differ sixteenfold.
curv_of <- function(d) {
  f <- d[d$inFit == 1, ]
  if (nrow(f) < 5) return(NA_real_)
  p2 <- stats::coef(stats::lm(logC ~ poly(logEps, 2, raw = TRUE), data = f))
  unname(p2[3] * diff(range(f$logEps))^2 / diff(range(f$logC)))
}

lab_cd <- do.call(rbind, lapply(split(cd, cd$system), function(d) data.frame(
  systemLabel = nice(d$system[1]),
  logEps = min(d$logEps) + 0.02 * diff(range(d$logEps)),
  logC   = max(d$logC)   - 0.06 * diff(range(d$logC)),
  txt    = sprintf("R² = %.4f   curvature = %+.3f\nD₂ = %.2f  (ref %.2f)",
                   d$fitR2[1], curv_of(d), d$slope[1], d$referenceD2[1]),
  stringsAsFactors = FALSE)))

# Ordered by curvature, worst first: that is the quantity that says whether the
# fitted stretch was straight, and ordering by R^2 would put a bent region in
# the middle of the panel.
# na.last keeps a system whose window is too short for a quadratic -- LCG fits
# four points -- in the panel rather than dropping its label to NA.
ord <- names(sort(sapply(split(cd, cd$system), function(d) -abs(curv_of(d))),
                  na.last = TRUE))
cd$systemLabel     <- factor(nice(cd$system),     levels = nice(ord))
lab_cd$systemLabel <- factor(lab_cd$systemLabel, levels = nice(ord))

p_cd <- ggplot(cd, aes(logEps, logC)) +
  geom_point(aes(colour = inFit == 1), size = 0.55) +
  geom_smooth(data = subset(cd, inFit == 1), method = "lm", formula = y ~ x,
              se = FALSE, colour = "black", linewidth = 0.45) +
  geom_text(data = lab_cd, aes(label = txt), hjust = 0, vjust = 1,
            size = 3.0, lineheight = 0.95) +
  scale_colour_manual(values = c(`FALSE` = rest_col, `TRUE` = fit_col),
                      labels = c("Outside Fit", "Fitted Region"),
                      breaks = c(FALSE, TRUE)) +
  facet_wrap(~ systemLabel, ncol = 2, scales = "free") +
  labs(x = "ln(ε)", y = "ln(C(ε))") +
  theme_pub +
  theme(strip.background = element_blank(), strip.text = element_text(size = 14))

save_fig(p_cd, "fig_scaling_region_corr_dim.png", w = 8.2, h = 10.4)

# ---- Rosenstein divergence curves -------------------------------------------

ly <- utils::read.csv(res("curves_lyapunov.csv"), stringsAsFactors = FALSE)
ly <- ly[is.finite(ly$divergence), ]

# The divergence curve is long and its tail is pure saturation; showing all of it
# would compress the fitted rise into the first few pixels. Each panel is cut at
# four times the end of its own fitted window, so the fit and the saturation it
# stops at are both visible.
ly <- do.call(rbind, lapply(split(ly, ly$system), function(d) {
  hi <- max(d$sample[d$inFit == 1], na.rm = TRUE)
  d[d$sample <= max(4 * hi, 40), ]
}))

curv_ly <- function(d) {
  f <- d[d$inFit == 1, ]
  if (nrow(f) < 5) return(NA_real_)
  p2 <- stats::coef(stats::lm(divergence ~ poly(sample, 2, raw = TRUE), data = f))
  unname(p2[3] * diff(range(f$sample))^2 / diff(range(f$divergence)))
}

lab_ly <- do.call(rbind, lapply(split(ly, ly$system), function(d) data.frame(
  systemLabel = nice(d$system[1]),
  sample = min(d$sample) + 0.03 * diff(range(d$sample)),
  divergence = max(d$divergence) - 0.04 * diff(range(d$divergence)),
  txt = sprintf("R² = %.4f   curvature = %s\n%d samples fitted",
                d$fitR2[1],
                ifelse(is.na(curv_ly(d)), "n/a (window too short)",
                       sprintf("%+.3f", curv_ly(d))),
                sum(d$inFit == 1)),
  stringsAsFactors = FALSE)))

ord2 <- names(sort(sapply(split(ly, ly$system), function(d) -abs(curv_ly(d))),
                   na.last = TRUE))
ly$systemLabel     <- factor(nice(ly$system),     levels = nice(ord2))
lab_ly$systemLabel <- factor(lab_ly$systemLabel, levels = nice(ord2))

p_ly <- ggplot(ly, aes(sample, divergence)) +
  geom_line(colour = rest_col, linewidth = 0.4) +
  geom_point(data = subset(ly, inFit == 1), colour = fit_col, size = 0.8) +
  geom_smooth(data = subset(ly, inFit == 1), method = "lm", formula = y ~ x,
              se = FALSE, colour = "black", linewidth = 0.45) +
  geom_text(data = lab_ly, aes(label = txt), hjust = 0, vjust = 1,
            size = 3.0, lineheight = 0.95) +
  facet_wrap(~ systemLabel, ncol = 2, scales = "free") +
  labs(x = "Sample", y = "Mean Log Divergence") +
  theme_pub +
  theme(strip.background = element_blank(), strip.text = element_text(size = 14),
        legend.position = "none")

save_fig(p_ly, "fig_scaling_region_lyapunov.png", w = 8.2, h = 10.4)
