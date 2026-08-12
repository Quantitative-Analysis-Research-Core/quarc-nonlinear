# _common.R — shared setup for per-figure scripts in R/figures/ (house style).
# Sourced by the master R/make_figures.R before each fig_*.R, and auto-sourced when a figure
# script is run standalone.
#
# Two deviations from the standard layout, both structural. There is no results/ directory:
# the tidy CSV interface is tests/reports/characterization_full.csv, written by
# quarctest.write_characterization. And there is no paper/ directory, because this is a
# library rather than a manuscript, so figures land in tests/reports/figures/.
suppressPackageStartupMessages({ library(ggplot2) })

if (!exists("repo")) {
  .args <- commandArgs(FALSE)
  .f <- sub("--file=", "", grep("--file=", .args, value = TRUE)[1])
  .here <- if (length(.f) && !is.na(.f)) dirname(normalizePath(.f)) else normalizePath("R/figures")
  repo <- normalizePath(file.path(.here, "..", ".."))
}
res <- function(f) file.path(repo, "tests", "reports", f)
fig <- function(f) {
  d <- file.path(repo, "tests", "reports", "figures")
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
  file.path(d, f)
}
has_repel <- requireNamespace("ggrepel", quietly = TRUE)

## ---- shared style ---------------------------------------------------------
theme_pub <- theme_classic(base_size = 15) +
  theme(axis.line = element_line(linewidth = 0.4),
        axis.ticks = element_line(linewidth = 0.4),
        axis.ticks.length = unit(-1.5, "mm"),                  # ticks point INWARD
        axis.text = element_text(size = 14, color = "black"),
        axis.text.x = element_text(margin = margin(t = 5)),    # keep labels off inward ticks
        axis.text.y = element_text(margin = margin(r = 5)),
        axis.title = element_text(size = 16),
        legend.position = "top", legend.title = element_blank(),
        legend.text = element_text(size = 13),
        legend.key.size = unit(0.9, "lines"),
        plot.margin = margin(8, 12, 8, 8))

# Six nominal categories, so a qualitative palette rather than an ordered one.
# Okabe-Ito, minus the yellow, which does not hold up on white.
category_cols <- c(
  "Noninvertible Map" = "#E69F00",
  "Dissipative Map"   = "#D55E00",
  "Conservative Map"  = "#CC79A7",
  "Driven Flow"       = "#009E73",
  "Autonomous Flow"   = "#0072B2",
  "Conservative Flow" = "#56B4E9")

CATLBL <- c(noninvertible_map = "Noninvertible Map",
            dissipative_map   = "Dissipative Map",
            conservative_map  = "Conservative Map",
            driven_flow       = "Driven Flow",
            autonomous_flow   = "Autonomous Flow",
            conservative_flow = "Conservative Flow")

# System ids that a generic title-caser would get wrong.
LBL <- c(lcg = "LCG", act = "ACT", windmi = "WINDMI",
         driven_vdp = "Driven VdP", shaw_vdp = "Shaw VdP", duffing_vdp = "Duffing-VdP",
         henon = "Henon", henon_area = "Henon Area-Preserving", henon_heiles = "Henon-Heiles",
         lorenz3d_map = "Lorenz 3D Map", diffusionless_lorenz = "Diffusionless Lorenz",
         rabinovich_fabrikant = "Rabinovich-Fabrikant", nose_hoover = "Nose-Hoover",
         burke_shaw = "Burke-Shaw", moore_spiegel = "Moore-Spiegel",
         duffing_two_well = "Duffing Two-Well", rayleigh_duffing = "Rayleigh-Duffing",
         delayed_logistic = "Delayed Logistic", dissipative_standard = "Dissipative Standard",
         predator_prey = "Predator-Prey", chaotic_web = "Chaotic Web",
         gingerbreadman = "Gingerbreadman", sine_circle = "Sine-Circle",
         simplest_quadratic = "Simplest Quadratic", simplest_cubic = "Simplest Cubic",
         simplest_piecewise = "Simplest Piecewise", simplest_driven = "Simplest Driven",
         double_scroll = "Double Scroll", complex_butterfly = "Complex Butterfly",
         damped_pendulum = "Damped Pendulum", driven_pendulum = "Driven Pendulum")

.title_case <- function(x) {
  x <- gsub("_", " ", x)
  gsub("(^|\\s)(\\w)", "\\1\\U\\2", x, perl = TRUE)
}
nice <- function(x) ifelse(x %in% names(LBL), LBL[x], .title_case(x))
nice_cat <- function(x) ifelse(x %in% names(CATLBL), CATLBL[x], .title_case(x))

# 95% CI half-width (t-based) from a summary row that already carries n, mean and sd.
# The characterization tables are aggregates, so the per-observation form of agg_ci does not
# apply; this is the same quantity computed from the moments that were stored.
ci_half <- function(sd, n, conf = 0.95) {
  ifelse(is.finite(sd) & n > 1, qt(1 - (1 - conf) / 2, n - 1) * sd / sqrt(n), NA_real_)
}

read_characterization <- function(tag = "full") {
  f <- res(paste0("characterization_", tag, ".csv"))
  stopifnot(file.exists(f))
  d <- utils::read.csv(f, stringsAsFactors = FALSE)
  d$categoryLabel <- nice_cat(d$category)
  d$systemLabel   <- nice(d$system)
  d
}

save_fig <- function(p, f, w = 7, h = 4.6) {
  ggsave(fig(f), p, width = w, height = h, dpi = 300); cat("  ", f, "\n")
}

.COMMON_LOADED <- TRUE
