# make_figures.R — master script. Rebuilds every figure from the characterization tables.
#
#   Rscript R/make_figures.R
#
# Each fig_*.R also runs standalone. The input is tests/reports/characterization_*.csv,
# written by quarctest.write_characterization; regenerate those with quarctest.characterize
# before rebuilding if the underlying run has changed.
repo <- normalizePath(file.path(dirname(sub(
  "--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])), ".."))

source(file.path(repo, "R", "figures", "_common.R"))

cat("figures ->", file.path(repo, "tests", "reports", "figures"), "\n")
for (f in c("fig_corr_dim_reference.R",
            "fig_lyapunov_reference.R",
            "fig_corr_dim_reference_dysts.R",
            "fig_lyapunov_reference_dysts.R",
            "fig_metric_spearman_dysts.R",
            "fig_scaling_regions.R",
            "fig_evolve_sweep.R")) {
  source(file.path(repo, "R", "figures", f))
}
