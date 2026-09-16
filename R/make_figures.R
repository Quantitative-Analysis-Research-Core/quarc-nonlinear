# make_figures.R — master script. Rebuilds every figure from the characterization tables.
#
#   Rscript R/make_figures.R
#
# Each fig_*.R also runs standalone. The inputs are tests/reports/characterization_*.csv,
# written by quarctest.write_characterization, and one CSV per sweep (evolve_sweep.csv,
# noise_sweep.csv) written from the table its driver returns; regenerate those with
# quarctest.characterize or the sweep's own driver if the underlying run has changed.
repo <- normalizePath(file.path(dirname(sub(
  "--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])), ".."))

source(file.path(repo, "R", "figures", "_common.R"))

cat("figures ->", file.path(repo, "tests", "reports", "figures"), "\n")
for (f in c("fig_corr_dim_reference.R",
            "fig_lyapunov_reference.R",
            "fig_corr_dim_reference_dysts.R",
            "fig_lyapunov_reference_dysts.R",
            "fig_scaling_regions.R",
            "fig_evolve_sweep.R",
            "fig_noise_sweep.R")) {
  source(file.path(repo, "R", "figures", f))
}
