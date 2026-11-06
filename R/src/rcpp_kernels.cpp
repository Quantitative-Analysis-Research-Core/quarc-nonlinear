// R wrappers for the C++ core's line_hist, corr_dim, fnn and lye_r kernels.
// See cpp/include/quarc/quarc.hpp.
//
// Copyright (c) 2021-2026 Quantitative Analysis Research Core,
// Center for Human Movement Variability, University of Nebraska at Omaha.
// MIT licence. See LICENSE.

#include <Rcpp.h>
#include "quarc/quarc.hpp"

#include <cstdint>
#include <vector>

// Diagonal and vertical line lengths of a square 0/1 matrix.
// [[Rcpp::export(.line_hist_core)]]
Rcpp::List line_hist_core(Rcpp::NumericMatrix a, double theiler, bool cross) {
    const std::size_t n = a.nrow();
    if (static_cast<std::size_t>(a.ncol()) != n) Rcpp::stop("a must be square.");
    std::vector<std::uint8_t> m(n * n);
    for (std::size_t k = 0; k < n * n; ++k) m[k] = a[k] != 0;
    auto h = quarc::line_hist(m.data(), n, theiler > 0 ? static_cast<std::size_t>(theiler) : 0, cross);
    return Rcpp::List::create(Rcpp::Named("diag") = Rcpp::wrap(h.diag),
                              Rcpp::Named("vert") = Rcpp::wrap(h.vert));
}

// Smallest and largest pair distance, c(min, max).
// [[Rcpp::export(.corr_dim_range_core)]]
Rcpp::NumericVector corr_dim_range_core(Rcpp::NumericVector x, int dim, double delay) {
    double lo = 0, hi = 0;
    quarc::corr_dim_range(x.begin(), x.size(), dim, static_cast<std::size_t>(delay), lo, hi);
    return Rcpp::NumericVector::create(lo, hi);
}

// histc counts of all pair distances against edges.
// [[Rcpp::export(.corr_dim_counts_core)]]
Rcpp::NumericVector corr_dim_counts_core(Rcpp::NumericVector x, int dim, double delay,
                                         Rcpp::NumericVector edges) {
    return Rcpp::wrap(quarc::corr_dim_counts(x.begin(), x.size(), dim,
                                             static_cast<std::size_t>(delay),
                                             edges.begin(), edges.size()));
}

// Fraction of false nearest neighbours when dimension j + 1 is added.
// [[Rcpp::export(.fnn_core)]]
double fnn_core(Rcpp::NumericVector x, double delay, double j, double n,
                double rtol, double atol, double ra) {
    return quarc::fnn_false_fraction(x.begin(), x.size(), static_cast<std::size_t>(delay),
                                     static_cast<std::size_t>(j), static_cast<std::size_t>(n),
                                     rtol, atol, ra);
}

// Rosenstein's nearest neighbours (1-based) and mean log divergence.
// [[Rcpp::export(.lye_r_core)]]
Rcpp::List lye_r_core(Rcpp::NumericMatrix Y, double theiler) {
    auto r = quarc::lye_r_divergence(Y.begin(), Y.nrow(), Y.ncol(),
                                     theiler > 0 ? static_cast<std::size_t>(theiler) : 0);
    for (auto& v : r.nn) v += 1;
    return Rcpp::List::create(Rcpp::Named("nn") = Rcpp::wrap(r.nn),
                              Rcpp::Named("avg") = Rcpp::wrap(r.avg),
                              Rcpp::Named("count") = static_cast<double>(r.count));
}
