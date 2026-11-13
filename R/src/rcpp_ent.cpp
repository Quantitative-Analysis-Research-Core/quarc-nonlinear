// R wrappers for the C++ core's entropy kernels. See cpp/include/quarc/quarc.hpp.
//
// Copyright (c) 2021-2026 Quantitative Analysis Research Core,
// Center for Human Movement Variability, University of Nebraska at Omaha.
// MIT licence. See LICENSE.

#include <Rcpp.h>
#include "quarc/quarc.hpp"

// [[Rcpp::export(.ent_samp_core)]]
double ent_samp_core(Rcpp::NumericVector x, int dim, double radius) {
    return quarc::ent_samp(x.begin(), x.size(), dim, radius);
}
