// Sample entropy. See quarc.hpp and matlab/ent_samp.m.
//
// Copyright (c) 2021-2026 Quantitative Analysis Research Core,
// Center for Human Movement Variability, University of Nebraska at Omaha.
// MIT licence. See LICENSE.txt.

#include "quarc/quarc.hpp"

#include <cmath>
#include <limits>
#include <stdexcept>

namespace quarc {

double ent_samp(const double* x, std::size_t n, int dim, double radius) {
    if (dim < 1) throw std::invalid_argument("dim must be at least 1.");
    if (n < static_cast<std::size_t>(dim) + 1)
        throw std::invalid_argument("The time series must have at least m+1 data points.");

    // Both counts run over the same N-dim templates (see ent_samp.m). One pass
    // does both: a pair that matches over dim samples is extended by one more.
    const std::size_t m = static_cast<std::size_t>(dim);
    const std::size_t T = n - m;
    double B = 0, A = 0;
    for (std::size_t i = 0; i + 1 < T; ++i) {
        for (std::size_t j = i + 1; j < T; ++j) {
            std::size_t k = 0;
            while (k < m && std::fabs(x[i + k] - x[j + k]) <= radius) ++k;
            if (k < m) continue;
            B += 1;
            if (std::fabs(x[i + m] - x[j + m]) <= radius) A += 1;
        }
    }
    if (B == 0) return std::numeric_limits<double>::infinity();
    return -std::log(A / B);
}

}  // namespace quarc
