// Correlation-sum pair loop. See quarc.hpp and matlab/corr_dim.m.
//
// For speed the distances to vector i are computed for all later j at once
// (a loop the compiler vectorises) and compared as squared distances. Both
// steps are exact: squares are summed in corr_dim.m's order, sqrt is
// monotonic, and each histogram edge e is replaced by the smallest s with
// sqrt(s) >= e, so a squared distance lands in the bin its distance would.
//
// Copyright (c) 2021-2026 Quantitative Analysis Research Core,
// Center for Human Movement Variability, University of Nebraska at Omaha.
// MIT licence. See LICENSE.txt.

#include "quarc/quarc.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace quarc {

namespace {

struct Pairs {
    const double* x;
    std::size_t dim, delay, n, k;
};

Pairs pairs(const double* x, std::size_t len, int dim, std::size_t delay) {
    if (dim < 1) throw std::invalid_argument("dim must be at least 1.");
    if (delay < 1) throw std::invalid_argument("delay must be at least 1.");
    const std::size_t span = static_cast<std::size_t>(dim - 1) * delay;
    if (len <= span) throw std::invalid_argument("series too short for this embedding.");
    const std::size_t d = static_cast<std::size_t>(dim);
    return {x, d, delay, len - span, d * delay};
}

// For i = 1..n-k-1 (1-based), squared distances to j = i+k+1..n into s, then f(s, m).
// embed.m puts the newest sample first, so row r of vector i is x[i + (dim-1-r)*delay];
// rows are summed in that order, as sum(...,1) in corr_dim.m.
template <class F>
void for_rows(const Pairs& p, std::vector<double>& s, F&& f) {
    for (std::size_t i = 0; i + p.k + 1 < p.n; ++i) {
        const std::size_t j0 = i + p.k + 1, m = p.n - j0;
        s.assign(m, 0.0);
        double* __restrict sp = s.data();
        for (std::size_t r = 0; r < p.dim; ++r) {
            const std::size_t off = (p.dim - 1 - r) * p.delay;
            const double xi = p.x[i + off];
            const double* __restrict xj = p.x + j0 + off;
            for (std::size_t t = 0; t < m; ++t) {
                const double e = xj[t] - xi;
                sp[t] += e * e;
            }
        }
        f(sp, m);
    }
}

// Smallest s >= 0 with sqrt(s) >= e (sqrt is correctly rounded and monotonic).
double sq_threshold(double e) {
    if (!(e > 0)) return e <= 0 ? -std::numeric_limits<double>::infinity() : e;
    double s = e * e;
    while (s > 0 && std::sqrt(s) >= e) s = std::nextafter(s, 0.0);
    while (std::sqrt(s) < e) s = std::nextafter(s, std::numeric_limits<double>::infinity());
    return s;
}

// Largest s with sqrt(s) <= e.
double sq_ceiling(double e) {
    double s = e * e;
    while (std::sqrt(s) > e) s = std::nextafter(s, 0.0);
    while (std::sqrt(std::nextafter(s, std::numeric_limits<double>::infinity())) <= e)
        s = std::nextafter(s, std::numeric_limits<double>::infinity());
    return s;
}

}  // namespace

void corr_dim_range(const double* x, std::size_t len, int dim, std::size_t delay,
                    double& dmin, double& dmax) {
    const Pairs p = pairs(x, len, dim, delay);
    double smax = 0, smin = std::numeric_limits<double>::infinity();
    std::vector<double> s;
    // NaN compares false, so it is skipped, as max() and min() skip it.
    for_rows(p, s, [&](const double* sp, std::size_t m) {
        for (std::size_t t = 0; t < m; ++t) {
            if (sp[t] > smax) smax = sp[t];
            if (sp[t] < smin) smin = sp[t];
        }
    });
    dmax = std::sqrt(smax);
    dmin = std::sqrt(smin);
}

std::vector<double> corr_dim_counts(const double* x, std::size_t len, int dim,
                                    std::size_t delay, const double* edges, std::size_t nb) {
    const Pairs p = pairs(x, len, dim, delay);
    std::vector<double> c(nb, 0.0);
    if (nb == 0) return c;
    // histc: bin b holds edges[b] <= d < edges[b+1]; the last holds d == edges[nb-1].
    std::vector<double> th(nb);
    for (std::size_t b = 0; b < nb; ++b) th[b] = sq_threshold(edges[b]);
    const double lo = th[0], last = th[nb - 1], hi = sq_ceiling(edges[nb - 1]);
    const double* tb = th.data();
    std::vector<double> s;
    for_rows(p, s, [&](const double* sp, std::size_t m) {
        for (std::size_t t = 0; t < m; ++t) {
            const double v = sp[t];
            if (!(v >= lo) || v > hi) continue;  // below range, above it, or NaN
            if (v >= last) { c[nb - 1] += 1; continue; }
            // Last threshold <= v, by a branch-free binary search.
            std::size_t base = 0, n = nb - 1;
            while (n > 1) {
                const std::size_t half = n / 2;
                base = (tb[base + half] <= v) ? base + half : base;
                n -= half;
            }
            c[base] += 1;
        }
    });
    return c;
}

}  // namespace quarc
