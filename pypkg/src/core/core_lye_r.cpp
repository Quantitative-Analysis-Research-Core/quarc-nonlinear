// Rosenstein divergence curve. See quarc.hpp and matlab/lye_r.m.
//
// The distances from point i to every point are computed column by column
// into one buffer (a loop the compiler vectorises), in lye_r.m's order, and
// compared squared. The nearest neighbour is then the first j whose distance
// equals the smallest, as min() picks it: equal square roots of different
// squares are possible, so ties are settled on sqrt, not on the squares.
//
// Copyright (c) 2021-2026 Quantitative Analysis Research Core,
// Center for Human Movement Variability, University of Nebraska at Omaha.
// MIT licence. See LICENSE.txt.

#include "quarc/quarc.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace quarc {

LyeR lye_r_divergence(const double* Y, std::size_t M, std::size_t D, std::size_t theiler) {
    LyeR out;
    out.nn.assign(M, 0.0);
    out.avg.assign(M, 0.0);
    if (M == 0) return out;
    const double inf = std::numeric_limits<double>::infinity();

    // Nearest neighbour outside the Theiler window. min() takes the first of
    // equal minima, skips NaN, and returns the first index when all are Inf
    // (or all NaN).
    std::vector<std::size_t> nn(M, 0);
    std::vector<double> s(M);
    for (std::size_t i = 0; i < M; ++i) {
        std::fill(s.begin(), s.end(), 0.0);
        double* __restrict sp = s.data();
        for (std::size_t c = 0; c < D; ++c) {
            const double yi = Y[i + c * M];
            const double* __restrict col = Y + c * M;
            for (std::size_t j = 0; j < M; ++j) {
                const double e = yi - col[j];
                sp[j] += e * e;
            }
        }
        const std::size_t lo = i > theiler ? i - theiler : 0;
        const std::size_t hi = std::min(i + theiler, M - 1);
        for (std::size_t j = lo; j <= hi; ++j) sp[j] = inf;

        double smin = inf;
        for (std::size_t j = 0; j < M; ++j) smin = sp[j] < smin ? sp[j] : smin;  // NaN skipped
        std::size_t arg = 0;
        if (smin < inf) {
            const double dmin = std::sqrt(smin);
            for (std::size_t j = 0; j < M; ++j) {
                if (sp[j] <= smin * (1 + 4 * std::numeric_limits<double>::epsilon()) &&
                    std::sqrt(sp[j]) == dmin) {
                    arg = j;
                    break;
                }
            }
        }
        nn[i] = arg;
        out.nn[i] = static_cast<double>(arg);
    }

    // Mean log divergence of each pair (i, nn[i]) followed forward together.
    const std::size_t T = M - 1;
    std::vector<double> sum(T, 0.0), cnt(T, 0.0), d(T);
    for (std::size_t i = 0; i + 1 < M; ++i) {
        const std::size_t steps = std::min(M - 1 - nn[i], M - 1 - i);  // EndITL
        double* __restrict dp = d.data();
        std::fill(dp, dp + steps, 0.0);
        for (std::size_t c = 0; c < D; ++c) {
            const double* __restrict a = Y + c * M + i + 1;
            const double* __restrict b = Y + c * M + nn[i] + 1;
            for (std::size_t t = 0; t < steps; ++t) {
                const double e = a[t] - b[t];
                dp[t] += e * e;
            }
        }
        for (std::size_t t = 0; t < steps; ++t) {
            if (dp[t] > 0) {
                sum[t] += std::log(std::sqrt(dp[t]));
                cnt[t] += 1;
            }
        }
    }
    for (std::size_t t = 0; t < T; ++t) {
        if (cnt[t] > 0) {
            out.avg[t] = sum[t] / cnt[t];
            out.count = t + 1;
        }
    }
    return out;
}

}  // namespace quarc
