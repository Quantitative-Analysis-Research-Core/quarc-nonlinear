// Line-length histograms of a recurrence plot. See quarc.hpp and matlab/line_hist.m.
//
// Copyright (c) 2021-2026 Quantitative Analysis Research Core,
// Center for Human Movement Variability, University of Nebraska at Omaha.
// MIT licence. See LICENSE.txt.

#include "quarc/quarc.hpp"

namespace quarc {

namespace {

// Append the lengths of the runs of nonzero values in a(i0 + t*step),
// t = 0..len-1, to out, in order.
void runs(const std::uint8_t* a, std::size_t i0, std::ptrdiff_t step, std::size_t len,
          std::vector<double>& out) {
    double run = 0;
    for (std::size_t t = 0; t < len; ++t) {
        if (a[i0 + t * step]) {
            run += 1;
        } else if (run > 0) {
            out.push_back(run);
            run = 0;
        }
    }
    if (run > 0) out.push_back(run);
}

}  // namespace

LineHist line_hist(const std::uint8_t* a, std::size_t n, std::size_t theiler, bool cross) {
    LineHist h;
    if (n == 0) {
        if (!cross) h.diag.push_back(0);
        return h;
    }
    const std::size_t w = cross ? 0 : theiler;
    const std::ptrdiff_t down = static_cast<std::ptrdiff_t>(n) + 1;  // (i,j) -> (i+1,j+1)

    // Offsets -(n-1)..-1: below the main diagonal, starting at row k, column 0.
    for (std::size_t k = n - 1; k >= 1; --k) {
        if (cross || k > w) runs(a, k, down, n - k, h.diag);
    }
    // Offsets 0..n-1: on and above it, starting at row 0, column k.
    for (std::size_t k = 0; k < n; ++k) {
        if (cross || k > w) runs(a, k * n, down, n - k, h.diag);
    }
    if (!cross && h.diag.empty()) h.diag.push_back(0);

    // Columns, with the Theiler band blanked (auto-recurrence only).
    for (std::size_t j = 0; j < n; ++j) {
        const std::uint8_t* col = a + j * n;
        if (cross || theiler == 0) {
            runs(col, 0, 1, n, h.vert);
            continue;
        }
        const std::size_t lo = j > theiler ? j - theiler : 0;
        const std::size_t hi = j + theiler < n ? j + theiler : n - 1;  // inclusive
        runs(col, 0, 1, lo, h.vert);
        if (hi + 1 < n) runs(col, hi + 1, 1, n - hi - 1, h.vert);
    }
    return h;
}

}  // namespace quarc
