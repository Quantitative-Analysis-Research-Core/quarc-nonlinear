// quarc core: the compiled kernels shared by the MATLAB, R and Python packages.
//
// Rules, so one copy of this code can serve all three:
//   - plain C++17 and the standard library; no MATLAB, R or Python headers
//   - inputs are (const double*, size_t) spans; outputs are returned values
//   - no global state, no printing; errors are std::invalid_argument, which
//     each wrapper turns into its language's own error
//   - the MATLAB function of the same name is the specification; every
//     kernel must reproduce tests/fixtures/matlab_reference.json
//
// Copyright (c) 2021-2026 Quantitative Analysis Research Core,
// Center for Human Movement Variability, University of Nebraska at Omaha.
// MIT licence. See LICENSE.txt.

#pragma once
#include <cstddef>
#include <cstdint>
#include <vector>

namespace quarc {

// Sample entropy (Richman & Moorman 2000), as matlab/ent_samp.m. `radius` is
// absolute here; the wrappers apply the 'prop' scaling by std(x) themselves.
// Returns +Inf when no template pair of length dim matches.
double ent_samp(const double* x, std::size_t n, int dim, double radius);

// Diagonal and vertical line lengths of a recurrence plot, as the loops in
// matlab/line_hist.m. `a` is an n-by-n column-major matrix, nonzero where a
// point recurs. Lines are listed in MATLAB's order: diagonals from offset
// -(n-1) to n-1, then columns 1..n, each from top to bottom.
//   auto-recurrence (cross = false): diagonals with |offset| <= theiler are
//     skipped, and, when theiler > 0, the band |row - col| <= theiler is
//     blanked in each column; no diagonal line at all gives diag = {0}.
//   cross-recurrence (cross = true): theiler is ignored; every diagonal and
//     every column is used, and empty results stay empty.
struct LineHist {
    std::vector<double> diag;
    std::vector<double> vert;
};
LineHist line_hist(const std::uint8_t* a, std::size_t n, std::size_t theiler, bool cross);

}  // namespace quarc
