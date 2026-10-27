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

// Correlation dimension, as matlab/corr_dim.m: the pair loop over embedded
// vectors i < j with j - i > dim*delay. corr_dim_range gives the smallest and
// largest pair distance; corr_dim_counts the histc counts of all pair
// distances against `edges` (nb ascending edges; the last bin counts
// distances equal to the last edge). Embedding as matlab/embed.m.
void corr_dim_range(const double* x, std::size_t len, int dim, std::size_t delay,
                    double& dmin, double& dmax);
std::vector<double> corr_dim_counts(const double* x, std::size_t len, int dim,
                                    std::size_t delay, const double* edges, std::size_t nb);

// False nearest neighbours at one dimension, as one pass of the j loop in
// matlab/fnn.m, including its k-d tree (512-point leaves) step for step so
// that ties resolve as in MATLAB. Returns sum(L)/n, the fraction of the n
// points that are false neighbours when dimension j+1 is added. ra is std(x).
double fnn_false_fraction(const double* x, std::size_t len, std::size_t delay, std::size_t j,
                          std::size_t n, double rtol, double atol, double ra);

// Rosenstein's divergence curve, as the two loops in matlab/lye_r.m. Y is the
// M-by-D column-major state matrix. nn[i] is the (0-based) nearest neighbour of
// point i outside the Theiler window; avg[i] the mean log divergence after i+1
// steps (0 where no pair contributes); count is the length of MATLAB's
// AveLnDiv, i.e. one past the last i with a contributing pair (0 if none).
struct LyeR {
    std::vector<double> nn;
    std::vector<double> avg;
    std::size_t count = 0;
};
LyeR lye_r_divergence(const double* Y, std::size_t M, std::size_t D, std::size_t theiler);

}  // namespace quarc
