// MEX gateway: [diag_hist, vertical_hist] = quarc_line_hist(a, theiler, cross)
// a is a square logical recurrence matrix. Built by build_mex.m into
// matlab/private, where line_hist.m calls it.
//
// Copyright (c) 2021-2026 Quantitative Analysis Research Core,
// Center for Human Movement Variability, University of Nebraska at Omaha.
// MIT licence. See LICENSE.txt.

#include "mex.h"
#include "quarc/quarc.hpp"

#include <algorithm>
#include <cstdint>

namespace {
// A 1-by-k row, or 0-by-0 when empty, as line_hist.m grows its outputs.
mxArray* row(const std::vector<double>& v) {
    if (v.empty()) return mxCreateDoubleMatrix(0, 0, mxREAL);
    mxArray* m = mxCreateDoubleMatrix(1, v.size(), mxREAL);
    std::copy(v.begin(), v.end(), mxGetDoubles(m));
    return m;
}
}  // namespace

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    if (nrhs != 3) mexErrMsgIdAndTxt("quarc:nargin", "quarc_line_hist(a, theiler, cross)");
    if (!mxIsLogical(prhs[0]) || mxGetM(prhs[0]) != mxGetN(prhs[0]))
        mexErrMsgIdAndTxt("quarc:type", "a must be a square logical matrix.");
    try {
        static_assert(sizeof(mxLogical) == 1, "mxLogical is one byte");
        const double th = mxGetScalar(prhs[1]);
        auto h = quarc::line_hist(reinterpret_cast<const std::uint8_t*>(mxGetLogicals(prhs[0])),
                                  mxGetM(prhs[0]), th > 0 ? static_cast<std::size_t>(th) : 0,
                                  mxIsLogicalScalarTrue(prhs[2]));
        plhs[0] = row(h.diag);
        if (nlhs > 1) plhs[1] = row(h.vert);
    } catch (const std::exception& e) {
        mexErrMsgIdAndTxt("quarc:core", "%s", e.what());
    }
}
