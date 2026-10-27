// MEX gateway for corr_dim.m's pair loop:
//   [dmin, dmax] = quarc_corr_dim(x, dim, delay)
//   counts       = quarc_corr_dim(x, dim, delay, edges)   (1-by-numel(edges))
//
// Copyright (c) 2021-2026 Quantitative Analysis Research Core,
// Center for Human Movement Variability, University of Nebraska at Omaha.
// MIT licence. See LICENSE.txt.

#include "mex.h"
#include "quarc/quarc.hpp"

#include <algorithm>

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    if (nrhs != 3 && nrhs != 4) mexErrMsgIdAndTxt("quarc:nargin", "quarc_corr_dim(x, dim, delay[, edges])");
    if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) || mxIsSparse(prhs[0]))
        mexErrMsgIdAndTxt("quarc:type", "x must be a real, full double array.");
    try {
        const double* x = mxGetDoubles(prhs[0]);
        const std::size_t n = mxGetNumberOfElements(prhs[0]);
        const int dim = static_cast<int>(mxGetScalar(prhs[1]));
        const std::size_t delay = static_cast<std::size_t>(mxGetScalar(prhs[2]));
        if (nrhs == 3) {
            double lo, hi;
            quarc::corr_dim_range(x, n, dim, delay, lo, hi);
            plhs[0] = mxCreateDoubleScalar(lo);
            if (nlhs > 1) plhs[1] = mxCreateDoubleScalar(hi);
        } else {
            const std::size_t nb = mxGetNumberOfElements(prhs[3]);
            auto c = quarc::corr_dim_counts(x, n, dim, delay, mxGetDoubles(prhs[3]), nb);
            plhs[0] = mxCreateDoubleMatrix(1, nb, mxREAL);
            std::copy(c.begin(), c.end(), mxGetDoubles(plhs[0]));
        }
    } catch (const std::exception& e) {
        mexErrMsgIdAndTxt("quarc:core", "%s", e.what());
    }
}
