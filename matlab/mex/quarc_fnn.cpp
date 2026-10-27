// MEX gateway for one dimension of fnn.m:
//   frac = quarc_fnn(x, delay, j, n, Rtol, Atol, RA)
//
// Copyright (c) 2021-2026 Quantitative Analysis Research Core,
// Center for Human Movement Variability, University of Nebraska at Omaha.
// MIT licence. See LICENSE.txt.

#include "mex.h"
#include "quarc/quarc.hpp"

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    if (nrhs != 7) mexErrMsgIdAndTxt("quarc:nargin", "quarc_fnn(x, delay, j, n, Rtol, Atol, RA)");
    if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) || mxIsSparse(prhs[0]))
        mexErrMsgIdAndTxt("quarc:type", "x must be a real, full double array.");
    try {
        auto sz = [&](int k) { return static_cast<std::size_t>(mxGetScalar(prhs[k])); };
        double f = quarc::fnn_false_fraction(mxGetDoubles(prhs[0]), mxGetNumberOfElements(prhs[0]),
                                             sz(1), sz(2), sz(3), mxGetScalar(prhs[4]),
                                             mxGetScalar(prhs[5]), mxGetScalar(prhs[6]));
        plhs[0] = mxCreateDoubleScalar(f);
    } catch (const std::exception& e) {
        mexErrMsgIdAndTxt("quarc:core", "%s", e.what());
    }
}
