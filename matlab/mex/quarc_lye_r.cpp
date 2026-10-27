// MEX gateway for lye_r.m's two loops:
//   [IND2, avg, count] = quarc_lye_r(Y, theiler)
// IND2 is 1-based (M-by-1); avg is out(:,3) (M-by-1); AveLnDiv = avg(1:count)'.
//
// Copyright (c) 2021-2026 Quantitative Analysis Research Core,
// Center for Human Movement Variability, University of Nebraska at Omaha.
// MIT licence. See LICENSE.txt.

#include "mex.h"
#include "quarc/quarc.hpp"

#include <algorithm>

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    if (nrhs != 2) mexErrMsgIdAndTxt("quarc:nargin", "quarc_lye_r(Y, theiler)");
    if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) || mxIsSparse(prhs[0]))
        mexErrMsgIdAndTxt("quarc:type", "Y must be a real, full double matrix.");
    try {
        const std::size_t M = mxGetM(prhs[0]), D = mxGetN(prhs[0]);
        const double th = mxGetScalar(prhs[1]);
        auto r = quarc::lye_r_divergence(mxGetDoubles(prhs[0]), M, D,
                                         th > 0 ? static_cast<std::size_t>(th) : 0);
        plhs[0] = mxCreateDoubleMatrix(M, 1, mxREAL);
        double* p = mxGetDoubles(plhs[0]);
        for (std::size_t i = 0; i < M; ++i) p[i] = r.nn[i] + 1;
        if (nlhs > 1) {
            plhs[1] = mxCreateDoubleMatrix(M, 1, mxREAL);
            std::copy(r.avg.begin(), r.avg.end(), mxGetDoubles(plhs[1]));
        }
        if (nlhs > 2) plhs[2] = mxCreateDoubleScalar(static_cast<double>(r.count));
    } catch (const std::exception& e) {
        mexErrMsgIdAndTxt("quarc:core", "%s", e.what());
    }
}
