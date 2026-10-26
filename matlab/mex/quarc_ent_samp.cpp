// MEX gateway: v = quarc_ent_samp(x, dim, radius), radius absolute.
// Built by build_mex.m into matlab/private, where ent_samp.m calls it.
//
// Copyright (c) 2021-2026 Quantitative Analysis Research Core,
// Center for Human Movement Variability, University of Nebraska at Omaha.
// MIT licence. See LICENSE.txt.

#include "mex.h"
#include "quarc/quarc.hpp"

#include <stdexcept>

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    if (nrhs != 3) mexErrMsgIdAndTxt("quarc:nargin", "quarc_ent_samp(x, dim, radius)");
    if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) || mxIsSparse(prhs[0]))
        mexErrMsgIdAndTxt("quarc:type", "x must be a real, full double array.");
    try {
        double v = quarc::ent_samp(mxGetDoubles(prhs[0]), mxGetNumberOfElements(prhs[0]),
                                   static_cast<int>(mxGetScalar(prhs[1])),
                                   mxGetScalar(prhs[2]));
        plhs[0] = mxCreateDoubleScalar(v);
    } catch (const std::exception& e) {
        mexErrMsgIdAndTxt("quarc:core", "%s", e.what());
    }
}
