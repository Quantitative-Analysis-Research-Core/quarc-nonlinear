# quarcnonlinear (Python)

Python bindings to the QUARC Nonlinear C++ core (`../cpp`), the same compiled
code the MATLAB (MEX) and R (Rcpp) packages call.

    pip install ./pypkg          # from the repository root
    pytest pypkg/tests

`src/core` is a copy of `../cpp` made by `tools/sync_core.sh`; edit `cpp/`,
never the copy.
