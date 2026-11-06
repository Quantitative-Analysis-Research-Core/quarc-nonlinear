# quarcnonlinear (Python)

Python bindings to the QUARC Nonlinear C++ core (`../cpp`), the same compiled
code the MATLAB (MEX) and R (Rcpp) packages call. So far: `ent_samp`, `fnn`,
`corr_dim`, `lye_r`, `line_hist` and `psr`, each checked against MATLAB's
answers on the shared fixtures. Indices in results are 0-based.

    pip install ./pypkg          # from the repository root
    pytest pypkg/tests

`src/core` is a copy of `../cpp` made by `tools/sync_core.sh`; edit `cpp/`,
never the copy.
