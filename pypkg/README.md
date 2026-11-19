# quarcnonlinear (Python)

The whole QUARC Nonlinear library in Python: the same 39 functions as the R
package (`../R`), with MATLAB's names, arguments and defaults. Every function
with deterministic output is checked against MATLAB's answers on the shared
fixtures (`../tests/fixtures`, relative tolerance 1e-10); functions that
draw random numbers (`fgn_sim`, the surrogates) take an `rng` argument (a
seed or `numpy.random.Generator`) and are tested by their properties.

The heavy kernels (`ent_samp`, `fnn`, `corr_dim`, `lye_r`, `line_hist`) call
the shared C++ core (`../cpp`), the same compiled code the MATLAB (MEX) and
R (Rcpp) packages use. Requires numpy and scipy.

Conventions: several MATLAB outputs come back as a tuple in MATLAB's order; a
MATLAB struct or table comes back as a dict. Indices in results are 0-based;
delays, dimensions and counts keep their MATLAB values. Plotting is not
ported.

    pip install ./pypkg          # from the repository root
    pytest pypkg/tests

`src/core` is a copy of `../cpp` made by `tools/sync_core.sh`; edit `cpp/`,
never the copy.
