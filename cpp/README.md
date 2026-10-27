# C++ core

One implementation of the slow kernels, called from all three packages:

| Package | Wrapper | Build |
|---|---|---|
| MATLAB | `matlab/mex/quarc_*.cpp` (MEX) | `addpath matlab/mex; build_mex` |
| R | `rpkg/src/rcpp_*.cpp` (Rcpp) | `devtools::install("rpkg")` |
| Python | `pypkg/src/bindings.cpp` (nanobind) | `pip install ./pypkg` |

Rules: plain C++17, built with `-ffp-contract=off` (no fused multiply-add, so
`e*e + s` rounds twice, as in MATLAB, on every compiler),, no MATLAB/R/Python headers, `(const double*, size_t)` in,
values out, errors as `std::invalid_argument`. The MATLAB function of the same
name is the specification; each kernel is tested in all three languages
against `tests/fixtures/matlab_reference.json`.

R and Python build from their own directories, so each holds a copy of this
folder. Edit here, then run `tools/sync_core.sh`; both test suites fail if a
copy is stale.

Kernels so far: `ent_samp`, `line_hist`, `corr_dim`, `fnn`, `lye_r` (MATLAB); `ent_samp` (R, Python).

`QUARC_NO_MEX=1` makes MATLAB run its own code; `tests/matlab/testCore.m` compares the two.
