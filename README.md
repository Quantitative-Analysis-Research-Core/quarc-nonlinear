<picture>
  <source media="(prefers-color-scheme: dark)" srcset="assets/quarc-logo-dark.png">
  <source media="(prefers-color-scheme: light)" srcset="assets/quarc-logo-light.png">
  <img alt="QUARC - Quantitative Analysis Research Core" src="assets/quarc-logo-light.png" width="520">
</picture>

# QUARC Nonlinear

Nonlinear time series analysis for MATLAB, from the Quantitative Analysis
Research Core (QUARC) at the Center for Human Movement Variability,
University of Nebraska at Omaha.

This supersedes the NONAN Library. QUARC is the Core's current name, and this
repository is where ongoing development happens, including changes that break
compatibility with the NONAN releases.

## INSTALLATION

Clone the repository and add `matlab/` to your MATLAB path:

```matlab
addpath(genpath('path/to/quarc-nonlinear/matlab'))
```

To keep pre-rename function names working, also add the shim folder:

```matlab
addpath('path/to/quarc-nonlinear/matlab/deprecated')
```

Each shim forwards its arguments unchanged and warns once per session.

## REQUIREMENTS

MATLAB R2019b or later. The `arguments` block and name-value syntax used by
the newer functions need R2019b; `ami`, `lyapunov` and the RQA family use it.

Most functions run on **base MATLAB with no toolboxes**. The exceptions are
noted per function in their help text. Functions rewritten during the audit
(`ami`, `ami_histogram`, `ami_kde`) had their Statistics Toolbox dependencies
removed and now run on base MATLAB.

## TESTS

```bash
matlab -batch "addpath('tests/matlab'); run_tests"
```

Headless, base MATLAB only, exits nonzero on failure and writes JUnit XML to
`tests/artifacts/`. Filter by name with `run_tests('Surr')`. See
`tests/README.md` for how the suite is organised.

The suite includes a benchmark of both Lyapunov estimators against the 62
systems of Sprott (2003) Appendix A; see `tests/fixtures/`.

## CHANGES FROM NONAN

Function and file names are now `lower_snake_case` with no date suffixes.
Old names remain available through `matlab/deprecated/`.

Corrected during the audit that preceded this repository:

- `surr_theiler` algorithm 1 now preserves the power spectrum exactly.
  Spectral error fell from ~0.59 to 3e-16 and the standard deviation ratio
  from 0.71 to 1.0000.
- `lye_r` is now scale invariant. A hard-coded exclusion marker of 1e5 meant
  the exponent collapsed to 8% of its correct value once distances exceeded
  that, silently.
- `lye_r` memory is now O(N) rather than O(N^2): 287 MB to 0.10 MB at
  N = 6000.
- `surr_find_rho` always returns a value. It previously failed to assign its
  output on 18-30% of calls depending on the series.
- `dbstop if error` removed from ten functions. It is global session state
  and made batch runs hang rather than fail.
- `waitbar` removed. It required a display and contributed nothing to the
  result.
- `ami` replaces `AMI_Stergiou` and `AMI_Thomas` with one entry point and an
  `Algorithm` argument. The histogram estimator's inverted `Bins` guard, bin
  off-by-one, and non-strict minimum test are fixed; the kernel estimator is
  7-10x faster and numerically identical.
- `lyapunov` provides the same wrapper pattern for the Lyapunov estimators
  and accepts a pre-built phase space.

The Python port has not been moved here. It has known defects and needs its
own rework; it remains in the NONAN repository for now.

### FILES

This is a list of the included functions and the full name of the methods.

All function and file names are lower_snake_case. The previous names still
work through shims in `matlab/deprecated/`; add that folder to your path if
you need them, and they will warn once per session.

| function | description |
|---|---|
| `ami` | Average mutual information versus lag, for choosing an embedding delay. Wrapper over `ami_histogram` and `ami_kde`. |
| `ami_histogram` | AMI by equal-width joint histogram. |
| `ami_kde` | AMI by Gaussian kernel density estimate. |
| `chaos_library` | Systems of differential equations that produce chaotic attractors. |
| `corr_dim` | Correlation dimension. |
| `crqa` | Cross recurrence quantification analysis. |
| `dfa` | Detrended fluctuation analysis. |
| `embed` | Delay embedding of a time series. |
| `ent_ap` | Approximate entropy. |
| `ent_ms_plus` | Refined composite multiscale, composite multiscale, multiscale, multiscale fuzzy, and generalized multiscale entropy. |
| `ent_permu` | Permutation entropy, log base 2. |
| `ent_samp` | Sample entropy. |
| `ent_symbolic` | Symbolic entropy. |
| `ent_weighted` | Weighted entropy of a recurrence plot. |
| `ent_xap` | Cross approximate entropy between two series. |
| `ent_xsamp` | Cross sample entropy between two series. |
| `fgn_sim` | Simulate fractional Gaussian noise at a specified Hurst exponent. |
| `fnn` | Embedding dimension by false nearest neighbours. |
| `jrqa` | Joint recurrence quantification analysis. |
| `line_hist` | Diagonal and vertical line histograms of a recurrence plot. |
| `lye_r` | Largest Lyapunov exponent, Rosenstein's method. Returns the divergence curve. |
| `lye_w` | Largest Lyapunov exponent, Wolf's method. Returns bits per unit time. |
| `mdrqa` | Multidimensional recurrence quantification analysis. |
| `psr` | Phase space reconstruction. |
| `rel_phase_cont` | Continuous relative phase between two cyclic series. |
| `rel_phase_disc` | Discrete relative phase between two series. |
| `rqa` | Recurrence quantification analysis. |
| `rqa_legacy` | Previous combined RQA/cRQA/jRQA/mdRQA entry point, kept for reproducibility. |
| `rqa_plot` | Plot a recurrence plot and its statistics. |
| `set_radius` | Find the radius giving a target percent recurrence. |
| `surr_find_rho` | Optimal noise radius for a pseudo-periodic surrogate. |
| `surr_pseudo_periodic` | Pseudo-periodic surrogate, using the radius from `surr_find_rho`. |
| `surr_theiler` | Theiler surrogates: shuffle, Fourier transform, and amplitude-adjusted Fourier transform. |

### TESTS

```
matlab -batch "addpath('tests/matlab'); run_tests"
python3 tests/python/run_tests.py
```

Headless, base MATLAB only, exits nonzero on failure. See `tests/README.md`.

## LICENCE

MIT. See `LICENSE.txt`.

Copyright (c) 2021-2026 Quantitative Analysis Research Core, Center for Human
Movement Variability, University of Nebraska at Omaha.

`matlab/embed.m` is Copyright (c) 1994 Kevin Judd and is not covered by the
MIT grant; its original notice is retained in the file. See
`THIRD-PARTY-NOTICES.txt`.

## CONTACT

Please contact bmchnonan@unomaha.edu regarding any questions or
troubleshooting.
