<picture>
  <source media="(prefers-color-scheme: dark)" srcset="assets/quarc-logo-dark.png">
  <source media="(prefers-color-scheme: light)" srcset="assets/quarc-logo-light.png">
  <img alt="QUARC - Quantitative Analysis Research Core" src="assets/quarc-logo-light.png" width="520">
</picture>

# QUARC Nonlinear

Nonlinear time series analysis for MATLAB, from the Quantitative Analysis
Research Core (QUARC) at the Center for Human Movement Variability,
University of Nebraska at Omaha.

This library supersedes the NONAN Library. Function names have changed; the
old names remain available through opt-in shims.

## INSTALLATION

Clone the repository and add `matlab/` to your MATLAB path:

```matlab
addpath(genpath('path/to/quarc-nonlinear/matlab'))
```

To keep the pre-rename function names working, also add the shim folder:

```matlab
addpath('path/to/quarc-nonlinear/deprecated')
```

It sits outside `matlab/` so that `genpath` does not pull it in by accident.
Each shim forwards its arguments unchanged and warns once per session.

## REQUIREMENTS

MATLAB R2019b or later. The `arguments` block and name-value syntax used by
the newer functions need R2019b; `ami`, `lyapunov` and the RQA family use it.

Most functions run on **base MATLAB with no toolboxes**. The exceptions are
noted per function in their help text.

## USAGE

A complete pass over a chaotic series: generate it, choose an embedding,
then quantify it.

```matlab
% Lorenz attractor, x component, sampled at 100 Hz
fs = 100;
[~, y] = chaos_library('Lorenz', 0:1/fs:100, [1 1 1], [10 28 8/3]);
x = y(:,1);

% Embedding parameters: delay from AMI, dimension from false nearest neighbours
delay = ami(x, 100);
dim   = fnn(x, delay, 10, 15, 2, 1);

% Recurrence quantification at 2.5% recurrence
[rp, results] = rqa(x, delay, dim, "rec", 2.5);
results.DET     % determinism
results.LAM     % laminarity

% Largest Lyapunov exponent, nats per unit time
lambda = lyapunov(x, fs, delay=delay, dim=dim);

% Sample entropy, and the multiscale family out to 20 scales
se   = ent_samp(x, dim, 0.2);
rcmse = ent_ms_plus(x, 20, dim, 0.2);
```

Every function documents its arguments and returns in its own help text
(`help rqa`, `doc lyapunov`).

## TESTS

```bash
matlab -batch "addpath('tests/matlab'); run_tests"
```

Headless, base MATLAB only, exits nonzero on failure and writes JUnit XML to
`tests/artifacts/`. Filter by name with `run_tests('Surr')`. See
`tests/README.md` for how the suite is organised.

The 62 systems of Sprott (2003) Appendix A are catalogued with their
published exponents in `tests/fixtures/sprott_appendix_a.md` and serve as the
reference set. One protocol is applied to every system, with nothing tuned
case by case.

`corr_dim` is estimated on every usable system as part of the suite.

The Lyapunov tests run in three layers rather than one sweep: invariances
that need no reference value at all, exact exponents from maps where lambda
is a theorem (skew tent, logistic), and deliberately loose order-of-magnitude
checks against published values for Henon, Rossler and Lorenz. A tight
assertion against a numerical reference for a flow would test the choice of
sampling rate, delay and scaling region as much as the code.

The full head-to-head comparison of the two estimators across the catalogue
is a separate benchmark, `quarctest.lye_benchmark`, run on demand rather than
in the suite; its per-system results are checked in as
`tests/fixtures/lye_benchmark_results.csv`. Over the 55 usable systems the
median ratio to the published exponent is **0.96 for Wolf's method** and
**0.85 for Rosenstein's**. Both are weakest on conservative systems, where
lambda is small and the noise floor dominates. Report which method produced
a published exponent.

## CHARACTERIZATION

Beyond pass/fail, the library is characterized across the whole catalogue:
for each system, an ensemble of initial conditions drawn around the published
one, the full metric battery run on each realization, and the distribution of
every metric reported with its spread.

```matlab
[S, per] = quarctest.characterize(Fast=true, R=10);
quarctest.write_characterization(S, per)
```

The report lands in `tests/reports/`. Two metrics carry published references
and are marked `validated` — the largest Lyapunov exponent, and the
correlation dimension with Sprott's own stated uncertainty. Everything else
is marked `characterized`: measured and reported, with no reference value in
existence to judge it against.

The catalogue is the starting point rather than the scope. The runner takes
any system that can produce a scalar observable, so higher-dimensional,
biological and biomechanical systems extend it without changing the design.

## WHAT IS INCLUDED

**Embedding.** Delay from average mutual information (`ami`, by histogram or
kernel density estimate), dimension from false nearest neighbours (`fnn`), and
phase space reconstruction (`psr`, `embed`).

**Recurrence quantification.** Single series (`rqa`), cross (`crqa`), joint
(`jrqa`), and multidimensional (`mdrqa`), with radius selection for a target
percent recurrence (`set_radius`), line-length histograms (`line_hist`),
recurrence plot entropy (`ent_weighted`), and plotting (`rqa_plot`).

**Entropy.** Sample, approximate, permutation, and symbolic entropy, their
cross-series forms, and the multiscale family — refined composite, composite,
multiscale, multiscale fuzzy, and generalized (`ent_samp`, `ent_ap`,
`ent_permu`, `ent_symbolic`, `ent_xsamp`, `ent_xap`, `ent_ms_plus`), under one
`ent` entry point.

**Divergence and dimension.** Largest Lyapunov exponent by Rosenstein's and
Wolf's methods, in nats per unit time and comparable across the two
(`lyapunov`, `lye_r`, `lye_w`), and correlation dimension (`corr_dim`).

**Scaling.** Detrended fluctuation analysis (`dfa`) and simulation of
fractional Gaussian noise at a chosen Hurst exponent (`fgn_sim`).

**Surrogates.** Theiler shuffle, Fourier, and amplitude-adjusted Fourier
surrogates (`surr_theiler`), and pseudo-periodic surrogates with automatic
noise radius (`surr_pseudo_periodic`, `surr_find_rho`).

**Coordination.** Continuous and discrete relative phase between two series
(`rel_phase_cont`, `rel_phase_disc`).

**Chaotic systems.** Named attractors integrated from their governing
equations, with the time scale, initial conditions and parameters under the
caller's control (`chaos_library`).

## ARGUMENT NAMES

One name per concept across the library. See `NAMING.md`.

| concept | name | was |
|---|---|---|
| time series | `x` (second series `y`) | `data`, `X`, `DATA`, `y`, `z`, `S1`/`S2`, `data1`/`data2` |
| embedding delay | `delay` | `tau` |
| embedding dimension | `dim` | `m`, `M`, `de`, `MaxDim` |
| sampling frequency | `fs` | `Fs`, `samprate` |
| tolerance / radius | `radius` | `r`, `R` |
| maximum lag | `maxlag` | `L` |

These are positional arguments, so **existing calls keep working** — only the
names in the documentation and in `arguments` blocks changed. The exception is
`lyapunov`, where the name-value option `Tau=` became `Delay=`.

`dim` and `radius` correspond to `m` and `r` in Richman & Moorman (2000); each
function header states the mapping so the code can be read alongside the
papers.

## LICENCE

MIT. See `LICENSE.txt`.

Copyright (c) 2021-2026 Quantitative Analysis Research Core, Center for Human
Movement Variability, University of Nebraska at Omaha.

`matlab/embed.m` is Copyright (c) 1994 Kevin Judd and is not covered by the
MIT grant; its original notice is retained in the file. See
`THIRD-PARTY-NOTICES.txt`.

## CONTACT

Please contact quarc@unomaha.edu regarding any questions or troubleshooting.
