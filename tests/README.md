# QUARC test suite

Headless. No GUI, no toolboxes beyond base MATLAB, no pytest.

```bash
matlab -batch "addpath('tests/matlab'); run_tests"
matlab -batch "addpath('tests/matlab'); run_tests('Surr')"   # name filter
```

JUnit XML lands in `tests/artifacts/results.xml`.

## The four kinds of test here

**Contract.** Write down what a function *claims* to preserve or return, then
assert it. `quarctest.surrogateContract` measures whether the spectrum, the
distribution, and the variance survived a surrogate generator; the caller
decides which of those the algorithm actually promised. Getting that
distinction right matters more than the measurement: Algorithm 1 owes you an
exact spectrum, Algorithm 2 owes you an exact distribution and only an
*approximate* spectrum. Holding AAFT to exactness would be filing its design as
a bug.

**Known-answer recovery.** Feed a signal whose answer is known analytically and
check the estimator returns it — white noise → DFA α = 0.5, Brownian → 1.5,
fGn at H → α = H. The generators in `quarctest.signals` are written from
scratch and are deliberately independent of the library: if `fgn_sim` were used
to test `dfa`, a matched pair of errors would cancel and the test would pass.

**Structural.** `testHeadless` scans the shipped source for things that make
the library unusable in batch — `dbstop`, `waitbar`, file/function name
mismatches, CR-only line endings. These are invisible to any test that calls
the function inside `try/catch`, which is why they need their own pass.

**Cross-language.** Not active here: the Python port was not moved into this
repository. `tests/fixtures/*.csv` and `matlab_reference.json` are retained
because they are the shared reference data, so the equivalence tests can be
restored when the port is reworked. Regenerate the reference with
`matlab -batch "addpath('tests/matlab'); make_reference"`.

## The characterization suite, which is not a test

`quarctest.characterize` is a separate instrument with a different job. For each
system in the Sprott catalogue it draws R initial conditions i.i.d. around the
published one — Sprott's own `x0` is always realization 1 — runs the full metric
battery on a series from each, and reports the distribution of every metric
across the ensemble.

```
matlab -batch "addpath('tests/matlab'); \
  [S,p] = quarctest.characterize(Fast=true, R=10); \
  quarctest.write_characterization(S, p, Tag='fast')"
```

**It asserts nothing and cannot fail**, which is why it is not named `test*.m`
and is not collected by `run_tests`. Sprott publishes exactly two of these
quantities — the largest Lyapunov exponent and the correlation dimension, the
latter with a stated uncertainty — and those rows are marked `validated` and
carry the reference. Every other row is marked `characterized`: a measurement of
an estimator on a system, with no reference value in existence to judge it
against. The report is the deliverable, and it is committed to `tests/reports/`
so that two runs can be diffed and read.

Three design points worth knowing before extending it:

- **The spread comes from initial conditions drawn on the attractor.** A
  deterministic system from a fixed `x0` has zero variance — the existing
  benchmark reproduces to 5e-15. The ensemble runs the transient from Sprott's
  `x0`, perturbs the resulting on-attractor point by a fraction of the
  attractor's own extent, and relaxes it back. Perturbing `x0` directly was the
  first design and it was wrong: a published `x0` is chosen to be convenient,
  not central, so a symmetric perturbation threw realizations out of the basin
  entirely — 47% of them on the delayed logistic map, 42% on the simplest
  quadratic flow, 30% on the simplest cubic. The survivors would have been a
  sample conditioned on not escaping. For a conservative system there is no
  attractor to relax onto and the realizations are neighbouring orbits, so it is
  a different quantity; `info.hasAttractor` and the category column carry the
  distinction.
- **Sampling is fixed per system.** `sprott_series` normally derives a flow's
  decimation from a pilot run, so a perturbed `x0` could shift the period
  estimate and change `fs` between realizations. `characterize` computes decim
  once from the published `x0` and passes it via `Decim=`, otherwise the scatter
  in every rate-dependent metric would be an artefact of the protocol.
- **Start-up transients are discarded twice, and the length was audited.**
  `sprott_series` drops 1000 iterations for a map and 20000 RK4 steps for a flow
  before collecting anything. `ensemble_ics` adds 5000 steps to reach the
  attractor plus 200 to relax the perturbed point, so a perturbed realization
  has ~25200 steps behind it and realization 1, which uses Sprott's `x0`
  verbatim, has 20000. Measured against each system's own timescale that is
  uneven: heavily decimated flows (`chua`, `double_scroll`, `labyrinth`,
  `rabinovich_fabrikant` and ten others) get only about 5 dominant periods,
  because decimation buys resolution of the period at the cost of transient
  length in integrator steps. Tested rather than assumed — raising the transient
  tenfold to 200000 steps moved `corr_dim` by at most 1.5% and `lyap_ros` by at
  most 3.5% across the six shortest cases, all inside the ensemble spread, so
  there is no evidence of transient contamination at this resolution. The test
  used 24 realizations per system, so it would not resolve an effect below
  roughly 1%. Note also that a conservative system has no attractor to settle
  onto, so for `nose_hoover`, `labyrinth` and `driven_pendulum` the question
  does not arise in the same form.
- **The embedding is a protocol, not a measurement** — see
  `quarctest.embed_policy`. Maps use delay 1 and dimension = state dimension + 1,
  because a map's iterates are its natural coordinates and AMI's first minimum on
  a map is not an embedding delay. Measured: on the logistic map, AMI's answer of
  6 with FNN's dimension gives a correlation dimension of 3.09 against a
  published D2 of exactly 1.0; at delay 1 the same estimator returns 0.94. AMI
  and FNN are still reported as metrics, next to the `embed_delay` and
  `embed_dim` actually used.

## Parameter sweeps

`quarctest.evolve_sweep` measures Wolf's exponent against its renormalisation
interval across every usable map plus eight flows, reusing one series per
realization for all interval values so the comparison is within-realization.

It answers a question the characterization report raised but could not settle:
Wolf's estimator is erratic on maps (0.33x on Arnold's cat, 0.63x on Ikeda at
the default) while behaving on flows. The cause is that `evolve` is fixed at 10
samples for every system. On a flow decimated to ~40 samples per orbit that is
a quarter turn; on a map it is 10 iterations, and a logistic-map pair separates
by 2^10 before the neighbour is replaced -- far outside the linear regime Wolf's
method assumes.

Median absolute log10 error against the published exponent, 200 realizations
per system:

| category | evolve = 1 | evolve = 10 |
|---|---|---|
| noninvertible maps | 0.007 | 0.093 |
| dissipative maps | 0.012 | 0.096 |
| conservative maps | 0.244 | 0.357 |
| flows | 0.033 | 0.036 |

The flows are the control, and they matter: without them a maps-only result
could not separate "maps need a different interval" from "the default is wrong
everywhere". Flows are flat across the whole range and mildly prefer larger
intervals; no flow's best is at 1.

Not universal, and the exceptions are recorded rather than smoothed over.
Conservative maps mostly do not follow the pattern -- Chirikov is best at 30,
the chaotic web returns 0.001 at every interval and is simply not estimable
this way, and Henon's area-preserving map peaks at 0.36. Their exponents are
small and the noise floor dominates, which is the same weakness the Lyapunov
benchmark already reports for conservative systems.

**The library default is unchanged at 10.** Changing it would silently alter
every exponent the library has ever produced, which is a decision for a release
rather than a test suite. The evidence is documented in `lyapunov`'s help so a
caller analysing a map can choose deliberately.

## Rules the harness follows

- **Base MATLAB only.** `corr` is Statistics Toolbox, so the suite uses
  `quarctest.pearson`. A test needing a toolbox must skip, not error.
- **No figures.** `quarctest.sideEffects` counts figures opened and fails the
  test if any survive.
- **`dbclear all` before and after every test.** Ten shipped functions execute
  `dbstop if error`, which is global session state. Under `matlab -batch` an
  uncaught error then hangs the process forever instead of failing it — this
  was measured, not assumed. The runner cannot let that leak between tests.
- **Fixtures and references are inputs, never outputs.** `make_fixtures.py` and
  `make_reference.m` are run by hand. A suite that can rewrite its own expected
  values will happily ratify a regression.

## Two traps that bit this harness during construction

Recorded because both produced *passing* tests over real defects:

- **MATLAB `regexp` has no `\b`.** It means backspace. `'^\s*dbstop\b'` matches
  nothing, so the `dbstop` scan passed while ten violations sat in the tree.
  MATLAB spells word boundaries `\<` and `\>`. `testGrepItselfWorks` now
  asserts the scanner can find something known-present, so a broken pattern
  shows up as a failure rather than a clean bill of health.
- **`strsplit` collapses consecutive delimiters by default**, silently dropping
  blank lines and shifting every reported line number. Pass
  `'CollapseDelimiters', false`.

Also: `functiontests(localfunctions)` collects every local function whose name
*ends* with "test", so a helper called `localSurrogateTest` is loaded as a test
case, fails to accept one argument, and takes the **whole file** out of the
suite with no warning. Helpers here are named `local*` and never `*Test`.

## Adding a function

1. Add a known-answer signal to `quarctest.signals` if the estimator has an
   analytic answer.
2. Write down the contract in the test file's header comment before writing any
   assertion. If you cannot state what the function promises, that is the
   finding.
3. Add a `quarctest.sideEffects` check — errors, figures, `dbstop`, runtime.
4. If the function will also exist in the Python port, add a fixture and a
   line to `make_reference.m`.
