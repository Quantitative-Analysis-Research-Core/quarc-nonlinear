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

- **The spread comes from initial conditions.** A deterministic system from a
  fixed `x0` has zero variance — the existing benchmark reproduces to 5e-15. For
  a dissipative system the realizations are finite windows of one attractor
  started at different phases, so the spread is the estimator's sampling
  variability. For a conservative system there is no attractor and the
  realizations are distinct orbits, so it is a different quantity; the category
  column carries the distinction.
- **Sampling is fixed per system.** `sprott_series` normally derives a flow's
  decimation from a pilot run, so a perturbed `x0` could shift the period
  estimate and change `fs` between realizations. `characterize` computes decim
  once from the published `x0` and passes it via `Decim=`, otherwise the scatter
  in every rate-dependent metric would be an artefact of the protocol.
- **The embedding is a protocol, not a measurement** — see
  `quarctest.embed_policy`. Maps use delay 1 and dimension = state dimension + 1,
  because a map's iterates are its natural coordinates and AMI's first minimum on
  a map is not an embedding delay. Measured: on the logistic map, AMI's answer of
  6 with FNN's dimension gives a correlation dimension of 3.09 against a
  published D2 of exactly 1.0; at delay 1 the same estimator returns 0.94. AMI
  and FNN are still reported as metrics, next to the `embed_delay` and
  `embed_dim` actually used.

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
