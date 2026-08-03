# Porting to other languages

The MATLAB library is the reference implementation. Ports to Python and R are
planned, and the standard they are held to is **identical function and
performance** — the same numbers, at comparable speed, not merely the same
algorithm described in the same words.

This document records how that is enforced, so a port can be checked rather
than assumed correct.

## MATLAB is the spec

Function names, signatures, argument conventions, default values and returned
quantities in `matlab/` define the interface. A port should match them unless
there is a language reason not to, and where it deviates it should say so.

That means changes to the MATLAB library are changes to the spec. Renaming an
argument or changing a default is not a local decision once ports exist.

## How agreement is checked

`tests/fixtures/` holds the mechanism:

| file | role |
|---|---|
| `*.csv` | input series, written at full float64 precision |
| `matlab_reference.json` | MATLAB's answers on those exact series |
| `make_fixtures.py` | regenerates the series (run deliberately, not per test) |
| `lye_benchmark_results.csv` | Lyapunov estimators over Sprott (2003) Appendix A |

Every language reads the **same samples** and is compared against the **same
frozen reference**. This is what makes an assertion like "agrees to 1e-10"
meaningful.

It is necessary rather than convenient: MATLAB, NumPy and R cannot draw the
same random numbers. A test that generates its own series in each language
compares two different realisations and can only check loose statistical
agreement, which is far too weak to catch a port that has drifted. Comparing
estimators over identical input isolates estimator disagreement from RNG
disagreement.

Regenerate the reference after an intentional change to a MATLAB result:

```bash
matlab -batch "addpath('tests/matlab'); make_reference"
```

Fixtures and references are **inputs** to the tests, never outputs. A suite
that can rewrite its own expected values will ratify a regression.

## What "identical" is allowed to mean

Some divergence is legitimate and should be stated rather than hidden:

- **Floating-point summation order.** A streaming accumulation cannot
  reproduce the exact order of a batch `mean`. Differences at 1e-15 are
  expected; differences at 1e-4 are a defect.
- **Library primitives.** `std` in MATLAB uses the N-1 denominator; NumPy's
  `np.std` defaults to N. This is a real source of divergence and did bite
  `ent_samp`, where it combined with a template-count difference to produce a
  1.9e-4 disagreement in which the two errors partly cancelled — correcting
  either one alone made the gap wider.

  Both are now resolved. `ent_samp` follows Richman & Moorman, counting A and
  B over the same N-dim templates, and a port matches the MATLAB reference
  exactly once it passes `ddof=1` to `np.std`. Measured on the shared fixture:

  | | SampEn |
  |---|---|
  | MATLAB reference | 1.7961940663 |
  | port with `ddof=1` | 1.7961940663 (exact) |
  | port with `ddof=0` | 1.7957514905 (4.4e-4 out) |
- **Performance.** Comparable, not equal. An order of magnitude apart on the
  same input is worth investigating.

## Before writing a port

Two lessons from auditing the existing Python port, recorded so they are not
repeated:

1. **Never swallow an exception.** `emd.py` wrapped its sifting loop in a bare
   `except:` and returned the input unchanged when anything went wrong, with
   no warning and a plausible-looking return value. A crash is recoverable; a
   wrong answer that looks right is not.

2. **Return the same shape regardless of the data.** The same function
   returned `(2, 512)` on one series and `(512,)` on another that differed by
   a single sample. A caller cannot defend against that.

Both have MATLAB analogues that were also fixed: `LyE_R` returned two columns
or three depending on the data.

## Status

- **MATLAB** — reference implementation, active.
- **Python** — port exists in the deprecated NONAN repository, has known
  defects, not moved here. Needs rework rather than a move.
- **R** — planned.
