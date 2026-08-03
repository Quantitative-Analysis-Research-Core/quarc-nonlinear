# Sprott (2003) Appendix A, reference values

Transcribed from J. C. Sprott, *Chaos and Time-Series Analysis*, Oxford
University Press (2003), Appendix A, "Common chaotic systems".

This is the human-readable record of the values encoded in
`tests/matlab/+nonantest/sprott_catalog.m`, and is generated from that
catalogue so the two cannot drift. Regenerate with:

```
matlab -batch "addpath('tests/matlab'); nonantest.write_sprott_markdown"
```

## What the appendix says about these numbers

> All Lyapunov exponents are base-e and were calculated using the methods
> in Chapter 5. The Kaplan-Yorke dimension is given by eqn (5.29). The
> correlation dimension uses at least 2 x 10^12 pairs, corresponding to
> data sets of over two million points for each case, and have been
> extrapolated to the limit of zero size scale.

and: "Most of the results are original calculations, and all have been
independently verified."

Two consequences for this library:

- lambda is in **nats** per unit time. `lye_w` accumulates log2 and returns
  bits, so multiply its output by ln 2 before comparing.
- D2 carries real uncertainty and is not a tight target. Values shown
  without a tolerance are given as exact in the appendix.

## Systems

| section | name | category | lambda (nats) | tier | D2 | usable |
|---|---|---|---|---|---|---|
| A.1.1 | `logistic` | noninvertible map | 0.6931 | exact | 1.000 (exact) | true |
| A.1.2 | `sine_map` | noninvertible map | 0.6891 | numerical | 1.000 (exact) | true |
| A.1.3 | `tent` | noninvertible map | 0.6931 | exact | 1.000 (exact) | false |
| A.1.4 | `lcg` | noninvertible map | 8.8736 | exact | 1.000 (exact) | true |
| A.1.5 | `cubic_map` | noninvertible map | 1.0986 | numerical | 1.000 (exact) | true |
| A.1.6 | `ricker` | noninvertible map | 0.3848 | numerical | 1.000 (exact) | true |
| A.1.7 | `gauss_map` | noninvertible map | 2.3734 | numerical | 1.000 (exact) | true |
| A.1.8 | `cusp` | noninvertible map | 0.5000 | exact | 1.000 (exact) | true |
| A.1.9 | `gauss_white` | noninvertible map | 0.6931 | exact | 1.000 (exact) | false |
| A.1.10 | `pinchers` | noninvertible map | 0.4679 | numerical | 1.000 (exact) | true |
| A.1.11 | `spence` | noninvertible map | infinite | exact | 1.000 (exact) | false |
| A.1.12 | `sine_circle` | noninvertible map | 0.3539 | numerical | 1.000 (exact) | true |
| A.2.1 | `henon` | dissipative map | 0.4192 | numerical | 1.220 +/- 0.036 | true |
| A.2.2 | `lozi` | dissipative map | 0.4702 | numerical | 1.384 +/- 0.053 | true |
| A.2.3 | `delayed_logistic` | dissipative map | 0.1831 | numerical | 1.144 +/- 0.034 | true |
| A.2.4 | `tinkerbell` | dissipative map | 0.1900 | numerical | 1.329 +/- 0.036 | true |
| A.2.5 | `burgers` | dissipative map | 0.1208 | numerical | 1.462 +/- 0.054 | true |
| A.2.6 | `holmes_cubic` | dissipative map | 0.5946 | numerical | 1.260 +/- 0.039 | true |
| A.2.7 | `kaplan_yorke` | dissipative map | 0.6931 | exact | 1.432 +/- 0.044 | false |
| A.2.8 | `dissipative_standard` | dissipative map | 1.4700 | numerical | 1.356 +/- 0.047 | true |
| A.2.9 | `ikeda` | dissipative map | 0.5076 | numerical | 1.690 +/- 0.073 | true |
| A.2.10 | `sinai` | dissipative map | 0.9595 | numerical | 1.779 +/- 0.063 | true |
| A.2.11 | `predator_prey` | dissipative map | 0.1966 | numerical | 1.903 +/- 0.079 | true |
| A.3.1 | `chirikov` | conservative map | 0.1050 | numerical | 1.954 +/- 0.077 | true |
| A.3.2 | `henon_area` | conservative map | 0.0064 | numerical | 2.200 +/- 0.063 | true |
| A.3.3 | `arnold_cat` | conservative map | 0.9624 | exact | 2.000 (exact) | true |
| A.3.4 | `gingerbreadman` | conservative map | 0.0734 | numerical | 2.171 +/- 0.078 | true |
| A.3.5 | `chaotic_web` | conservative map | 0.0485 | numerical | 1.779 +/- 0.059 | true |
| A.3.6 | `lorenz3d_map` | conservative map | 0.0746 | numerical | 1.745 +/- 0.057 | true |
| A.4.1 | `damped_pendulum` | driven flow | 0.1414 | numerical | 2.764 +/- 0.158 | true |
| A.4.2 | `driven_vdp` | driven flow | 0.1933 | numerical | 2.190 +/- 0.080 | true |
| A.4.3 | `shaw_vdp` | driven flow | 0.1180 | numerical | 2.007 +/- 0.091 | true |
| A.4.4 | `brusselator` | driven flow | 0.0140 | numerical | 2.224 +/- 0.095 | true |
| A.4.5 | `ueda` | driven flow | 0.1034 | numerical | 2.675 +/- 0.132 | true |
| A.4.6 | `duffing_two_well` | driven flow | 0.1572 | numerical | 2.334 +/- 0.114 | true |
| A.4.7 | `duffing_vdp` | driven flow | 0.0963 | numerical | 2.333 +/- 0.115 | true |
| A.4.8 | `rayleigh_duffing` | driven flow | 0.0912 | numerical | 2.194 +/- 0.120 | true |
| A.5.1 | `lorenz` | autonomous flow | 0.9056 | numerical | 2.068 +/- 0.086 | true |
| A.5.2 | `rossler` | autonomous flow | 0.0714 | numerical | 1.991 +/- 0.065 | true |
| A.5.3 | `diffusionless_lorenz` | autonomous flow | 0.2101 | numerical | 2.169 +/- 0.128 | true |
| A.5.4 | `complex_butterfly` | autonomous flow | 0.1690 | numerical | 2.491 +/- 0.131 | true |
| A.5.5 | `chen` | autonomous flow | 2.0272 | numerical | 2.147 +/- 0.117 | true |
| A.5.6 | `hadley` | autonomous flow | 0.1665 | numerical | 2.162 +/- 0.114 | true |
| A.5.7 | `act` | autonomous flow | 0.1634 | numerical | 2.039 +/- 0.106 | true |
| A.5.8 | `rabinovich_fabrikant` | autonomous flow | 0.1981 | numerical | 2.191 +/- 0.113 | true |
| A.5.9 | `rigid_body` | autonomous flow | 0.1421 | numerical | 2.069 +/- 0.121 | true |
| A.5.10 | `chua` | autonomous flow | 0.3271 | numerical | 2.125 +/- 0.098 | true |
| A.5.11 | `moore_spiegel` | autonomous flow | 0.1119 | numerical | 2.309 +/- 0.107 | true |
| A.5.12 | `thomas` | autonomous flow | 0.0349 | numerical | 1.843 +/- 0.075 | true |
| A.5.13 | `halvorsen` | autonomous flow | 0.7899 | numerical | 2.110 +/- 0.095 | true |
| A.5.14 | `burke_shaw` | autonomous flow | 2.2499 | numerical | 2.211 +/- 0.132 | true |
| A.5.15 | `rucklidge` | autonomous flow | 0.0643 | numerical | 2.108 +/- 0.095 | true |
| A.5.16 | `windmi` | autonomous flow | 0.0755 | numerical | 2.035 +/- 0.095 | true |
| A.5.17 | `simplest_quadratic` | autonomous flow | 0.0551 | numerical | 2.187 +/- 0.075 | true |
| A.5.18 | `simplest_cubic` | autonomous flow | 0.0837 | numerical | 2.174 +/- 0.083 | true |
| A.5.19 | `simplest_piecewise` | autonomous flow | 0.0362 | numerical | 2.131 +/- 0.072 | true |
| A.5.20 | `double_scroll` | autonomous flow | 0.0497 | numerical | 2.184 +/- 0.107 | true |
| A.6.1 | `driven_pendulum` | conservative flow | 0.1633 | numerical | 2.756 +/- 0.149 | true |
| A.6.2 | `simplest_driven` | conservative flow | 0.0971 | numerical | 2.634 +/- 0.160 | true |
| A.6.3 | `nose_hoover` | conservative flow | 0.0138 | numerical | 2.521 +/- 0.146 | true |
| A.6.4 | `labyrinth` | conservative flow | 0.1402 | numerical | 2.837 +/- 0.173 | true |
| A.6.5 | `henon_heiles` | conservative flow | 0.0450 | numerical | 2.706 +/- 0.126 | true |

## Systems excluded from benchmarking

**`tent`** (A.1.3)

The slope-2 tent map is exactly a binary shift, so in double precision it consumes one mantissa bit per iteration and reaches exactly 0 after ~52 steps: 52 distinct values in 2000 from Sprott's own x0 = 1/sqrt(2). Sprott specifies an irrational x0 for precisely this reason, but EVERY double is a dyadic rational, so no representable initial condition avoids it. Not an implementation fault and not fixable in double precision; it needs exact or extended arithmetic. Sprott flags the issue at his section 2.5.2.

**`gauss_white`** (A.1.9)

Transcribed literally from the appendix as A*erfinv(1 - 2*erf(x/A)), this leaves the domain on step 2: once x goes negative, erf(x/A) < 0, so 1 - 2*erf(x/A) > 1 and erfinv returns NaN. Marked unusable rather than guessed at, because the intended form could not be verified from the appendix page alone. Note also that lambda = ln 2 indicates conjugacy to the doubling map, so even a correct transcription would inherit the same mantissa exhaustion that disqualifies the tent map above.

**`spence`** (A.1.11)

lambda -> Inf. No finite reference exists, so it cannot be used as a benchmark. Retained so the catalogue is complete and the exclusion is explicit rather than silent.

**`kaplan_yorke`** (A.2.7)

The x equation is x -> 2x mod 1, the same binary shift as the tent map, so x reaches exactly 0 after ~51 iterations and the attractor collapses. Same root cause, same non-fixability in double precision; Sprott again refers to his section 2.5.2.

## Counts

| category | systems |
|---|---|
| noninvertible map | 12 |
| dissipative map | 11 |
| conservative map | 6 |
| driven flow | 8 |
| autonomous flow | 20 |
| conservative flow | 5 |
| **total** | **62** |

58 of 62 are usable as benchmarks. 8 have lambda in closed form, marked
tier `exact`; the rest are well-converged numerical results.
