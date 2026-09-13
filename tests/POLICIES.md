# Characterization policies

What the suite computes, on which systems, and why — a readable transcription
of `quarctest.metric_policy` and `quarctest.embed_policy`. **The code is the
source of truth**; if you edit a judgement here, change it there too (each
exclusion is one `applies` flag with a `reason` string). The point of writing
the exclusions down is that a reader can disagree with a judgement they can
see — this file is where to disagree.

Every metric in the library returns a number for any series. Whether the
number *means* anything depends on the system, and a table that quietly
computed everything everywhere would put undefined quantities next to
well-defined ones with no way to tell them apart.

## Protocol constants

| constant | value | where set |
|---|---|---|
| recurrence target | 2.5% of admissible pairs | `metric_battery` `RecTarget` |
| RQA minimum line lengths | dmin = vmin = 2 (library defaults; see the dmin/vmin sweeps for what the knob does) | `rqa` defaults |
| Theiler window for `rqa_maxL_tw` | 10 samples (validated: λ-correlation −0.49 at 10, stable to 40) | `metric_battery` `TheilerWindow` |
| Wolf renormalisation interval | `evolve` = 10 samples, the library default — measured to cost accuracy on maps (see `lyapunov` help); deliberately unchanged pending a release decision | `lyapunov` default |
| AMI lag cap | min(100, N/10) | `metric_battery` `MaxLag` |

## Embedding policy (`embed_policy`)

| kind | delay | dimension | why |
|---|---|---|---|
| map | 1 | state dimension + 1 | A map's iterates are its natural coordinates; AMI's first minimum on a map is not an embedding delay (measured: logistic map at AMI's delay of 6 gives D2 = 3.09 against a published 1.0; at delay 1, 0.94) |
| flow | AMI first minimum, per realization | FNN, per realization (floored at state dim + 1) | A flow observed through one variable has a genuinely unknown embedding, estimated the way a user would estimate it |

`ami_delay` and `fnn_dim` are always *reported* for both kinds, next to the
`embed_delay`/`embed_dim` actually used, so the estimate and the choice are
both visible.

## Metrics computed everywhere

`ami_delay`, `fnn_dim`, `embed_delay`, `embed_dim` · `lyap_wolf`, `lyap_ros`
(+ six Rosenstein fit-window diagnostics) · `corr_dim` (+ six scaling-region
diagnostics) · `rqa_radius`, `rqa_det`, `rqa_lam`, `rqa_meanL`, `rqa_maxL`,
`rqa_entL`, `rqa_entV`, `rqa_entW` · `ent_samp`, `ent_ap`, `ent_permu`.

Only two carry a published reference (λ for the Lyapunov estimators, D2 for
the correlation dimension); those rows are marked `validated`. Everything
else is `characterized`: a measurement of this library's estimator on this
system, with no external truth to judge it against.

The fit diagnostics exist because R² is not a test of linearity: the Thomas
system fits its D2 region at R² = 0.988 over a visibly curved stretch and
returns 2.48 against a published 1.84. Curvature, max deviation and residual
runs are what separate straight from bent — and none of them makes an
estimate *correct* (the LCG is straight by every measure and returns
D2 = 1.95 against exactly 1).

## Exclusions — maps only

| metric | reason as enforced | confidence |
|---|---|---|
| `rqa_maxL_tw` | A map has no tangential band — iterates are not samples of a continuous trajectory — so a Theiler window would delete genuine short-time recurrences rather than sampling artifacts | solid |
| `ent_ms_s1`, `ent_ms_ci` | Multiscale coarse-graining averages adjacent samples, which for a flow approximates the same signal at lower resolution; averaging adjacent *iterates* of a map produces a series that is not the dynamics of anything | solid |
| `dfa_alpha` | Same stated reason as the multiscale entropies ("scaling exponents assume samples of an underlying continuous-time signal") | **contested — see below** |

### The DFA-on-maps question (open)

The enforced rationale treats DFA like coarse-graining: no continuum of time
scales, therefore no scaling exponent. That argument is airtight for the
multiscale entropies (their coarse-graining step literally constructs a new
series that is not the system) but weaker for DFA, on two grounds:

1. **DFA is routinely and legitimately applied to discrete event series** —
   stride intervals, RR intervals — which are iterate-indexed sequences, not
   sampled flows. DFA's profile (the cumulative sum) is well-defined for map
   iterates, and α then measures long-range correlation structure *of the
   sequence*, a meaningful quantity even when "n iterates" is not a duration.
2. **The real failure axis is regime, not kind.** On a *chaotic* map with
   fast-decaying correlations, α ≈ 0.5 — well-defined, just uninformative
   (nearly every chaotic map looks white at long lags). On a *periodic*
   orbit, F(n) oscillates and any fitted α is an artifact of the window — but
   that is equally true of a periodic flow, which the current policy happily
   computes. Regime-dependence is not map-specific.

The counterargument for keeping the exclusion: this suite's catalogues are
chaotic by construction, so map-DFA rows would be a column of ~0.5s inviting
cross-kind comparison with flow α values whose time axis means something
different; excluding them costs little information. The counter-counter: the
same "column of near-constants" description fits map-`rqa_lam` and nobody
excludes that, and a measured column of 0.5s is itself evidence the estimator
behaves.

**Status: policy currently excludes; flip it by setting `applies = true` for
`dfa_alpha` in `metric_policy.m` (drop the `mapScaleReason`) and re-running
the Sprott characterization to populate the rows.** The multiscale-entropy
exclusion should stay either way — the two arguments are genuinely different,
even though the code currently reuses one reason string for both.

## Exclusions — per system

Not policy but measurement, recorded in the reports' "incomplete ensemble"
tables: realizations that degenerate are dropped (with the survivors flagged
as a conditioned sample), and systems that defeat a step entirely carry n=0
rows with the mechanism named — e.g. PanXuZhou cannot reach the 2.5%
recurrence target on any realization, so every RQA row for it is empty by
measurement, not by judgement.
