# Naming conventions

One name per concept, across the whole library. A reader who learns an
argument in one function should not have to relearn it in the next.

Files and functions are `lower_snake_case`. No date suffixes — git records
when things changed.

## Canonical argument names

| concept | name | notes |
|---|---|---|
| time series | `x` | a second series is `y` |
| embedding delay | `delay` | samples between successive coordinates |
| embedding dimension | `dim` | number of coordinates |
| sampling frequency | `fs` | Hz. Use 1 for maps and iteration counts |
| recurrence radius | `radius` | absolute distance threshold |
| tolerance | `radius` | entropy tolerance; see the note below |
| maximum lag | `maxlag` | upper bound of a lag sweep |
| phase space | `Y` | an N-by-D matrix, already embedded |

Name-value options are `PascalCase`: `PhaseSpace`, `TheilerWindow`,
`Algorithm`, `Bins`, `ChunkSize`, `Norm`, `Zscore`.

## Relation to the literature

Two renames depart from notation that is standard in the source papers. The
mapping is stated in each function's header so a reader can move between the
code and the papers.

| here | literature | source |
|---|---|---|
| `dim` | `m` | Richman & Moorman (2000), sample and approximate entropy |
| `radius` | `r` | Richman & Moorman (2000), tolerance |
| `delay` | `tau`, `τ` | Fraser & Swinney (1986) and most of the embedding literature |

`tau` was consistent within this library, but it is ambiguous outside it: it
is equally conventional for time constants, and it carries no unit. `delay`
says what it is and implies samples.

`m` and `r` are single letters that mean different things in different
subfields. Since the same two quantities appear as `dim` and `radius`
throughout the rest of the library, using `m` and `r` only in the entropy
functions would preserve the paper notation at the cost of the consistency
this document exists to establish.

## Wrappers

Where several methods compute the same quantity, provide one entry point with
an `Algorithm` argument, and keep each method in its own file:

```matlab
ami(x, maxlag, Algorithm="histogram")     % ami_histogram.m
ami(x, maxlag, Algorithm="kde")           % ami_kde.m
lyapunov(x, fs, Algorithm="rosenstein")   % lye_r.m
lyapunov(x, fs, Algorithm="wolf")         % lye_w.m
```

The method belongs in an argument, not in a function name: which estimator
produced a published number is a material detail, and a name like
`ami_thomas` does not tell a reader what distinguishes it.

## Phase space as an input

Any function that reconstructs a phase space should also accept one. This
lets a reconstruction be built once and shared, and lets an algorithm be
tested apart from the embedding.

```matlab
Y = psr(x, delay, dim);
lyapunov(Y, fs, Delay=delay);
rqa(Y, delay, dim, "rec", 2.5, PhaseSpace=true);
```

A supplied phase space carries no record of the delay that built it. Any
parameter normally derived from `delay` therefore cannot be defaulted — warn
rather than guess. This is not hypothetical: in `lyapunov` a silently
defaulted Theiler window changed the exponent from 0.8022 to 0.7833 on
identical data.

## Breaking changes

The library is updated to the best interface we can define, and old spellings
are kept working through shims in `deprecated/` that warn once per session
and forward unchanged. Add that folder to your path to enable them.

Renames are announced in the README with the mapping from old to new.
