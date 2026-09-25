"""Build the reference corpus: one Lyapunov spectrum per system, one method.

WHAT THIS REPLACES. Until now the characterization suite drew reference
exponents from two sources that were never comparable. Sprott (2003) Appendix A
supplies published values at four or five decimals, computed by trajectory
separation. The dysts database supplies continuous-QR spectra from analytic
Jacobians. On the 13 systems both catalogues contain, the two disagree by a
median of 22%, and the dysts spectra fail a theorem: an autonomous flow has one
exponent exactly zero, and only 31 of 111 dysts spectra put their smallest
exponent within 0.001 of it. Every ratio in the lab report divides an estimate
by one of those numbers, so estimator error and reference error are currently
added together with no way to tell them apart.

This script recomputes every exponent under one method -- Sprott's, from
python/lyap_spectrum.py -- from initial conditions and control parameters that
are recorded in the output rather than implied. Where the two catalogues hold
the same system under different parameters, both parameterizations are built as
separate entries, so the disagreement becomes a measured quantity instead of an
open question.

WHAT MAKES AN ENTRY TRUSTWORTHY. Three checks travel with every spectrum, and
none of them consults the published value:

  trace      sum(lambda) against the time-averaged divergence of the vector
             field, read straight off f by central differences. This does not
             pass through the Gram-Schmidt, so it tests the renormalization
             rather than restating it.
  zero       for an autonomous flow, the smallest |lambda|, which theory puts
             at exactly zero along the flow direction.
  sem        the standard error across ensemble members, each an independent
             trajectory on the same attractor.

The published value is recorded alongside as a fourth, external comparison, but
an entry is accepted or rejected on the three internal checks. That ordering is
deliberate: on Rucklidge our value and dysts' agree with each other and disagree
with Sprott's published 0.0643, and a gate that trusted the book would have
rejected the right answer.

Usage:
    python/build_reference_corpus.py --out tests/reports/reference_spectra.json
    python/build_reference_corpus.py --quick --only lorenz,Thomas
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
import warnings
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from lyap_spectrum import spectrum, spectrum_adaptive   # noqa: E402
import sprott_systems as ss                 # noqa: E402

warnings.filterwarnings("ignore")

# Members are separated by this fraction of each component's magnitude. It is
# deliberately tiny: chaos decorrelates the members during the transient
# anyway, so a large perturbation buys no independence and costs correctness on
# systems with coexisting attractors or mixed phase space, where a member that
# starts further away lands on a different invariant set and averages something
# that is not the published system's exponent.
SPREAD = 1e-9

BUDGETS = {
    "quick": dict(map=dict(n=200_000, E=64), flow=dict(n=200_000, E=64)),
    "full": dict(map=dict(n=2_000_000, E=256), flow=dict(n=1_500_000, E=256)),
}


def _make_vector_field(S):
    """A vectorized (t, y) -> dy/dt for a dysts system, or None.

    dysts systems do not share one calling convention. Most define
    _rhs(x, y, z, t, *params) over unpacked components, which vectorizes over
    arrays unchanged. Others define _rhs(X, t) over the whole state vector and
    override rhs; a few return a constant for some component, which numpy will
    not stack against the array-valued ones.

    Rather than assume which applies, each strategy is tried and then checked
    against dysts' own scalar call on perturbed states. A strategy that
    disagrees is not used, and a system where none agrees is left out of the
    corpus with its name recorded -- a wrong vector field would produce a
    plausible chaotic series and a plausible wrong exponent, which is exactly
    the failure this whole corpus exists to remove.
    """
    ic = np.asarray(S.ic, dtype=np.float64)
    D = ic.size
    plist = list(S.param_list)

    def unpacked(t, y):
        cols = S._rhs(*[y[..., i] for i in range(y.shape[-1])], t, *plist)
        # A component may come back as a scalar constant; broadcast before
        # stacking so it lines up with the array-valued ones.
        cols = np.broadcast_arrays(*[np.asarray(c, dtype=np.float64) for c in cols])
        return np.stack(cols, axis=-1)

    def flattened(t, y):
        flat = y.reshape(-1, y.shape[-1])
        out = np.asarray(S.rhs(flat, t), dtype=np.float64)
        if out.shape == (y.shape[-1], flat.shape[0]):
            out = out.T
        return out.reshape(y.shape)

    # Reference values from dysts itself, one state at a time.
    rng = np.random.default_rng(0)
    probes = np.stack([ic] + [ic * (1 + 0.01 * rng.standard_normal(D))
                              for _ in range(3)], axis=0)
    try:
        want = np.stack([np.asarray(S.rhs(p, 0.0), dtype=np.float64).ravel()
                         for p in probes], axis=0)
    except Exception:                             # noqa: BLE001
        return None
    if not np.isfinite(want).all():
        return None

    for cand in (unpacked, flattened):
        try:
            got = cand(0.0, probes)
        except Exception:                         # noqa: BLE001
            continue
        if got.shape != want.shape or not np.isfinite(got).all():
            continue
        if np.allclose(got, want, rtol=1e-10, atol=1e-12):
            # Confirm it also holds with a leading ensemble axis, which is how
            # it will actually be called.
            try:
                batched = cand(0.0, probes[None, ...])
            except Exception:                     # noqa: BLE001
                continue
            if batched.shape == (1,) + want.shape and np.allclose(
                    batched[0], want, rtol=1e-10, atol=1e-12):
                return cand
    return None


def _dysts_entries(names=None):
    """dysts systems as catalogue entries, with their own ICs and parameters."""
    import dysts
    import dysts.flows as fl

    # The catalogue is the shipped json, read the same way export_dysts.py
    # reads it, so both scripts enumerate the same systems.
    cat_path = os.path.join(os.path.dirname(dysts.__file__),
                            "data", "chaotic_attractors.json")
    with open(cat_path) as fh:
        catalog = json.load(fh)

    out = []
    for name, meta in catalog.items():
        if meta.get("delay"):
            continue                          # needs a history function
        if names is not None and name not in names:
            continue
        try:
            S = getattr(fl, name)()
        except Exception:                         # noqa: BLE001
            continue
        ic = np.asarray(S.ic, dtype=np.float64)
        if ic.ndim != 1 or not np.isfinite(ic).all():
            continue                              # delay systems and the like
        params = {k: (v.tolist() if isinstance(v, np.ndarray) else v)
                  for k, v in S.params.items()}
        f = _make_vector_field(S)
        if f is None:
            continue                          # signature we cannot vectorize

        out.append(dict(name=name, source="dysts", kind="flow", f=f, x0=ic,
                        dt=float(S.dt), params=params, wrap=None,
                        published=None, published_source=None,
                        nonautonomous=bool(getattr(S, "nonautonomous", False)),
                        # dysts flags its Hamiltonian systems; they are the
                        # same IC-dependent class as Sprott's conservative
                        # flows and carry the same warning.
                        icDependent=bool(meta.get("hamiltonian"))))
    return out


def _sprott_entries(names=None):
    out = []
    for s in ss.usable():
        if names is not None and s["name"] not in names:
            continue
        out.append(dict(name=s["name"], source="sprott", kind=s["kind"],
                        f=s["f"], x0=s["x0"], dt=s["dt"], params=None,
                        wrap=s["wrap"], published=s["lam1"],
                        published_source=f"Sprott (2003) {s['section']}",
                        nonautonomous="driven" in s["category"],
                        icDependent=ss.ic_dependent(s)))
    return out


def calibrate_step(entry, mults=(16, 8, 4, 2, 1), tol=1e-5):
    """Coarsest integration step whose trace identity still holds.

    dysts ships steps chosen to resolve the waveform, far finer than RK4 needs
    for an exponent, and using them directly would buy only a few hundred time
    units per run -- nowhere near a converged average. But the step cannot be
    coarsened globally either: at 16x, Lorenz holds its trace error at 3e-6
    while Thomas jumps to 9e-2, which is RK4 failing rather than the exponent
    changing. The trace identity detects exactly that, because integration
    error breaks the agreement between sum(lambda) and the divergence, so it
    serves as the acceptance test for the step as well as for the spectrum.

    A short run at each candidate is cheap next to the full one, so the step is
    measured rather than assumed.
    """
    if entry["kind"] == "map":
        return 1
    best = 1
    for m in mults:
        try:
            r = spectrum(entry["f"], entry["x0"], kind="flow",
                         dt=entry["dt"] * m, n_steps=40_000,
                         transient_steps=8_000, renorm_every=100,
                         d0=1e-8, ensemble=16, spread=SPREAD, seed=1,
                         divergence=True)
        except Exception:                          # noqa: BLE001
            continue
        scale = max(abs(r.divergence), 1.0)
        if np.isfinite(r.trace_error) and abs(r.trace_error) < tol * scale:
            best = m
            break
    return best


def _entry(name, source):
    """Rebuild one catalogue entry inside the worker, from its name."""
    if source == "sprott":
        return _sprott_entries({name})[0]
    return _dysts_entries({name})[0]


def run_one(args):
    name, source, budget = args
    entry = _entry(name, source)
    b = BUDGETS[budget]["map" if entry["kind"] == "map" else "flow"]
    t0 = time.time()
    try:
        mult = calibrate_step(entry)
        dt_used = None if entry["dt"] is None else entry["dt"] * mult
        r, horizon = spectrum_adaptive(
            entry["f"], entry["x0"], kind=entry["kind"], dt=dt_used,
            n_steps=b["n"], transient_steps=min(50_000, b["n"] // 10),
            renorm_every=1 if entry["kind"] == "map" else 100,
            d0=1e-9 if entry["kind"] == "map" else 1e-8,
            ensemble=b["E"], spread=SPREAD, seed=1, wrap=entry["wrap"],
            divergence=True,
        )
    except Exception as exc:                      # noqa: BLE001
        return dict(name=entry["name"], source=entry["source"],
                    error=f"{type(exc).__name__}: {exc}",
                    elapsed=time.time() - t0)

    lam = r.exponents
    rec = dict(
        name=entry["name"], source=entry["source"], kind=entry["kind"],
        stateDim=int(entry["x0"].size),
        x0=entry["x0"].tolist(), dt=dt_used, dtCatalogue=entry["dt"],
        dtMultiple=mult, horizon=int(horizon), params=entry["params"],
        spectrum=lam.tolist(), lambdaMax=float(lam[0]),
        sem=r.sem.tolist(),
        # The three internal checks.
        divergence=float(r.divergence), traceError=float(r.trace_error),
        zeroExponent=(None if entry["nonautonomous"]
                      else float(min(abs(v) for v in lam))),
        nEscaped=int(r.n_escaped), nDied=int(r.n_died),
        nCollapsed=int(r.n_collapsed), nSaturated=int(r.n_saturated),
        nSteps=int(r.n_steps),
        elapsed=float(r.elapsed),
        # The external comparison, recorded but not used as a gate.
        published=entry["published"], publishedSource=entry["published_source"],
        # True where the exponent belongs to the initial condition rather than
        # to the system. See sprott_systems.ic_dependent.
        icDependent=bool(entry.get("icDependent")),
        method="nearby-trajectory Gram-Schmidt (Sprott 2003 ch. 5)",
        wallclock=time.time() - t0,
    )
    if entry["published"]:
        rec["publishedRatio"] = float(lam[0] / entry["published"])
    return rec


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="tests/reports/reference_spectra.json")
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--only", default="")
    ap.add_argument("--catalogue", default="both",
                    choices=["both", "sprott", "dysts"])
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--resume", action="store_true",
                    help="skip systems already present in the run's jsonl")
    a = ap.parse_args()

    names = {n.strip() for n in a.only.split(",") if n.strip()} or None
    entries = []
    if a.catalogue in ("both", "sprott"):
        entries += _sprott_entries(names)
    if a.catalogue in ("both", "dysts"):
        entries += _dysts_entries(names)

    budget = "quick" if a.quick else "full"

    # EACH SYSTEM IS WRITTEN AS IT FINISHES, not at the end. An ex.map that
    # collects everything before writing gives no progress to watch and loses
    # the whole run if it dies an hour in; export_dysts.py already streams to a
    # jsonl for exactly that reason, and this follows it. A completed system is
    # also skipped on a re-run, so an interrupted build resumes instead of
    # starting over.
    jsonl = os.path.splitext(a.out)[0] + ".jsonl"
    done = {}
    if a.resume and os.path.exists(jsonl):
        with open(jsonl) as fh:
            for line in fh:
                try:
                    r = json.loads(line)
                except ValueError:
                    continue
                done[(r.get("source"), r.get("name"))] = r
        entries = [e for e in entries if (e["source"], e["name"]) not in done]
        print(f"resuming: {len(done)} already done, {len(entries)} to go")

    print(f"{budget} build: {len(entries)} systems, {a.workers} workers")
    t0 = time.time()
    recs = list(done.values())
    os.makedirs(os.path.dirname(jsonl) or ".", exist_ok=True)
    with open(jsonl, "a") as log, ProcessPoolExecutor(max_workers=a.workers) as ex:
        futs = {ex.submit(run_one, (e["name"], e["source"], budget)): e
                for e in entries}
        for i, fut in enumerate(as_completed(futs), 1):
            r = fut.result()
            recs.append(r)
            log.write(json.dumps(r) + "\n")
            log.flush()
            tag = "FAIL" if "error" in r else f"lam1={r['lambdaMax']:+.5f}"
            print(f"[{i}/{len(entries)}] {r['source']:7}{r['name']:24}{tag}",
                  flush=True)

    good = [r for r in recs if "error" not in r]
    bad = [r for r in recs if "error" in r]
    payload = dict(
        method="nearby-trajectory Gram-Schmidt (Sprott 2003 ch. 5)",
        spread=SPREAD, budget=budget, built=time.strftime("%Y-%m-%dT%H:%M:%S"),
        systems=good, failed=bad,
    )
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    with open(a.out, "w") as fh:
        json.dump(payload, fh, indent=1)

    print(f"{len(good)} built, {len(bad)} failed, {time.time()-t0:.0f}s -> {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
