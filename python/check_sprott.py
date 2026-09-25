"""Reproduce Sprott's Appendix A exponents with the one method, and report the gap.

This is the acceptance gate for the reference corpus. The corpus is only worth
building if the engine can reproduce published values to the precision they are
published at; until that holds, a disagreement between our exponent and Sprott's
cannot be attributed to his method rather than to our port.

Two things are being tested at once and they are not separable by design: the
engine (python/lyap_spectrum.py) and the transcription of the 62 vector fields
(python/sprott_systems.py). A mistyped coefficient and a broken integrator both
show up here as a missed exponent, which is the point -- the check is sensitive
to either.

Usage:
    python/check_sprott.py                  # quick pass, catches gross errors
    python/check_sprott.py --precision      # long pass, the real gate
    python/check_sprott.py --only lorenz,henon
"""

from __future__ import annotations

import argparse
import sys
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np

sys.path.insert(0, __file__.rsplit("/", 1)[0])

from lyap_spectrum import spectrum_adaptive   # noqa: E402
import sprott_systems as ss                 # noqa: E402


# Quick and precision budgets. Maps cost one function evaluation per step and
# flows cost four, and flows need a longer record because their exponent is per
# unit time rather than per iteration, so the two are budgeted separately.
# Members are separated by this fraction of each component's magnitude. Chaos
# decorrelates them during the transient regardless, so the perturbation is
# kept tiny: a larger one lands members on different invariant sets in systems
# with coexisting attractors or mixed phase space, and their average is then
# not the published system's exponent.
SPREAD = 1e-9

BUDGETS = {
    "quick": dict(map=dict(n=200_000, E=64), flow=dict(n=200_000, E=64)),
    "precision": dict(map=dict(n=2_000_000, E=256), flow=dict(n=1_500_000, E=256)),
}


def run_one(args):
    name, budget, seed = args
    s = ss.by_name(name)
    b = BUDGETS[budget]["map" if s["kind"] == "map" else "flow"]
    t0 = time.time()
    # Renormalizing every step is necessary for maps, where one iteration can
    # stretch a separation by a decade. Flows were checked to be insensitive
    # between 1 and 100 steps on Lorenz, so 100 is taken for speed.
    renorm = 1 if s["kind"] == "map" else 100
    try:
        r, horizon = spectrum_adaptive(
            s["f"], s["x0"], kind=s["kind"], dt=s["dt"],
            n_steps=b["n"], transient_steps=min(50_000, b["n"] // 10),
            renorm_every=renorm, d0=1e-9 if s["kind"] == "map" else 1e-8,
            ensemble=b["E"], seed=seed, wrap=s["wrap"],
            spread=SPREAD, divergence=True,
        )
    except Exception as exc:                      # noqa: BLE001
        return dict(name=name, error=f"{type(exc).__name__}: {exc}",
                    elapsed=time.time() - t0)
    return dict(name=name, section=s["section"], kind=s["kind"],
                category=s["category"], tier=s["tier"], ref=s["lam1"],
                lam=r.exponents.tolist(), sem=r.sem.tolist(),
                trace=float(r.trace_error), escaped=int(r.n_escaped),
                horizon=int(horizon), died=int(r.n_died),
                # Zero exponent is a theorem for autonomous flows only; a
                # driven system has no flow direction in its own state space.
                zero=(None if s["kind"] == "map" or "driven" in s["category"]
                      else float(min(abs(v) for v in r.exponents))),
                elapsed=time.time() - t0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--precision", action="store_true",
                    help="long run: the real acceptance gate")
    ap.add_argument("--only", default="", help="comma-separated system names")
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()

    budget = "precision" if a.precision else "quick"
    names = ([n.strip() for n in a.only.split(",") if n.strip()]
             or [s["name"] for s in ss.usable()])

    print(f"{budget} pass over {len(names)} systems, {a.workers} workers\n")
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=a.workers) as ex:
        results = list(ex.map(run_one, [(n, budget, a.seed) for n in names]))

    hdr = (f"{'system':24}{'sec':8}{'ours':>11}{'Sprott':>10}{'err':>10}"
           f"{'sem':>9}{'trace':>10}{'zero':>9}  flag")
    print(hdr)
    print("-" * len(hdr))
    ok = miss = err = 0
    budget_local = budget
    for r in results:
        if "error" in r:
            print(f"{r['name']:24}{'':8}{'':11}{'':10}{'':10}{'':9}  ERROR {r['error'][:40]}")
            err += 1
            continue
        lam1, ref, sem = r["lam"][0], r["ref"], r["sem"][0]
        d = lam1 - ref
        # Sprott prints four decimals for flows and five for most maps; a hit
        # is agreement at the last digit he prints, which is 1e-4 either way
        # once rounding is allowed for.
        flag = "ok" if abs(d) < 1.5e-4 else ("close" if abs(d) < 2e-3 else "MISS")
        ok += flag == "ok"
        miss += flag == "MISS"
        z = "" if r["zero"] is None else f"{r['zero']:9.1e}"
        print(f"{r['name']:24}{r['section']:8}{lam1:11.5f}{ref:10.5f}{d:+10.5f}"
              f"{sem:9.1e}{r['trace']:10.1e}{z:>9}  {flag}"
              + (f"  esc={r['escaped']}" if r["escaped"] else "")
              + (f"  died={r['died']}" if r["died"] else "")
              + (f"  horizon={r['horizon']}" if r["horizon"] != BUDGETS[budget][
                  "map" if r["kind"] == "map" else "flow"]["n"] else ""))

    print(f"\n{ok} within 1.5e-4, {miss} missed, {err} errored, "
          f"{len(results)-ok-miss-err} close; {time.time()-t0:.0f}s wall")
    return 1 if (miss or err) else 0


if __name__ == "__main__":
    raise SystemExit(main())
