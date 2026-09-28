"""Audit the reference corpus, and measure how far the published values were off.

Two jobs, in this order.

FIRST, is the corpus itself sound? Every entry carries three checks that never
consult a published number: the trace identity, the zero exponent for an
autonomous flow, and the standard error across ensemble members. Entries that
fail them are named here rather than quietly averaged into a table. This has to
come first, because the second job assumes the new values are trustworthy.

SECOND, how much of the lab report's ratio error was never the estimator's? The
report divides each estimate by a published reference, so a reference that is
itself 20% high inflates every ratio built on it by the same 20%. Comparing the
recomputed spectra against the values they replace turns that from a worry into
a number, per system and per catalogue.

Usage:
    python/audit_corpus.py                      # audit + comparison
    python/audit_corpus.py --json summary.json  # also write the summary
"""

from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

# Bands for the internal checks. The trace tolerance is relative to the
# divergence, because a system whose volumes contract at 14 per unit time
# cannot be held to the same absolute error as one contracting at 0.2.
TRACE_REL = 1e-4
TRACE_ABS = 1e-6
ZERO_ABS = 1e-3
SEM_REL = 1e-2


def band(rec):
    """Which internal checks an entry passes, and why it fails if it does."""
    why = []
    lam1 = rec["lambdaMax"]
    tr, div = rec.get("traceError"), rec.get("divergence")
    if tr is None or not np.isfinite(tr):
        why.append("trace:missing")
    else:
        tol = max(TRACE_ABS, TRACE_REL * max(abs(div or 0.0), 1.0))
        # A map with a critical point or a discontinuity defeats the
        # finite-difference determinant -- log|f'| diverges where f' vanishes,
        # and an integer map has no meaningful derivative at all -- so the
        # trace check is reported for maps but not enforced on them.
        if abs(tr) > tol and rec["kind"] != "map":
            why.append(f"trace:{tr:+.1e}")
    # The zero exponent is a theorem about autonomous FLOWS: it is the
    # direction along the trajectory. A map has no such direction, so the
    # field is not a check on one.
    z = rec.get("zeroExponent") if rec["kind"] == "flow" else None
    if z is not None and z > ZERO_ABS:
        why.append(f"zero:{z:.1e}")
    sem = rec["sem"][0]
    if lam1 > 0 and sem > SEM_REL * abs(lam1):
        why.append(f"sem:{sem/abs(lam1)*100:.1f}%")
    for k, label in (("nDied", "died"), ("nEscaped", "escaped"),
                     ("nSaturated", "saturated")):
        if rec.get(k):
            why.append(f"{label}:{rec[k]}")
    return why


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--corpus",
                    default=os.path.join(ROOT, "tests/reports/reference_spectra.json"))
    ap.add_argument("--manifest",
                    default=os.path.join(ROOT, "tests/reports/dysts/manifest.json"))
    ap.add_argument("--json", default="")
    a = ap.parse_args()

    with open(a.corpus) as fh:
        corpus = json.load(fh)
    recs = corpus["systems"]
    print(f"corpus: {len(recs)} systems, {len(corpus.get('failed', []))} failed, "
          f"method = {corpus['method']}\n")

    # ---- 1. internal audit
    flagged = [(r, band(r)) for r in recs]
    bad = [(r, w) for r, w in flagged if w]
    print(f"INTERNAL CHECKS: {len(recs) - len(bad)} of {len(recs)} clean")
    if bad:
        print(f"{'system':26}{'src':8}{'lambda1':>10}   flags")
        for r, w in sorted(bad, key=lambda x: -len(x[1]))[:25]:
            print(f"{r['name']:26}{r['source']:8}{r['lambdaMax']:10.5f}   {', '.join(w)}")
    print()

    n_icdep = sum(1 for r in recs if r.get("icDependent"))
    n_cons = sum(1 for r in recs if r.get("conservative"))
    notch = [r["name"] for r in recs if r.get("chaoticAtIC") is False]
    print(f"REFERENCE FLAGS: {n_icdep} icDependent ({n_cons} by divergence), "
          f"{len(notch)} not chaotic at their initial condition: {', '.join(notch)}\n")

    # ---- 2. against the published values. A system whose recorded initial
    # condition is not chaotic has a reference of essentially zero, and a ratio
    # against it says nothing about either value; those are left out here.
    scorable = [r for r in recs if r.get("chaoticAtIC", True)]
    sp = [r for r in scorable if r["source"] == "sprott" and r.get("published")]
    if sp:
        ratio = np.array([r["publishedRatio"] for r in sp])
        err = np.array([r["lambdaMax"] - r["published"] for r in sp])
        print(f"vs SPROTT published ({len(sp)} systems)")
        print(f"  |error|: median {np.median(np.abs(err)):.2e}, "
              f"90th pct {np.percentile(np.abs(err), 90):.2e}")
        print(f"  within 1.5e-4: {(np.abs(err) < 1.5e-4).sum()}   "
              f"within 2e-3: {(np.abs(err) < 2e-3).sum()}")
        worst = sorted(sp, key=lambda r: -abs(r["lambdaMax"] - r["published"]))[:6]
        for r in worst:
            print(f"    {r['name']:24}{r['lambdaMax']:10.5f} vs {r['published']:9.5f}"
                  f"  ratio {r['publishedRatio']:6.3f}")
        print()

    dy = {r["name"]: r for r in scorable if r["source"] == "dysts"}
    if dy and os.path.exists(a.manifest):
        with open(a.manifest) as fh:
            man = {s["system"]: s for s in json.load(fh)["systems"]}
        pairs = [(n, dy[n]["lambdaMax"], man[n]["lambdaMax"])
                 for n in dy if n in man and man[n].get("lambdaMax")]
        if pairs:
            ours = np.array([p[1] for p in pairs])
            theirs = np.array([p[2] for p in pairs])
            good = ours > 0
            ratio = theirs[good] / ours[good]
            print(f"vs DYSTS published spectra ({len(pairs)} systems in common)")
            print(f"  published / recomputed: median {np.median(ratio):.3f}, "
                  f"IQR {np.percentile(ratio, 25):.3f}-{np.percentile(ratio, 75):.3f}")
            print(f"  within 10%: {(np.abs(ratio - 1) < 0.10).sum()} of {good.sum()}"
                  f"   within 2x: {((ratio > 0.5) & (ratio < 2)).sum()}")
            idx = np.argsort(-np.abs(np.log(np.maximum(ratio, 1e-9))))[:8]
            names = [p[0] for p, g in zip(pairs, good) if g]
            o2 = ours[good]
            t2 = theirs[good]
            print(f"    {'system':24}{'ours':>10}{'dysts':>10}{'ratio':>8}")
            for i in idx:
                print(f"    {names[i]:24}{o2[i]:10.5f}{t2[i]:10.5f}{ratio[i]:8.2f}")
            print()

    # The headline: the two catalogues hold 13 systems in common, and the
    # question that started this is whether one method makes them agree.
    SHARED = ["lorenz", "rossler", "chen", "hadley", "rabinovich_fabrikant",
              "chua", "moore_spiegel", "thomas", "halvorsen", "burke_shaw",
              "rucklidge", "nose_hoover", "henon_heiles"]
    bysp = {r["name"]: r for r in recs if r["source"] == "sprott"}
    rows = []
    for n in SHARED:
        d = dy.get(n.replace("_", "").lower()) or next(
            (v for k, v in dy.items() if k.lower() == n.replace("_", "")), None)
        if n in bysp and d:
            rows.append((n, bysp[n]["lambdaMax"], d["lambdaMax"]))
    if rows:
        print(f"SHARED SYSTEMS under one method ({len(rows)})")
        print(f"  {'system':24}{'sprott params':>14}{'dysts params':>14}{'ratio':>8}")
        for n, s_, d_ in rows:
            print(f"  {n:24}{s_:14.5f}{d_:14.5f}{(d_/s_ if s_ else float('nan')):8.2f}")

    if a.json:
        with open(a.json, "w") as fh:
            json.dump(dict(n=len(recs),
                           flagged={r["name"]: w for r, w in bad}), fh, indent=1)
        print(f"\nwrote {a.json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
