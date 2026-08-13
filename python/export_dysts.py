"""Export ensembles of scalar observables from the dysts catalogue.

The characterization suite validates estimators against systems whose exponents
were computed from the linearized equations rather than from a time series.
Sprott's Appendix A gives 62 such systems; the dysts database of Gilpin (2021)
gives 135 more, each carrying a Lyapunov spectrum computed by integrating the
tangent space with Gram-Schmidt/QR orthonormalization (Christiansen & Rugh 1997)
and a correlation dimension. Those are the same two reference quantities the
Sprott catalogue supplies, obtained the same way, so the two sets can be
analysed under one protocol.

WHY THIS IS PYTHON. The dysts vector fields live in Python. Porting 129 of them
to MATLAB by hand would introduce transcription errors precisely where they are
hardest to notice -- a mistyped coefficient still produces a plausible chaotic
series. Instead the trajectories are generated here, by the same code that
produced the published reference exponents, and written to disk for the MATLAB
battery to read. The reference values and the data then share a provenance.

PROTOCOL, MATCHING quarctest.ensemble_ics AND quarctest.sprott_series.
  - realization 1 is the published initial condition, verbatim
  - realizations 2..R are drawn on the attractor: integrate from the published
    IC to land on it, measure the attractor's per-component standard deviation,
    perturb by Spread times that, and relax briefly so the point sits on the
    attractor rather than near it
  - sampling is fixed per system at 40 points per dominant period, the same
    target the Sprott protocol decimates towards, so fs is comparable
  - the observable is the first component not listed in unbounded_indices,
    because delay embedding requires a bounded recurrent observable

Excluded: the six delay-differential systems, which need a history function
rather than an initial condition and do not fit this protocol.

Usage:
    python/export_dysts.py --out tests/reports/dysts --R 100 --N 16000

N is the ceiling for every downstream experiment: truncating a stored series is
free, lengthening one means regenerating the whole export. It is set high for
that reason, so series-length experiments cost nothing later.
"""

import argparse
import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np

TARGET_PTS_PER_PERIOD = 40
SPREAD = 0.01           # IC perturbation, as a fraction of attractor extent
ANCHOR_PERIODS = 50     # integration to land on the attractor
RELAX_PERIODS = 5       # relaxation after perturbing


def system_names(catalog):
    """Systems eligible for this protocol, in catalogue order."""
    out = []
    for name, meta in catalog.items():
        if meta.get("delay"):
            continue                       # delay-DE: needs a history function
        if meta.get("maximum_lyapunov_estimated") in (None, ""):
            continue                       # no reference exponent
        out.append(name)
    return out


def observable_index(meta, dim):
    """First component that is not flagged unbounded."""
    unb = set(meta.get("unbounded_indices") or [])
    for i in range(dim):
        if i not in unb:
            return i
    return None


def export_one(name, meta, R, N, seed, outdir):
    """Write one system's ensemble. Returns a manifest entry or None."""
    import dysts.flows as flows

    try:
        model = getattr(flows, name)()
    except AttributeError:
        return {"system": name, "skipped": "not in dysts.flows"}

    ic0 = np.asarray(model.ic, dtype=float)
    dim = ic0.shape[-1]
    obs = observable_index(meta, dim)
    if obs is None:
        return {"system": name, "skipped": "every component unbounded"}

    rng = np.random.default_rng(seed + abs(hash(name)) % (2**31))

    def traj(ic, npts):
        model.ic = np.asarray(ic, dtype=float)
        out = model.make_trajectory(npts, pts_per_period=TARGET_PTS_PER_PERIOD,
                                    resample=True)
        return None if out is None else np.asarray(out)

    # Land on the attractor and measure its extent.
    probe = traj(ic0, ANCHOR_PERIODS * TARGET_PTS_PER_PERIOD)
    if probe is None or not np.all(np.isfinite(probe)):
        return {"system": name, "skipped": "anchor integration failed"}
    anchor = probe[-1]
    scale = np.std(probe, axis=0)
    flat = ~np.isfinite(scale) | (scale <= 0)
    if np.any(flat):
        good = scale[~flat]
        scale[flat] = good.max() if good.size else 1.0

    relax = RELAX_PERIODS * TARGET_PTS_PER_PERIOD
    # float64 throughout. The metrics themselves are insensitive to the
    # difference, but corr_dim works on inter-point distances at the smallest
    # radii, where float32's ~7 digits is the only thing that would stand
    # between an estimate and a question about it. Disk is cheaper than that
    # question.
    series = np.full((R, N), np.nan, dtype=np.float64)
    degenerate = 0

    for r in range(R):
        if r == 0:
            y = traj(ic0, N)                      # published IC, verbatim
        else:
            ic = anchor + SPREAD * scale * rng.standard_normal(dim)
            y = traj(ic, relax + N)
            if y is not None:
                y = y[relax:]
        if y is None or y.shape[0] < N:
            degenerate += 1
            continue
        v = y[:N, obs]
        if not np.all(np.isfinite(v)) or np.std(v) <= 0:
            degenerate += 1
            continue
        series[r, :] = v.astype(np.float64)

    if degenerate >= R:
        return {"system": name, "skipped": "every realization degenerate"}

    path = os.path.join(outdir, f"{name}.bin")
    series.tofile(path)

    period = float(meta.get("period") or 1.0)
    spec = meta.get("lyapunov_spectrum_estimated") or []
    return {
        "system": name,
        "file": os.path.basename(path),
        "dtype": "float64",
        "R": R, "N": N,
        "obsIndex": int(obs),
        "stateDim": int(dim),
        "period": period,
        "fs": TARGET_PTS_PER_PERIOD / period,   # samples per unit time
        "dt": float(meta.get("dt") or np.nan),
        "lambdaMax": float(meta["maximum_lyapunov_estimated"]),
        "lyapunovSpectrum": [float(v) for v in spec],
        "corrDim": (float(meta["correlation_dimension"])
                    if meta.get("correlation_dimension") is not None else None),
        "kaplanYorke": (float(meta["kaplan_yorke_dimension"])
                        if meta.get("kaplan_yorke_dimension") is not None else None),
        "nonautonomous": bool(meta.get("nonautonomous")),
        "hamiltonian": bool(meta.get("hamiltonian")),
        "degenerate": int(degenerate),
        "citation": (meta.get("citation") or "")[:300],
        "doi": meta.get("doi") or "",
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="tests/reports/dysts")
    ap.add_argument("--R", type=int, default=100)
    ap.add_argument("--N", type=int, default=16000)
    ap.add_argument("--seed", type=int, default=20260810)
    ap.add_argument("--limit", type=int, default=0, help="first K systems only")
    ap.add_argument("--only", default="", help="comma-separated system names")
    ap.add_argument("--workers", type=int, default=os.cpu_count() - 2)
    args = ap.parse_args()

    import dysts

    cat_path = os.path.join(os.path.dirname(dysts.__file__),
                            "data", "chaotic_attractors.json")
    catalog = json.load(open(cat_path))

    names = system_names(catalog)
    if args.only:
        want = {s.strip() for s in args.only.split(",")}
        names = [n for n in names if n in want]
    if args.limit:
        names = names[: args.limit]

    os.makedirs(args.out, exist_ok=True)
    print(f"exporting {len(names)} systems, R={args.R}, N={args.N}, "
          f"{args.workers} workers", flush=True)

    entries, t0 = [], time.time()
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        futs = {ex.submit(export_one, n, catalog[n], args.R, args.N,
                          args.seed, args.out): n for n in names}
        for k, fut in enumerate(as_completed(futs), 1):
            e = fut.result()
            entries.append(e)
            tag = e.get("skipped", f"degenerate {e.get('degenerate', 0)}")
            print(f"[{k}/{len(names)}] {e['system']:<28} {tag}  "
                  f"{time.time()-t0:.0f}s", flush=True)

    entries = [e for e in entries if "skipped" not in e]
    entries.sort(key=lambda e: e["system"])
    man = os.path.join(args.out, "manifest.json")
    json.dump({"spread": SPREAD, "ptsPerPeriod": TARGET_PTS_PER_PERIOD,
               "seed": args.seed, "systems": entries}, open(man, "w"), indent=1)
    print(f"\nwrote {len(entries)} systems to {man} in "
          f"{(time.time()-t0)/60:.1f} min")


if __name__ == "__main__":
    sys.exit(main())
