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
import signal
import sys
import time
import zlib
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np

TARGET_PTS_PER_PERIOD = 40
SPREAD = 0.01           # IC perturbation, as a fraction of attractor extent
ANCHOR_PERIODS = 50     # integration to land on the attractor
RELAX_PERIODS = 5       # relaxation after perturbing
TRAJ_TIMEOUT_S = 900    # wall-clock cap on one solve; ~4x a normal realization


class _IntegrationTimeout(Exception):
    pass


def _raise_timeout(signum, frame):
    raise _IntegrationTimeout()


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


def observable_indices(meta, dim):
    """Every component that is not flagged unbounded.

    The first export kept one observable per system, which cost more than it
    saved: a learned estimator's channel augmentation became a no-op, and a
    system's 100 initial conditions are statistically alike in a way that two
    different state variables are not. Delay embedding still needs each channel
    to be bounded and recurrent, so the unbounded ones stay excluded.
    """
    unb = set(meta.get("unbounded_indices") or [])
    return [i for i in range(dim) if i not in unb]


def export_one(name, meta, R, N, seed, outdir, solver="", rtol=1e-9, atol=1e-9):
    """Write one system's ensemble. Returns a manifest entry or None."""
    import dysts.flows as flows

    try:
        model = getattr(flows, name)()
    except AttributeError:
        return {"system": name, "skipped": "not in dysts.flows"}

    ic0 = np.asarray(model.ic, dtype=float)
    dim = ic0.shape[-1]
    obs_list = observable_indices(meta, dim)
    if not obs_list:
        return {"system": name, "skipped": "every component unbounded"}

    # crc32, not hash(): Python salts str hashes per process, so hash(name)
    # gave every run -- and every pool worker -- a different draw sequence
    # under the same --seed.
    rng = np.random.default_rng(seed + zlib.crc32(name.encode()) % (2**31))

    def traj(ic, npts):
        # A perturbed IC can wedge the adaptive integrator indefinitely: one
        # Sakarya realization burned 110 CPU-hours without returning while a
        # normal one takes ~200 s. A solve that has not come back within the
        # alarm is returned as None, which the callers already record as a
        # degenerate realization.
        model.ic = np.asarray(ic, dtype=float)
        kw = dict(pts_per_period=TARGET_PTS_PER_PERIOD, resample=True)
        if solver:
            kw.update(method=solver, rtol=rtol, atol=atol)
        signal.signal(signal.SIGALRM, _raise_timeout)
        signal.alarm(TRAJ_TIMEOUT_S)
        try:
            out = model.make_trajectory(npts, **kw)
            if out is None and solver:
                # dysts defaults to Radau because some of these systems are
                # stiff. An explicit solver is ~60x faster where it works; where
                # it does not, fall back rather than lose the system. The
                # fallback shares the alarm, so one call is still capped.
                out = model.make_trajectory(npts, pts_per_period=TARGET_PTS_PER_PERIOD,
                                            resample=True)
        except _IntegrationTimeout:
            return None
        finally:
            signal.alarm(0)
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
    C = len(obs_list)
    series = np.full((R, C, N), np.nan, dtype=np.float64)
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
        v = y[:N, obs_list]                      # (N, C)
        if not np.all(np.isfinite(v)) or np.min(np.std(v, axis=0)) <= 0:
            degenerate += 1
            continue
        series[r] = v.T.astype(np.float64)       # (C, N)

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
        "obsIndices": [int(i) for i in obs_list],
        "nChannels": int(C),
        "stateDim": int(dim),
        "period": period,
        "fs": TARGET_PTS_PER_PERIOD / period,   # samples per unit time
        "dt": float(meta.get("dt") or np.nan),
        # LABEL PROVENANCE. lambdaMax is max(lyapunov_spectrum_estimated), the
        # Christiansen-Rugh continuous-QR exponent from the analytic Jacobian.
        # It is NOT maximum_lyapunov_estimated, which dysts computes by
        # trajectory divergence -- a different estimator with a different bias.
        # Measured across this corpus the two agree within 10% on only 30 of
        # 123 systems and differ by more than 50% on 51 of them, ranging from
        # 0.013x to 23.7x. Training on the divergence field would teach a model
        # to reproduce an estimator rather than the quantity, and it is carried
        # here under its own name so the two can never be confused.
        "lambdaMax": (max(float(v) for v in spec) if spec
                      else float(meta["maximum_lyapunov_estimated"])),
        "lambdaMaxSource": ("qr_spectrum" if spec else "trajectory_divergence"),
        "lambdaDivergence": float(meta["maximum_lyapunov_estimated"]),
        "lyapunovSpectrum": [float(v) for v in spec],
        "corrDim": (float(meta["correlation_dimension"])
                    if meta.get("correlation_dimension") is not None else None),
        "pesinEntropy": (float(meta["pesin_entropy"])
                         if meta.get("pesin_entropy") is not None else None),
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
    ap.add_argument("--solver", default="",
                    help="scipy method, e.g. RK45. Empty uses dysts' Radau default.")
    ap.add_argument("--rtol", type=float, default=1e-9)
    ap.add_argument("--atol", type=float, default=1e-9)
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

    # Append each entry as it lands.  A previous run wrote the manifest only on
    # completion, was interrupted at hour 44, and left 125 binaries that no
    # longer described themselves.  The .bin files are the expensive part; the
    # few hundred bytes that make them readable should never depend on the run
    # finishing.
    jsonl = os.path.join(args.out, "manifest.jsonl")
    entries, t0 = [], time.time()
    with open(jsonl, "a") as log, ProcessPoolExecutor(max_workers=args.workers) as ex:
        futs = {ex.submit(export_one, n, catalog[n], args.R, args.N,
                          args.seed, args.out, args.solver, args.rtol,
                          args.atol): n for n in names}
        for k, fut in enumerate(as_completed(futs), 1):
            e = fut.result()
            entries.append(e)
            log.write(json.dumps(e) + "\n")
            log.flush()
            os.fsync(log.fileno())
            tag = e.get("skipped", f"degenerate {e.get('degenerate', 0)}")
            print(f"[{k}/{len(names)}] {e['system']:<28} {tag}  "
                  f"{time.time()-t0:.0f}s", flush=True)

    # Build the manifest from the whole jsonl, not from this run's entries, so a
    # follow-up run with --only adds to an export instead of replacing its
    # manifest with the handful of systems it just redid. The latest line for a
    # system wins.
    latest = {}
    with open(jsonl) as log:
        for line in log:
            if line.strip():
                e = json.loads(line)
                latest[e["system"]] = e
    entries = [e for e in latest.values() if "skipped" not in e]
    entries.sort(key=lambda e: e["system"])
    man = os.path.join(args.out, "manifest.json")
    json.dump({"spread": SPREAD, "ptsPerPeriod": TARGET_PTS_PER_PERIOD,
               "seed": args.seed, "systems": entries}, open(man, "w"), indent=1)
    print(f"\nwrote {len(entries)} systems to {man} in "
          f"{(time.time()-t0)/60:.1f} min")


if __name__ == "__main__":
    sys.exit(main())
