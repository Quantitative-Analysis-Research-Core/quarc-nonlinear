"""Rebuild a dysts manifest from the .bin files an interrupted export left behind.

export_dysts.py writes manifest.json only after every system finishes, so a run
that is stopped partway leaves a directory of .bin files that nothing can read:
the manifest is what carries R, N, the dtype, the sampling rate, the observable
index and the reference values. That is a real failure mode -- the first full
export ran 44 hours and was stopped with one pathological system still
integrating -- and this recovers from it without regenerating anything.

Nothing here is guesswork. Every manifest field is a pure function of the dysts
catalogue plus the export parameters, and the observable index is computed by
importing export_dysts.observable_index rather than reimplementing it, so the
two cannot drift. The only field that is not derivable from the catalogue is the
per-system count of degenerate realizations, which is read back from the export
log.

A system is included only if its .bin exists AND is exactly R*N*itemsize bytes.
A truncated file is excluded rather than trusted: the reader seeks to computed
offsets, so a short file yields plausible numbers from the wrong part of the
array instead of an error.

Usage:
    python/rebuild_manifest.py --dir tests/reports/dysts --R 100 --N 16000 \
        --log ~/dysts_export.log
"""

import argparse
import json
import os
import re
import sys

DTYPE_SIZE = {"float64": 8, "float32": 4}


def degenerate_counts(log_path):
    """Per-system degenerate counts, from the export log if it is available."""
    out = {}
    if not log_path or not os.path.isfile(os.path.expanduser(log_path)):
        return out
    pat = re.compile(r"^\[\d+/\d+\]\s+(\S+)\s+degenerate (\d+)")
    with open(os.path.expanduser(log_path)) as fh:
        for line in fh:
            m = pat.match(line)
            if m:
                out[m.group(1)] = int(m.group(2))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default="tests/reports/dysts")
    ap.add_argument("--R", type=int, default=100)
    ap.add_argument("--N", type=int, default=16000)
    ap.add_argument("--dtype", default="float64")
    ap.add_argument("--seed", type=int, default=20260810)
    ap.add_argument("--log", default="")
    args = ap.parse_args()

    here = os.path.dirname(os.path.abspath(__file__))
    sys.path.insert(0, here)
    import export_dysts as ex          # reuse the exporter's own logic
    import dysts

    cat_path = os.path.join(os.path.dirname(dysts.__file__),
                            "data", "chaotic_attractors.json")
    catalog = json.load(open(cat_path))

    itemsize = DTYPE_SIZE[args.dtype]
    expected = args.R * args.N * itemsize
    degen = degenerate_counts(args.log)

    entries, skipped = [], []
    for name in sorted(catalog):
        path = os.path.join(args.dir, f"{name}.bin")
        if not os.path.isfile(path):
            continue
        size = os.path.getsize(path)
        if size != expected:
            skipped.append((name, f"size {size} != expected {expected}"))
            continue

        meta = catalog[name]
        dim = len(meta["initial_conditions"])
        obs = ex.observable_index(meta, dim)
        if obs is None:
            skipped.append((name, "every component unbounded"))
            continue

        period = float(meta.get("period") or 1.0)
        spec = meta.get("lyapunov_spectrum_estimated") or []
        entries.append({
            "system": name,
            "file": f"{name}.bin",
            "dtype": args.dtype,
            "R": args.R, "N": args.N,
            "obsIndex": int(obs),
            "stateDim": int(dim),
            "period": period,
            "fs": ex.TARGET_PTS_PER_PERIOD / period,
            "dt": float(meta.get("dt") or float("nan")),
            # See export_dysts.py: the label is the continuous-QR spectrum
            # maximum, not the trajectory-divergence estimate.
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
            "degenerate": int(degen.get(name, 0)),
            "citation": (meta.get("citation") or "")[:300],
            "doi": meta.get("doi") or "",
        })

    man = os.path.join(args.dir, "manifest.json")
    json.dump({"spread": ex.SPREAD, "ptsPerPeriod": ex.TARGET_PTS_PER_PERIOD,
               "seed": args.seed, "rebuilt": True,
               "systems": entries}, open(man, "w"), indent=1)

    print(f"wrote {man}")
    print(f"  systems: {len(entries)}")
    print(f"  with degenerate realizations: "
          f"{sum(1 for e in entries if e['degenerate'])}")
    if skipped:
        print(f"  excluded: {len(skipped)}")
        for n, why in skipped:
            print(f"    {n}: {why}")


if __name__ == "__main__":
    sys.exit(main())
