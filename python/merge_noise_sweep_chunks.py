"""Assemble the chunked RQA noise sweep into one summary table.

The RQA family was swept in chunks of six systems, one MATLAB call each, so
that a forty-hour run could be stopped and resumed at chunk granularity. Each
call wrote its own summary (tests/reports/noise_sweep_rqa/chunk_NN.csv) and
per-realization table (tests/artifacts/noise_sweep_rqa/, gitignored). This
concatenates the summaries into tests/reports/noise_sweep_rqa.csv, the
committed record alongside noise_sweep.csv from the Lyapunov pass, and checks
that every manifest system is present exactly once before writing.

The chunks are the same schema as quarctest.noise_sweep returns, so the merged
table is what a single call over all 126 systems would have produced.

Usage:
    python/merge_noise_sweep_chunks.py [--allow-partial]
"""

from __future__ import annotations

import argparse
import csv
import glob
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
CHUNKS = os.path.join(ROOT, "tests", "reports", "noise_sweep_rqa")
OUT = os.path.join(ROOT, "tests", "reports", "noise_sweep_rqa.csv")
MANIFEST = os.path.join(ROOT, "tests", "reports", "dysts", "manifest.json")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--allow-partial", action="store_true",
                    help="write even if some manifest systems are missing")
    a = ap.parse_args()

    files = sorted(glob.glob(os.path.join(CHUNKS, "chunk_*.csv")))
    if not files:
        print("no chunk files found", file=sys.stderr)
        return 1
    header, rows = None, []
    for f in files:
        with open(f, newline="") as fh:
            r = csv.reader(fh)
            h = next(r)
            if header is None:
                header = h
            elif h != header:
                print(f"schema differs in {os.path.basename(f)}", file=sys.stderr)
                return 1
            rows.extend(list(r))

    isys = header.index("system")
    seen = {}
    for row in rows:
        seen[row[isys]] = seen.get(row[isys], 0) + 1
    with open(MANIFEST) as fh:
        manifest = [s["system"] for s in json.load(fh)["systems"]]
    missing = [s for s in manifest if s not in seen]
    extra = [s for s in seen if s not in manifest]
    per_system = {n for n in seen.values()}
    print(f"{len(files)} chunks, {len(rows)} rows, {len(seen)} systems; "
          f"{len(missing)} manifest systems missing, {len(extra)} unexpected; "
          f"rows per system: {sorted(per_system)}")
    if extra:
        print(f"unexpected systems: {extra}", file=sys.stderr)
        return 1
    if missing and not a.allow_partial:
        print(f"missing: {', '.join(missing)}", file=sys.stderr)
        return 1
    if len(per_system) != 1:
        print("systems have differing row counts; a chunk is incomplete",
              file=sys.stderr)
        return 1

    order = {s: i for i, s in enumerate(manifest)}
    rows.sort(key=lambda r: (order[r[isys]], r[header.index("metric")],
                             r[header.index("arm")], float(r[header.index("noise")])))
    with open(OUT, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        w.writerows(rows)
    print(f"wrote {os.path.relpath(OUT, ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
