"""Re-divide the committed Lyapunov estimates by the recomputed references.

The characterization summaries record, per system and estimator, the ensemble
median of the estimate and the reference it was divided by. The estimates do
not depend on the reference, so replacing the reference is a post-processing
step: every ratio in the lab report can be recomputed from the committed CSVs
without re-running the ten-hour battery. This script does that, writes the
result beside the originals, and prints what moved.

That answers the question the corpus was built for. If the estimators' apparent
bias against the dysts catalogue shrinks once the references are corrected,
that share of the reported error was never the estimators'.

Only lyap_wolf and lyap_ros are rescored: they are the metrics whose reference
is lambda1. Correlation dimension keeps its published D2, which the corpus does
not recompute.

Usage:
    python/rescore_ratios.py
"""

from __future__ import annotations

import csv
import json
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
REPORTS = os.path.join(ROOT, "tests", "reports")

LYAP_METRICS = ("lyap_wolf", "lyap_ros")


def _stats(r):
    r = np.asarray([x for x in r if np.isfinite(x)])
    if r.size == 0:
        return "n=0"
    q1, med, q3 = np.percentile(r, [25, 50, 75])
    within = np.mean((r > 0.5) & (r < 2.0)) * 100
    return f"n={r.size:3d}  median {med:.3f}  IQR {q1:.3f}-{q3:.3f}  within 2x {within:5.1f}%"


def main():
    with open(os.path.join(REPORTS, "reference_spectra.json")) as fh:
        corpus = json.load(fh)
    # Keyed on (catalogue, lower-case name): the Sprott CSV and the corpus
    # both use snake_case, the dysts CSV and the corpus both use CamelCase,
    # and lower-casing makes the join indifferent to which.
    ref = {(s["source"], s["name"].lower()): s for s in corpus["systems"]}

    out_rows = []
    for source, fname in (("sprott", "characterization_full.csv"),
                          ("dysts", "characterization_dysts.csv")):
        path = os.path.join(REPORTS, fname)
        with open(path) as fh:
            rows = list(csv.DictReader(fh))
        for r in rows:
            if r["metric"] not in LYAP_METRICS:
                continue
            med = float(r["median"]) if r["median"] not in ("", "NaN") else np.nan
            pub = float(r["reference"]) if r["reference"] not in ("", "NaN") else np.nan
            s = ref.get((source, r["system"].lower()))
            new = s["lambdaMax"] if s else np.nan
            out_rows.append(dict(
                catalogue=source, system=r["system"], metric=r["metric"],
                median=med, referencePublished=pub, referenceRecomputed=new,
                ratioPublished=(med / pub if np.isfinite(pub) and pub else np.nan),
                ratioRecomputed=(med / new if np.isfinite(new) and new else np.nan),
                referenceShift=(new / pub if np.isfinite(pub) and pub else np.nan),
                icDependent=bool(s and s.get("icDependent")),
                # False where the corpus found the recorded initial condition
                # not chaotic: the reference is essentially zero there and the
                # ratio is noise.
                chaoticAtIC=bool(s and s.get("chaoticAtIC", True)),
                inCorpus=bool(s),
            ))

    out = os.path.join(REPORTS, "reference_rescoring.csv")
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out_rows[0].keys()))
        w.writeheader()
        w.writerows(out_rows)

    print(f"wrote {os.path.relpath(out, ROOT)}: {len(out_rows)} rows\n")
    for source in ("sprott", "dysts"):
        for metric in LYAP_METRICS:
            sel = [r for r in out_rows
                   if r["catalogue"] == source and r["metric"] == metric
                   and r["inCorpus"] and not r["icDependent"] and r["chaoticAtIC"]]
            print(f"{source:7}{metric:10}  (dissipative, chaotic-at-IC systems in corpus)")
            print(f"   published  refs: {_stats([r['ratioPublished'] for r in sel])}")
            print(f"   recomputed refs: {_stats([r['ratioRecomputed'] for r in sel])}")
        shift = [r["referenceShift"] for r in out_rows
                 if r["catalogue"] == source and r["metric"] == "lyap_wolf"
                 and r["inCorpus"] and not r["icDependent"] and r["chaoticAtIC"]]
        print(f"   reference itself (recomputed/published): {_stats(shift)}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
