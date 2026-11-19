"""The RQA family against MATLAB's answers on the shared fixtures."""

import json
import warnings
from contextlib import contextmanager
from pathlib import Path

import numpy as np
import pytest

from quarcnonlinear.crqa import crqa
from quarcnonlinear.jrqa import jrqa
from quarcnonlinear.mdrqa import mdrqa
from quarcnonlinear.rqa import rqa
from quarcnonlinear.set_radius import set_radius

REPO = Path(__file__).resolve().parents[2]
FIX = REPO / "tests" / "fixtures"
TOL = 1e-10  # as the R package's tests: summation order passes, a defect does not

needs_fixtures = pytest.mark.skipif(not (FIX / "matlab_reference.json").exists(),
                                    reason="shared fixtures not found (not in the repository)")


def ref():
    return json.loads((FIX / "matlab_reference.json").read_text())


@contextmanager
def _quiet():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        yield


def fixture(name):
    return np.loadtxt(FIX / f"{name}.csv", delimiter=",")


def _cases():
    lz, hn, ar = fixture("lorenz_2048"), fixture("henon_512"), fixture("ar1_phi70_512")
    with _quiet():
        return {
            "rqa": {
                "lorenz_rec": rqa(lz[:600], 10, 3, "rec", 2.5),
                "lorenz_rad_euc_theiler": rqa(lz[:600], 10, 3, "rad", 0.3, norm="euc", theiler=3),
                "henon_rec_max": rqa(hn, 1, 2, "rec", 5, norm="max"),
                "ar1_rad_min": rqa(ar, 1, 1, "rad", 0.5, norm="min", dmin=3, vmin=3, zscore=0)},
            "crqa": {
                "ar1_henon_rec": crqa(np.column_stack([ar[:400], hn[:400]]), 1, 2, "rec", 3),
                "lorenz_rad_euc": crqa(np.column_stack([lz[:400], lz[5:405]]), 4, 3, "rad", 0.4,
                                       norm="euc")},
            "jrqa": {
                "lorenz_rec": jrqa(np.column_stack([lz[:300], lz[10:310]]), [5, 5], [2, 2],
                                   "rec", 5),
                "ar1_henon_rad": jrqa(np.column_stack([ar[:300], hn[:300]]), [1, 1], [1, 1],
                                      "rad", 0.5, norm="max")},
            "mdrqa": {
                "lorenz_henon_rec": mdrqa(np.column_stack([lz[:400], hn[:400]]), 1, 1, "rec", 4),
                "lorenz_henon_rad": mdrqa(np.column_stack([lz[:400], hn[:400]]), 2, 2, "rad", 0.5,
                                          norm="max")}}


FIELDS = ["REC", "RADIUS", "Size", "DET", "MeanL", "MaxL", "EntrL", "LAM", "MeanV",
          "MaxV", "EntrV", "EntrW"]


@needs_fixtures
def test_rqa_family_matches_matlab():
    r = ref()
    for fam, cases in _cases().items():
        for s, out in cases.items():
            want = r[fam][s]
            lab = f"{fam} {s}"
            np.testing.assert_array_equal(out["RP"].sum(axis=1), want["rowsum"], err_msg=lab)
            np.testing.assert_array_equal(out["RP"].sum(axis=0), want["colsum"], err_msg=lab)
            for f in FIELDS:
                got = out["RESULTS"].get(f, np.nan)            # jrqa has no EntrW
                w = want.get(f)
                w = np.nan if w is None else w
                assert got == pytest.approx(w, rel=TOL, nan_ok=True), f"{lab} {f}"


def test_rqa_hits_recurrence_target_and_guards_inputs():
    t = np.arange(1, 401)
    x = np.sin(t / 6) + 0.2 * np.cos(t / 1.9)
    with _quiet():
        out = rqa(x, 6, 3, "rec", 2.5)
    assert out["RESULTS"]["REC"] == pytest.approx(2.5, rel=0.02)
    assert out["RP"].shape == (400 - 12, 400 - 12)
    with pytest.raises(ValueError, match="single column"):
        rqa(np.column_stack([x, x]))
    with pytest.raises(ValueError, match="two column"):
        crqa(x)
    with pytest.raises(ValueError, match="at least two"):
        jrqa(x)
    with pytest.raises(ValueError, match="differs between"):
        jrqa(np.column_stack([x, x]), [1, 1], [1, 2], "rad", 0.5, norm="euc")


def test_set_radius_reports_and_brackets():
    t = np.arange(1, 201)
    x = np.sin(t / 5)
    a = -np.abs(np.subtract.outer(x, x))
    with pytest.warns(UserWarning, match="% recurrence"):
        sr = set_radius(x, a, 0.01, 0.5, 5, "rqa", 20)
    assert sr["rec"] == pytest.approx(5, rel=0.05)
    assert set(sr) == {"rec", "diag_hist", "vertical_hist", "rad_final", "A"}


def test_rqa_bad_param_raises():
    with pytest.raises(ValueError, match="param"):
        rqa(np.sin(np.arange(50)), param="radius", threshold=0.5)
