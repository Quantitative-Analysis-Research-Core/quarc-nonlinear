"""The Python bindings against MATLAB's answers on the shared fixtures."""

import json
import warnings
from contextlib import contextmanager
from pathlib import Path

import numpy as np
import pytest

import quarcnonlinear as qn

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


@needs_fixtures
def test_ent_samp_matches_matlab():
    r = ref()
    assert qn.ent_samp(fixture("ar1_phi70_512"), 2, 0.2) == pytest.approx(
        r["ent_samp"]["ar1_phi70_512"], rel=TOL)
    for name, want in r["ent_samp_dims"].items():
        x = fixture(name)
        got = [qn.ent_samp(x, 1, 0.2), qn.ent_samp(x, 2, 0.2), qn.ent_samp(x, 3, 0.15)]
        assert got == pytest.approx(want, rel=TOL), name


def test_ent_samp_no_match_is_inf():
    assert qn.ent_samp([1, 100, 1e4, 1e6, 1e8], 2, 1e-9, flag="abs") == np.inf


def test_ent_samp_too_short_raises():
    with pytest.raises(ValueError):
        qn.ent_samp([1.0, 2.0], 2)


@needs_fixtures
def test_fnn_matches_matlab():
    for name, r in ref()["fnn"].items():
        x = fixture(name)
        with _quiet():
            d0, e0 = qn.fnn(x, r["delay"], r["maxdim"], 15, 2, 0)
            d1, e1 = qn.fnn(x, r["delay"], r["maxdim"], 10, 2, 1)
        assert e0 == pytest.approx(r["dE_speed0"], rel=TOL, abs=1e-15), name
        assert d0 == r["dim_speed0"], name
        assert e1 == pytest.approx(r["dE_speed1"], rel=TOL, abs=1e-15), name
        assert d1 == r["dim_speed1"], name


@needs_fixtures
def test_corr_dim_matches_matlab():
    for name, r in ref()["corr_dim"].items():
        cod, extra = qn.corr_dim(fixture(name), r["delay"], r["dim"])
        want = np.array([-np.inf if v is None else v for v in r["logC"]])  # jsonencode: -Inf -> null
        assert extra["logEps"] == pytest.approx(r["logEps"], rel=TOL), name
        assert extra["logC"] == pytest.approx(want, rel=TOL), name
        assert list(extra["idx"] + 1) == r["idx"], name
        assert cod == pytest.approx(r["CoD"], rel=TOL), name
        assert extra["r2"] == pytest.approx(r["r2"], rel=TOL), name


@needs_fixtures
def test_lye_r_matches_matlab():
    r = ref()["lye_r"]
    lz = fixture("lorenz_2048")
    cases = {"lorenz": qn.lye_r(lz, 100, 10, 3),
             "lorenz_theiler25": qn.lye_r(lz, 100, 10, 3, theiler=25),
             "henon": qn.lye_r(fixture("henon_512"), 1, 1, 2)}
    for name, out in cases.items():
        assert list(out[:, 1] + 1) == r[name]["pairs"], name
        assert out[:, 2] == pytest.approx(r[name]["divergence"], rel=TOL, abs=1e-15), name
    es, el, _ = qn.lye_r(lz, 100, 10, 3, slope=[0, 0.4, 1, 3], mean_period=0.75)
    assert es == pytest.approx(r["lorenz_slopes"]["LyES"], rel=TOL)
    assert el == pytest.approx(r["lorenz_slopes"]["LyEL"], rel=TOL)


def test_line_hist_kernel_matches_numpy_reference():
    from quarcnonlinear.line_hist import _lines
    rng = np.random.default_rng(2)
    for n in (1, 2, 7, 60):
        A = (rng.random((n, n)) < 0.4).astype(float)
        A = A * A.T
        for theiler in (0, 1, 3):
            for cross in (False, True):
                th = 0 if cross else theiler
                got = qn.line_hist(None, A - 1, 0.5, "crqa" if cross else "rqa", th)
                diag, vert = _lines(A, th, cross)
                assert np.array_equal(got["diag_hist"], diag), (n, theiler, cross)
                assert np.array_equal(got["vertical_hist"], vert), (n, theiler, cross)


def test_line_hist_recurrence_by_hand():
    A = np.array([[1, 1, 0], [1, 1, 1], [0, 1, 1]], float)
    got = qn.line_hist(None, A - 1, 0, "rqa")
    assert got["recurrence"] == pytest.approx(100 * 4 / 6)   # off-diagonal 1s / off-diagonal cells
    assert list(got["diag_hist"]) == [2, 2]
    assert list(got["vertical_hist"]) == [2, 3, 2]


def test_psr_layout():
    x = np.arange(1.0, 8.0)
    assert np.array_equal(qn.psr(x, 2, 3), np.column_stack([x[0:3], x[2:5], x[4:7]]))


@pytest.mark.skipif(not (REPO / "cpp").exists(), reason="cpp/ not found")
def test_core_copy_is_current():
    core, copy = REPO / "cpp", REPO / "pypkg" / "src" / "core"
    for f in (core / "src").glob("*.cpp"):
        assert (copy / f"core_{f.name}").read_text() == f.read_text(), f.name
    for f in (core / "include" / "quarc").glob("*.hpp"):
        assert (copy / "quarc" / f.name).read_text() == f.read_text(), f.name
