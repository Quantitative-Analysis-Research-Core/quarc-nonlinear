"""The entropy family against MATLAB's answers on the shared fixtures."""

import json
from pathlib import Path

import numpy as np
import pytest

from quarcnonlinear.ent import ent
from quarcnonlinear.ent_ap import ent_ap
from quarcnonlinear.ent_ms_plus import ent_ms_plus
from quarcnonlinear.ent_permu import ent_permu
from quarcnonlinear.ent_samp import ent_samp
from quarcnonlinear.ent_symbolic import ent_symbolic
from quarcnonlinear.ent_weighted import ent_weighted
from quarcnonlinear.ent_xap import ent_xap
from quarcnonlinear.ent_xsamp import ent_xsamp
from quarcnonlinear.psr import psr

REPO = Path(__file__).resolve().parents[2]
FIX = REPO / "tests" / "fixtures"
TOL = 1e-10  # as the R package's tests: summation order passes, a defect does not

needs_fixtures = pytest.mark.skipif(not (FIX / "matlab_reference.json").exists(),
                                    reason="shared fixtures not found (not in the repository)")


def ref():
    return json.loads((FIX / "matlab_reference.json").read_text())


def fixture(name):
    return np.loadtxt(FIX / f"{name}.csv", delimiter=",")


def _nan(v):
    return np.array([np.nan if e is None else e for e in v], dtype=float)  # jsonencode: NaN -> null


@needs_fixtures
def test_ent_ap_matches_matlab():
    for name, want in ref()["ent_ap"].items():
        x = fixture(name)
        got = [ent_ap(x, 1, 0.2), ent_ap(x, 2, 0.2), ent_ap(x, 3, 0.15)]
        assert got == pytest.approx(want, rel=TOL), name


@needs_fixtures
def test_ent_permu_matches_matlab_value_and_histogram():
    r = ref()["ent_permu"]
    x = fixture("lorenz_2048")
    a, ha = ent_permu(x, 3, 1)
    b, hb = ent_permu(x, 4, 2)
    assert a == pytest.approx(r["lorenz_2048"]["dim3_delay1"], rel=TOL)
    assert ha.tolist() == r["lorenz_2048"]["hist3"]
    assert b == pytest.approx(r["lorenz_2048"]["dim4_delay2"], rel=TOL)
    assert hb.tolist() == r["lorenz_2048"]["hist4"]
    w, hw = ent_permu(fixture("white_4096"), 4, 1)
    assert w == pytest.approx(r["white_4096"]["dim4_delay1"], rel=TOL)
    assert hw.tolist() == r["white_4096"]["hist4"]


@needs_fixtures
def test_cross_symbolic_and_weighted_entropies_match_matlab():
    r = ref()
    a = fixture("ar1_phi70_512")
    h = fixture("henon_512")
    got = [ent_xsamp(a, h, 2, 0.2, 0), ent_xsamp(a, h, 2, 0.2, 1),
           ent_xsamp(a, h, 2, 0.2, 2), ent_xsamp(a, h, 1, 0.3, 2)]
    assert got == pytest.approx(r["ent_xsamp"]["ar1_henon"], rel=TOL)
    got = [ent_xap(a, h, 2, 0.2, 1), ent_xap(a, h, 1, 0.3, 1), ent_xap(a, h, 2, 0.25, 2)]
    assert got == pytest.approx(r["ent_xap"]["ar1_henon"], rel=TOL)
    b = (a > np.median(a)).astype(float)
    got = [ent_symbolic(b, 2), ent_symbolic(b, 3), ent_symbolic(b, 4)]
    assert got == pytest.approx(r["ent_symbolic"]["ar1_median"], rel=TOL)
    y = fixture("lorenz_2048")
    Y = psr(y[:400], 10, 3)
    D = np.sqrt(((Y[:, None, :] - Y[None, :, :]) ** 2).sum(axis=2))
    assert ent_weighted(np.exp(-D / np.std(y, ddof=1))) == pytest.approx(
        r["ent_weighted"]["lorenz_400"], rel=TOL)


@needs_fixtures
def test_ent_ms_plus_matches_matlab_on_all_five_outputs():
    for name, r in ref()["ent_ms_plus"].items():
        out = ent_ms_plus(fixture(name), r["scale"], r["dim"], r["radius"])
        for f, got in zip(["rcmse", "cmse", "mse", "msfe", "gmse"], out):
            assert got == pytest.approx(_nan(r[f]), rel=TOL, nan_ok=True), f"{name} {f}"


def test_ent_ms_plus_gmse_scale1_is_nan():
    x = np.sin(np.arange(1, 301) / 4) + 0.3 * np.cos(np.arange(1, 301) / 1.3)
    assert np.isnan(ent_ms_plus(x, 3, 2, 0.2)[4][0])


def test_ent_dispatches_to_each_estimator_and_guards_its_inputs():
    t = np.arange(1, 301)
    x = np.sin(t / 4) + 0.3 * np.cos(t / 1.3)
    y = np.cos(t / 4)
    assert ent(x)[0] == ent_samp(x, 2, 0.2)
    assert ent(x, "approximate")[0] == ent_ap(x, 2, 0.2)
    v, extra = ent(x, "permutation", dim=3)
    assert v == ent_permu(x, 3, 1)[0]
    assert np.array_equal(extra["patternCounts"], ent_permu(x, 3, 1)[1])
    assert ent(x, y=y)[0] == ent_xsamp(x, y, 2, 0.2, 1)
    assert ent(x, "approximate", y=y)[0] == ent_xap(x, y, 2, 0.2, 1)
    v, extra = ent(x, "multiscale", scale=3)
    assert np.array_equal(v, ent_ms_plus(x, 3, 2, 0.2)[0])
    assert extra["estimator"] == "refined composite multiscale sample"
    b = (x > np.median(x)).astype(float)
    assert ent(b, "symbolic", dim=3)[0] == ent_symbolic(b, 3)
    with pytest.raises(ValueError, match="no cross-entropy form"):
        ent(x, "permutation", y=y)
    with pytest.raises(ValueError, match="binary"):
        ent(x, "symbolic")
    with pytest.raises(ValueError, match="equal lengths"):
        ent(x, y=y[1:])
    with pytest.raises(ValueError, match="contains NaN"):
        ent([1, np.nan, 3, 4, 5])
    with pytest.raises(ValueError, match="Unknown algorithm"):
        ent(x, "fuzzy")


def test_ent_symbolic_rejects_non_binary():
    with pytest.raises(ValueError, match="binary"):
        ent_symbolic([0, 1, 2, 1, 0], 2)


def test_ent_xsamp_rejects_unequal_lengths():
    with pytest.raises(ValueError, match="same length"):
        ent_xsamp([1, 2, 3, 4], [1, 2, 3], 2, 0.2, 1)


def test_ent_weighted_is_zero_when_column_sums_are_equal():
    assert ent_weighted(np.ones((5, 5))) == 0
