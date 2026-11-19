"""The AMI family against MATLAB's answers on the shared fixtures, and its contracts."""

import json
import warnings
from contextlib import contextmanager
from pathlib import Path

import numpy as np
import pytest

from quarcnonlinear.ami import ami
from quarcnonlinear.ami_copula import ami_copula
from quarcnonlinear.ami_histogram import ami_histogram
from quarcnonlinear.ami_kde import ami_kde
from quarcnonlinear.ami_knn import ami_knn
from quarcnonlinear.ami_partition import ami_partition

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


def _arr(v):
    """A reference vector; jsonencode writes NaN/Inf as null."""
    return np.array([np.nan if e is None else e for e in v], dtype=np.float64)


def _same_tau(got, want):
    if want is None:
        return np.isnan(got)
    return got == want


@needs_fixtures
def test_ami_histogram_matches_matlab():
    r = ref()
    for s, want in r["ami_histogram"].items():
        x = fixture(s)
        with _quiet():
            tau, curve, info = ami_histogram(x, 20)
            tau16, curve16, _ = ami_histogram(x, 20, bins=16)
        assert info["bins"] == want["bins"], s
        assert curve[:, 1] == pytest.approx(_arr(want["curve"]), rel=TOL), s
        assert _same_tau(tau, want["tau"]), s
        w16 = r["ami_histogram_bins16"][s]
        assert curve16[:, 1] == pytest.approx(_arr(w16["curve"]), rel=TOL), s
        assert _same_tau(tau16, w16["tau"]), s


@needs_fixtures
def test_ami_kde_matches_matlab():
    r = ref()
    for s, want in r["ami_kde"].items():
        tau, curve, _ = ami_kde(fixture(s), len(want["curve"]) - 1)
        assert curve[:, 1] == pytest.approx(_arr(want["curve"]), rel=TOL), s
        assert _same_tau(tau, want["tau"]), s
    # The older single-value reference from the original AMI_Thomas.
    assert ami(fixture("lorenz_2048"), 20, algorithm="kde")[0] == r["ami_thomas"]["lorenz_2048"]


@needs_fixtures
def test_ami_kde_adaptive_matches_matlab():
    for s, want in ref()["ami_kde_adaptive"].items():
        tau, curve, _ = ami_kde(fixture(s), len(want["curve"]) - 1, bandwidth="adaptive")
        assert curve[:, 1] == pytest.approx(_arr(want["curve"]), rel=TOL), s
        assert _same_tau(tau, want["tau"]), s


@needs_fixtures
def test_ami_partition_knn_copula_match_matlab():
    r = ref()
    cases = {
        "ami_partition": lambda x: ami_partition(x, 15),
        "ami_partition_a01": lambda x: ami_partition(x, 15, alpha=0.01),
        "ami_knn": lambda x: ami_knn(x, 10),
        "ami_knn_k3": lambda x: ami_knn(x, 6, k=3),
        "ami_copula": lambda x: ami_copula(x, 30),
    }
    for key, f in cases.items():
        for s, want in r[key].items():
            with _quiet():
                tau, curve, _ = f(fixture(s))
            assert curve[:, 1] == pytest.approx(_arr(want["curve"]), rel=TOL), (key, s)
            assert _same_tau(tau, want["tau"]), (key, s)


def test_ami_refuses_missing_data_and_impossible_lags():
    x = np.sin(np.arange(1, 201) / 5)
    x[16] = np.nan
    with pytest.raises(ValueError, match="contains NA"):
        ami_histogram(x, 10)
    with pytest.raises(ValueError, match="contains NA"):
        ami_kde(x, 10)
    with pytest.raises(ValueError, match="must be smaller"):
        ami_histogram(np.random.default_rng(1).standard_normal(50), 60)


def test_ami_dispatches_aliases_and_rejects_unknown():
    t = np.arange(1, 601)
    x = np.sin(t / 7) + 0.1 * np.cos(t / 3)
    assert ami(x, 10, algorithm="stergiou")[2]["estimator"] == "histogram"
    assert ami(x, 10, algorithm="thomas")[2]["estimator"] == "kde"
    with pytest.raises(ValueError, match="Unknown algorithm"):
        ami(x, 10, algorithm="wavelet")
    with pytest.warns(UserWarning, match="ignored"):
        ami(x, 10, algorithm="kde", bins=8)


def test_ami_returns_scalar_tau_and_fixed_shape_curve():
    x = np.sin(np.arange(1, 601) / 7)
    tau, curve, _ = ami(x, 15)
    assert np.isscalar(tau)
    assert curve.shape == (16, 2)
    assert tau in curve[:, 0]


def test_new_ami_estimators_dispatch_with_aliases():
    t = np.arange(1, 601)
    x = np.sin(t / 7) + 0.1 * np.cos(t / 3)
    pairs = {"kde-adaptive": "kde-adaptive", "adaptive": "adaptive", "partition": "adaptive",
             "knn": "knn", "ksg": "knn", "copula": "copula", "gcmi": "copula"}
    for a, est in pairs.items():
        with _quiet():
            assert ami(x, 8, algorithm=a)[2]["estimator"] == est, a
    assert ami(x, 8, algorithm="knn", k=5)[2]["k"] == 5
    with pytest.warns(UserWarning, match="ignored"):
        ami(x, 8, algorithm="histogram", k=5)


def test_rank_based_estimators_ignore_monotone_transforms():
    t = np.arange(1, 801)
    x = np.sin(t / 9) + 0.3 * np.sin(t / 2.3)
    for f in (ami_partition, ami_copula):
        with _quiet():
            a = f(np.exp(x), 10)[1]
            b = f(x, 10)[1]
        np.testing.assert_allclose(a, b, rtol=1e-12)


def test_ami_knn_leaves_global_random_state_alone():
    np.random.seed(42)
    before = np.random.get_state()[1].copy()
    t = np.arange(1, 301)
    with _quiet():
        ami_knn(np.sin(t / 5) + 0.2 * np.cos(t / 1.7), 15)
    assert np.array_equal(np.random.get_state()[1], before)


def test_ami_histogram_accepts_bin_edges():
    x = np.cumsum(np.random.default_rng(3).standard_normal(800))
    with _quiet():
        _, c24, i24 = ami(x, 30, bins=24)
        _, ce, ie = ami(x, 30, bins=i24["edges"])
        np.testing.assert_array_equal(ce, c24)
        assert ie["bins"] == 24
        q = np.quantile(x, np.linspace(0, 1, 25))   # R's default type 7
        assert np.max(np.abs(ami(x, 30, bins=q)[1][:, 1] - c24[:, 1])) > 1e-3
    with pytest.raises(ValueError, match="do not span"):
        ami(x, 30, bins=q[1:])
    with pytest.raises(ValueError, match="strictly increasing"):
        ami(x, 30, bins=[1, 1, 2])
