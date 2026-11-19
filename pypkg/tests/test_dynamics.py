"""Dynamics, phase, surrogates: against MATLAB's answers and by contract."""

import json
import warnings
from contextlib import contextmanager
from pathlib import Path

import numpy as np
import pytest

from quarcnonlinear.chaos_library import chaos_library
from quarcnonlinear.dfa import dfa
from quarcnonlinear.embed import embed
from quarcnonlinear.fgn_sim import fgn_sim
from quarcnonlinear.hilbert import hilbert
from quarcnonlinear.lyapunov import lyapunov
from quarcnonlinear.lye_w import lye_w
from quarcnonlinear.ode23 import ode23
from quarcnonlinear.ode45 import ode45
from quarcnonlinear.psr import psr
from quarcnonlinear.rel_phase_cont import rel_phase_cont
from quarcnonlinear.rel_phase_disc import rel_phase_disc
from quarcnonlinear.surr_find_rho import surr_find_rho
from quarcnonlinear.surr_pseudo_periodic import surr_pseudo_periodic
from quarcnonlinear.surr_theiler import surr_theiler

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


def _num(v):
    """JSON list to float array; jsonencode writes NaN/Inf as null."""
    return np.array([np.nan if e is None else e for e in np.atleast_1d(v)], dtype=np.float64)


# ---------------------------------------------------------------- MATLAB agreement

@needs_fixtures
def test_dfa_matches_matlab_at_every_scale():
    for name, r in ref()["dfa"].items():
        _, fl, alpha = dfa(fixture(name), r["scales"], 1)
        assert fl == pytest.approx(r["fluctuation"], rel=TOL), name
        assert alpha == pytest.approx(r["alpha"], rel=TOL), name


@needs_fixtures
def test_lye_w_matches_matlab_every_replacement_step():
    r = ref()["lye_w"]
    lz = fixture("lorenz_2048")
    cases = {"lorenz_e10": lye_w(lz, 100, 10, 3, 10),
             "henon_e1": lye_w(fixture("henon_512"), 1, 1, 2, 1),
             "lorenz_custom": lye_w(lz[:1500], 100, 8, 4, 5, 2, 20 * np.pi / 180, 1)}
    for name, (out, LyE) in cases.items():
        want = _num(r[name]["out"]).reshape(-1, r[name]["ncol"])
        want[:, 1:3] -= np.where(want[:, 0:1] > 0, 1, 0)   # point, neighbour 0-based
        assert out.shape == want.shape, name
        assert out == pytest.approx(want, rel=TOL), name
        assert LyE == pytest.approx(r[name]["LyE"], rel=TOL), name


@needs_fixtures
def test_lyapunov_matches_matlab_for_both_methods():
    r = ref()["lyapunov"]
    lz = fixture("lorenz_2048")
    lam, extra = lyapunov(lz, 100, delay=10, dim=3)
    assert lam == pytest.approx(r["lorenz_rosenstein"]["lambda"], rel=TOL)
    assert np.array_equal(extra["scalingRegion"], np.array(r["lorenz_rosenstein"]["scalingRegion"]) - 1)
    assert extra["fitR2"] == pytest.approx(r["lorenz_rosenstein"]["fitR2"], rel=TOL)
    lam, extra = lyapunov(fixture("henon_512"), 1, delay=1, dim=2)
    assert lam == pytest.approx(r["henon_rosenstein"]["lambda"], rel=TOL)
    assert np.array_equal(extra["scalingRegion"], np.array(r["henon_rosenstein"]["scalingRegion"]) - 1)
    lam, _ = lyapunov(lz, 100, algorithm="wolf", delay=10, dim=3, evolve=10)
    assert lam == pytest.approx(r["lorenz_wolf"]["lambda"], rel=TOL)


@needs_fixtures
def test_rel_phase_cont_matches_matlab_hilbert_included():
    r = ref().get("rel_phase_cont")
    if r is None:
        pytest.skip("MATLAB reference made without the Signal Processing Toolbox")
    lz = fixture("lorenz_2048")
    crp, crpH = rel_phase_cont(lz[:800], lz[20:820], 100)
    assert crp == pytest.approx(_num(r["lorenz"]["crp"]), rel=TOL, nan_ok=True)
    assert crpH == pytest.approx(_num(r["lorenz"]["crpH"]), rel=TOL, nan_ok=True)


@needs_fixtures
def test_chaos_library_matches_matlab_maps_exactly_flows_to_rounding():
    r = ref()["chaos_library"]

    def want(nm, k):
        return np.array(r[nm]["y"], dtype=np.float64).reshape(-1, k)

    assert np.array_equal(chaos_library("Hennon", 500, [0.1, 0.1], [1.4, 0.3])[1], want("henon", 2))
    assert np.array_equal(chaos_library("Logistic", 300, 0.2, 3.9)[1], want("logistic", 1))
    # Flows: ode45/ode23 follow MATLAB's step control, so the steps coincide
    # and the trajectories agree to rounding over this short horizon.
    tt = np.arange(201) * 0.01
    flows = {"lorenz": ("Lorenz", [1, 1, 1], [10, 28, 8 / 3]),
             "rossler": ("Rossler", [1, 1, 0], [0.2, 0.2, 5.7]),
             "aizawa": ("Aizawa", [0.1, 0, 0], [0.7, 0.95, 0.6, 3.5, 0.25, 0.1]),
             "nosehoover": ("NoseHoover", [0, 5, 0], 1),
             "arneodo": ("Arneodo", [0.2, 0.2, -0.75], [-5.5, 3.5, -1])}
    for nm, (S, ic, p) in flows.items():
        t, y = chaos_library(S, tt, ic, p)
        assert t == pytest.approx(tt, abs=1e-15), nm
        assert np.max(np.abs(y - want(nm, 3))) < 1e-10, nm


@needs_fixtures
def test_rel_phase_disc_matches_matlab():
    r = ref().get("rel_phase_disc")
    if r is None:
        pytest.skip("MATLAB reference made without the Signal Processing Toolbox")
    tt = np.arange(2000)
    cases = {"ratio2": (np.sin(2 * np.pi * tt / 50),
                        np.sin(2 * np.pi * tt / 25 + 0.4 * np.sin(2 * np.pi * tt / 300))),
             "ratio1_lag": (np.sin(2 * np.pi * tt / 40),
                            np.sin(2 * np.pi * tt / 40 - 0.8 + 0.6 * np.sin(2 * np.pi * tt / 500))),
             "ratio1": (np.sin(2 * np.pi * tt / 40),
                        np.sin(2 * np.pi * tt / 40 + 0.8 + 0.6 * np.sin(2 * np.pi * tt / 500)))}
    for nm, (x, y) in cases.items():
        w = r[nm]
        FR, PC, DRP = rel_phase_disc(x, y)
        for f in ("t", "n", "T", "DRP"):
            assert DRP[f] == pytest.approx(_num(w[f]), rel=TOL), (nm, f)
        for f in ("Ratio", "Occurances", "Percentage_of_Total"):
            assert FR[f] == pytest.approx(_num(w[f]), rel=TOL), (nm, f)
        assert np.array_equal(PC["n"], _num(w["n_pc"])), nm
        assert PC["PC"] == pytest.approx(_num(w["PC"]), rel=TOL), nm


# ---------------------------------------------------------------- contracts

def test_lyapunov_requires_evolve_for_wolf_and_guards_inputs():
    i = np.arange(1, 801)
    x = np.sin(i / 7) + 0.3 * np.sin(i / 2.1)
    with pytest.raises(ValueError, match="needs evolve"):
        lyapunov(x, 1, algorithm="wolf", delay=5, dim=3)
    with pytest.raises(ValueError, match="contains NA"):
        lyapunov(np.append(x[:10], np.nan), 1)
    with pytest.raises(ValueError, match="Unknown algorithm"):
        lyapunov(x, 1, algorithm="ftle")
    with pytest.warns(UserWarning, match="without theiler"):
        lyapunov(psr(x, 5, 3), 1)
    assert lyapunov(x, 1, algorithm="w", delay=5, dim=3, evolve=3)[1]["estimator"] == "wolf"


def test_chaos_library_names_failures_and_solvers_reach_tfinal():
    with pytest.raises(ValueError, match="Unknown system"):
        chaos_library("Chua", np.arange(1, 11), 1, 1)
    with pytest.raises(ValueError, match="at least one"):
        chaos_library("Logistic", 0, 0.2, 3.9)
    t, y = ode45(lambda t, y: -y, [0, 2], 1, rtol=1e-8, atol=1e-10)
    assert t[-1] == 2
    assert y[-1, 0] == pytest.approx(np.exp(-2), rel=1e-7)
    ts = np.arange(11) * 0.1
    t3, y3 = ode23(lambda t, y: -y, ts, 1, rtol=1e-6, atol=1e-9)
    assert y3[:, 0] == pytest.approx(np.exp(-ts), rel=1e-5)


def test_embed_and_hilbert_follow_matlab():
    z = np.arange(1.0, 11.0)
    e = embed(z, [0, 2])
    assert e.shape == (2, 8)
    assert np.array_equal(e[0], z[2:]) and np.array_equal(e[1], z[:-2])
    assert np.array_equal(embed(z, 3, 2), embed(z, [0, 2, 4]))
    assert np.array_equal(embed(z, [0, 3]), psr(z, 3, 2)[:, ::-1].T)
    with pytest.raises(ValueError, match="too small"):
        embed(z[:3], [0, 5])
    for n in (64, 65):
        t = np.arange(n)
        x = np.cos(2 * np.pi * 4 * t / n)
        h = hilbert(x)
        assert np.real(h) == pytest.approx(x, abs=1e-12)
        assert np.imag(h) == pytest.approx(np.sin(2 * np.pi * 4 * t / n), abs=1e-12)


# ---------------------------------------------------------------- surrogates and fGn
# Random, and numpy cannot draw MATLAB's numbers, so tested by what each
# algorithm promises to preserve (as R's test-surrogates.R).

def _amp(v):
    return np.abs(np.fft.fft(v - v.mean()))


def test_algorithm1_keeps_amplitude_spectrum_exactly():
    rng = np.random.default_rng(1)
    i = np.arange(1, 513)
    x = np.cumsum(rng.standard_normal(512)) + np.sin(i / 5)
    for _ in range(5):
        z = surr_theiler(x, 1, rng=rng)
        assert _amp(z) == pytest.approx(_amp(x), rel=1e-9, abs=1e-9)
        assert z.mean() == pytest.approx(x.mean(), rel=1e-12, abs=1e-12)
        assert not np.allclose(z, x)


def test_algorithms_0_and_2_keep_distribution_exactly():
    rng = np.random.default_rng(2)
    x = rng.exponential(size=400) ** 2
    for alg in (0, 2):
        z = surr_theiler(x, alg, rng=rng)
        assert np.array_equal(np.sort(z), np.sort(x)), alg
        assert not np.array_equal(z, x), alg


def test_algorithm2_spectrum_close_not_exact():
    rng = np.random.default_rng(3)
    i = np.arange(1, 1025)
    x = np.sin(i / 8) + 0.3 * rng.standard_normal(1024)
    z = surr_theiler(x, 2, rng=rng)
    assert np.corrcoef(_amp(z)[1:200], _amp(x)[1:200])[0, 1] > 0.9


def test_surrogates_keep_shape_and_reject_unknown_algorithms():
    rng = np.random.default_rng(4)
    m = rng.standard_normal((100, 2))
    for alg in (0, 1, 2):
        assert surr_theiler(m, alg, rng=rng).shape == (100, 2)
    assert surr_theiler(rng.standard_normal(99), 1, rng=rng).shape == (99,)
    with pytest.raises(ValueError, match="Unknown algorithm"):
        surr_theiler(rng.standard_normal(50), 3)


def test_seeded_surrogates_repeat():
    x = np.sin(np.arange(200) / 5)
    assert np.array_equal(surr_theiler(x, 2, rng=9), surr_theiler(x, 2, rng=9))
    assert np.array_equal(fgn_sim(64, 0.7, rng=9), fgn_sim(64, 0.7, rng=9))


def test_pseudo_periodic_surrogates_copy_values_from_series():
    rng = np.random.default_rng(5)
    i = np.arange(1, 601)
    x = np.sin(2 * np.pi * i / 40) + 0.05 * rng.standard_normal(600)
    ys, yi = surr_pseudo_periodic(x, 10, 3, 0.2, rng=rng)
    assert ys.size == 600 - 20
    assert np.array_equal(ys, x[yi])
    assert np.all((yi >= 0) & (yi <= 579))


def test_surr_find_rho_returns_rho_inside_bounds():
    rng = np.random.default_rng(6)
    i = np.arange(1, 401)
    x = np.sin(2 * np.pi * i / 40) + 0.05 * rng.standard_normal(400)
    with _quiet():
        rho, out = surr_find_rho(x, 10, 3, 0.05, 1, rng=rng)
    assert 0.05 <= rho <= 1
    assert out.shape[1] == 3
    with pytest.raises(ValueError, match="0 < rhoL < rhoH"):
        surr_find_rho(x, 10, 3, 1, 0.5)


def test_fgn_sim_unit_variance_and_dfa_recovers_H():
    rng = np.random.default_rng(7)
    for H in (0.3, 0.5, 0.85):
        x = fgn_sim(8192, H, rng=rng)
        assert x.size == 8192
        assert np.std(x, ddof=1) == pytest.approx(1, rel=0.1)
        sc = np.unique(np.floor(2 ** np.arange(4, 10.25, 0.5) + 0.5))
        assert dfa(x, sc)[2] == pytest.approx(H, rel=0.08 / H), H
