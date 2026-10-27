"""The Python bindings against MATLAB's answers on the shared fixtures."""

import json
from pathlib import Path

import numpy as np
import pytest

import quarcnonlinear as qn

REPO = Path(__file__).resolve().parents[2]
FIX = REPO / "tests" / "fixtures"
TOL = 1e-10  # as rpkg's tests: summation order passes, a defect does not

needs_fixtures = pytest.mark.skipif(not (FIX / "matlab_reference.json").exists(),
                                    reason="shared fixtures not found (not in the repository)")


def ref():
    return json.loads((FIX / "matlab_reference.json").read_text())


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


@pytest.mark.skipif(not (REPO / "cpp").exists(), reason="cpp/ not found")
def test_core_copy_is_current():
    core, copy = REPO / "cpp", REPO / "pypkg" / "src" / "core"
    for f in (core / "src").glob("*.cpp"):
        assert (copy / f"core_{f.name}").read_text() == f.read_text(), f.name
    for f in (core / "include" / "quarc").glob("*.hpp"):
        assert (copy / "quarc" / f.name).read_text() == f.read_text(), f.name
