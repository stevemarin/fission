"""The state matrix (M2): branches in the model file and the driver."""

import numpy as np
import pytest
from test_model import toy

from fission import branches
from fission.fewgroup.homogenize import FineConstants
from fission.model import BRANCH_AXES, ModelError, parse
from fission.processing.groups import GroupStructure

MATRIX = {
    "base": "hzp",
    "fuel_temperature": [565.0, 900.0, 1200.0],
    "boron_ppm": [0.0, 1300.0],
}


def test_branches_are_the_tensor_product_of_the_axes():
    m = parse(toy(branches={"doppler": MATRIX}))
    b = m.branches["doppler"]
    assert b.shape == (3, 1, 1, 2)
    states = b.states()
    assert len(states) == 6
    # C order, boron fastest; unlisted axes keep the base state's value
    assert [(s.fuel_temperature, s.boron_ppm) for s in states[:3]] == [
        (565.0, 0.0), (565.0, 1300.0), (900.0, 0.0)
    ]  # fmt: skip
    assert {s.moderator_density for s in states} == {0.743}
    assert states[-1].name == "doppler-2-0-0-1"


@pytest.mark.parametrize(
    ("patch", "message"),
    [
        ({"base": "hfp"}, "branches.b.base: unknown state 'hfp'"),
        ({"boron_ppm": [1300.0, 0.0]}, "values must increase"),
        ({"boron_ppm": [0.0, 2000.0]}, "must contain the base state's value 1300"),
        ({"boron_ppm": []}, "expected a non-empty list"),
        ({"fuel_temperature": [-1.0, 565.0]}, "values must be positive"),
        ({"fuel_temp": [565.0]}, "unknown key(s) ['fuel_temp']"),
        ({"lattice": "B"}, "unknown lattice 'B'"),
    ],
)  # fmt: skip
def test_branch_mistakes_are_reported(patch, message):
    with pytest.raises(ModelError) as err:
        parse(toy(branches={"b": {**MATRIX, **patch}}))
    assert message in str(err.value)


def _medium(state) -> tuple[float, FineConstants]:
    """A two-group stand-in for a lattice solver: thermal absorption grows
    with boron, resonance absorption with the fuel temperature."""
    t, ppm = state.fuel_temperature, state.boron_ppm
    total = np.array([0.6, 1.4])
    s0 = np.array([[0.55 - 1e-5 * t, 0.03], [0.0, 1.3 - 2e-5 * ppm]])
    nu = np.array([0.006, 0.16])
    c = FineConstants(
        state.name, GroupStructure("two", (2e7, 0.625, 1e-5)), np.ones(2), total,
        np.array([s0, 0.1 * s0]), nu, nu / 2.4, nu / 2.4 * 2e8, np.array([1.0, 0.0]),
    )  # fmt: skip
    return c.k_infinity(), c


@pytest.fixture
def fake_solver(monkeypatch):
    calls = []

    def solve(model, lattice, state, solver, options):
        calls.append(state.name)
        return _medium(state)

    monkeypatch.setattr(branches, "lattice_constants", solve)
    monkeypatch.setattr(branches, "_drop_shaped_libraries", lambda name: None)
    return calls


def test_driver_fills_the_grid_and_resumes(tmp_path, fake_solver):
    m = parse(toy(branches={"doppler": MATRIX}))
    path = tmp_path / "table.json"
    table = branches.run_matrix(m, "doppler", cuts=(0.625,), path=path)
    assert len(fake_solver) == 6 and table.complete
    k = table.array("k_infinity")
    assert k.shape == (3, 1, 1, 2)
    assert table.array("D").shape == (3, 1, 1, 2, 2)
    assert table.array("sigma_s").shape == (3, 1, 1, 2, 2, 2)
    assert k[0, 0, 0, 0] > k[2, 0, 0, 0] and k[0, 0, 0, 0] > k[0, 0, 0, 1]
    # B1 critical spectrum: B^2 > 0 where k-infinity > 1, < 0 where below
    for rec in table.points.values():
        assert (rec["k_infinity"] - 1) * rec["buckling_cm-2"] > 0
    # a second run reads the table and solves nothing; --fresh solves all again
    again = branches.run_matrix(m, "doppler", cuts=(0.625,), path=path)
    assert len(fake_solver) == 6
    np.testing.assert_array_equal(again.array("k_infinity"), k)
    branches.run_matrix(m, "doppler", cuts=(0.625,), path=path, fresh=True)
    assert len(fake_solver) == 12
    # other settings: the old table is not reused
    branches.run_matrix(m, "doppler", spectrum="infinite", cuts=(0.625,), path=path)
    assert len(fake_solver) == 18


def test_coefficients_are_differences_around_the_base(fake_solver):
    m = parse(toy(branches={"doppler": MATRIX}))
    table = branches.run_matrix(m, "doppler", cuts=(0.625,))
    rho = 1e5 * (1 - 1 / table.array("k_infinity"))
    found = {c.axis: c for c in branches.coefficients(table)}
    assert set(found) == {"fuel_temperature", "boron_ppm"}
    # base: 565 K (first node, one-sided), 1300 ppm (last node, one-sided)
    doppler = (rho[1, 0, 0, 1] - rho[0, 0, 0, 1]) / (900.0 - 565.0)
    boron = (rho[0, 0, 0, 1] - rho[0, 0, 0, 0]) / 1300.0
    assert found["fuel_temperature"].value == pytest.approx(doppler, rel=1e-12)
    assert found["boron_ppm"].value == pytest.approx(boron, rel=1e-12)
    assert doppler < 0 and boron < 0
    assert set(BRANCH_AXES) >= set(found)


def test_infinite_medium_constants_keep_the_eigenvalue():
    from fission.data import EndfLibrary

    if not (
        EndfLibrary().root / "archives" / "neutrons-version.VIII.1.tar.gz"
    ).exists():
        pytest.skip("ENDF/B-VIII.1 neutron data not downloaded")
    from fission.fewgroup.homogenize import from_mixture
    from fission.pincell.run import run

    sol = run("1A", "none")
    c = from_mixture(sol, "1A")
    assert c.k_infinity() == pytest.approx(sol.k, rel=1e-12)
    assert c.chi.sum() == pytest.approx(1.0, rel=1e-9)
