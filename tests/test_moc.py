"""The MOC solver (step B2) and the C5G7 model."""

import numpy as np
import pytest

from fission.model import MODELS_DIR, ModelError, load_model, parse
from fission.transport import Geometry, cyclic_tracks
from fission.transport.moc import Polar, build, solve
from fission.transport.run import fuel_pins


def _pin_model(fuel="uo2", moderator="moderator", rings=None):
    return parse(
        {
            "model": {"format": 1, "name": "pin", "multigroup": "c5g7-xs.toml"},
            "materials": {
                "uo2": {"macroscopic": "UO2"},
                "mox": {"macroscopic": "MOX87"},
                "moderator": {"macroscopic": "moderator"},
            },
            "pins": {"p": {"radii": [0.54], "materials": [fuel, moderator]}},
            "lattices": {"c": {"pitch": 1.26, "legend": {"P": "p"}, "map": "P"}},
        },
        MODELS_DIR,
    )


def _problem(model, azim=8, spacing=0.1, boundary="reflective", rings=None, n_polar=2):
    g = Geometry.from_lattice(model, "c", rings=rings, sectors=4)
    return build(
        model, g, cyclic_tracks(g.width, g.height, azim, spacing, boundary), n_polar
    )


def _k_infinity(xs):
    a = np.diag(xs.total) - xs.scatter.T
    return max(
        abs(np.linalg.eigvals(np.linalg.solve(a, np.outer(xs.chi, xs.nu_fission))))
    )


def test_polar_quadrature():
    q = Polar.gauss(3)
    assert q.weight.sum() == pytest.approx(1.0, rel=1e-15)
    mu = np.sqrt(1 - q.sin**2)
    assert np.all((0 < mu) & (mu < 1))
    assert q.weight @ mu == pytest.approx(0.5, rel=1e-14)  # mean cosine, exact
    assert q.weight @ mu**4 == pytest.approx(0.2, rel=1e-14)  # degree 4 < 2n


def test_one_material_everywhere_gives_k_infinity():
    """With one material in every region the flat flux is exact, whatever the
    tracking: MOC must reproduce the infinite-medium eigenvalue."""
    model = _pin_model("uo2", "uo2")
    sol = solve(_problem(model), tol_k=1e-11, tol_source=1e-10)
    assert sol.k == pytest.approx(_k_infinity(model.multigroup["UO2"]), rel=1e-9)
    spread = np.ptp(sol.flux, axis=0) / sol.flux.max(axis=0)
    assert spread.max() < 1e-9


def test_rust_sweep_matches_python():
    model = _pin_model()
    problem = _problem(model, rings={"uo2": 2})
    py = solve(problem, use_rust=False, max_iterations=25)
    rs = solve(problem, use_rust=True, max_iterations=25)
    assert rs.k == pytest.approx(py.k, rel=1e-13)
    np.testing.assert_allclose(rs.flux, py.flux, rtol=1e-12)
    assert rs.leakage == py.leakage == 0.0


@pytest.mark.parametrize("boundary", ["reflective", "vacuum"])
def test_neutron_balance(boundary):
    """At convergence k = production / (absorption + leakage)."""
    model = _pin_model("mox", "moderator")
    sol = solve(_problem(model, boundary=boundary), tol_k=1e-10, tol_source=1e-9)
    assert sol.balance() == pytest.approx(sol.k, rel=1e-7)
    if boundary == "vacuum":
        assert sol.leakage > 0 and sol.k < 0.5  # one bare pin leaks most neutrons
    else:
        assert sol.leakage == 0.0 and sol.k > 1.0


def test_pin_cell_converges_with_tracking():
    model = _pin_model()
    ks = [solve(_problem(model, azim, spacing, rings={"uo2": 3}, n_polar=3),
                tol_k=1e-9, tol_source=1e-8).k
          for azim, spacing in [(8, 0.1), (32, 0.02), (64, 0.01)]]  # fmt: skip
    assert abs(ks[2] - ks[1]) < abs(ks[1] - ks[0])
    assert abs(ks[2] - ks[1]) < 1e-3


# ---------------------------------------------------------------- C5G7
@pytest.fixture(scope="module")
def c5g7():
    return load_model("c5g7")


def test_c5g7_cross_sections_are_consistent(c5g7):
    for name, xs in c5g7.multigroup.items():
        # transport = absorption + out-scattering, to the 6 printed digits
        rebuilt = xs.absorption + xs.scatter.sum(axis=1)
        assert rebuilt == pytest.approx(xs.total, rel=5e-6), name
        if xs.fissile:
            assert xs.chi.sum() == pytest.approx(1.0, rel=1e-14)
            assert np.all(xs.fission <= xs.absorption)


def test_c5g7_layout(c5g7):
    g = Geometry.from_core(c5g7, sectors=1)
    assert (g.nx, g.ny, g.width) == (51, 51, pytest.approx(64.26))
    assert g.areas().sum() == pytest.approx(64.26**2, rel=1e-13)
    fuel = fuel_pins(c5g7, g)
    assert fuel.sum() == 4 * 264
    centre = [c.materials[0] for c in g.cells]
    assert centre.count("guide") == 4 * 24 and centre.count("chamber") == 4
    assert centre.count("mox43") == 2 * 64 and centre.count("mox87") == 2 * 100
    materials = g.materials()
    top = 64.26
    assert (
        materials[g.region(8.5 * 1.26, top - 8.5 * 1.26)] == "chamber"
    )  # inner UO2 centre
    assert materials[g.region(0.63, top - 21.42 - 0.63)] == "mox43"  # MOX below it
    assert materials[g.region(60.0, 5.0)] == "moderator"  # reflector corner
    assert c5g7.core is not None
    assert c5g7.core.boundary == {"left": "reflective", "right": "vacuum",
                                  "bottom": "vacuum", "top": "reflective"}  # fmt: skip


def test_model_options_for_macroscopic_data():
    base = {"model": {"format": 1, "name": "t", "multigroup": "c5g7-xs.toml"},
            "materials": {"w": {"macroscopic": "moderator"}},
            "pins": {"water": {"radii": [], "materials": ["w"]}},
            "lattices": {"r": {"pitch": 1.26, "fill": "water", "size": 3}}}  # fmt: skip
    m = parse(base, MODELS_DIR)
    assert m.lattices["r"].pins == (("water",) * 3,) * 3
    bad = {**base, "materials": {"w": {"macroscopic": "steel"}}}
    with pytest.raises(ModelError, match="'steel' is not in the multigroup library"):
        parse(bad, MODELS_DIR)
    bad = {**base, "materials": {"w": {"macroscopic": "moderator", "density": 1.0}}}
    with pytest.raises(ModelError, match="unknown key"):
        parse(bad, MODELS_DIR)
