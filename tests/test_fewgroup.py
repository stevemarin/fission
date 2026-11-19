"""Homogenization, B1 and condensation (step C3)."""

import math

import numpy as np
import pytest

from fission.data import EndfLibrary
from fission.fewgroup import FineConstants, condense, critical, mode
from fission.fewgroup.b1 import alpha
from fission.processing.groups import GroupStructure

HAVE_DATA = (
    EndfLibrary().root / "archives" / "neutrons-version.VIII.1.tar.gz"
).exists()


def test_alpha_limits_and_series():
    sigma = np.array([1.0])
    assert alpha(0.0, sigma)[0] == 1.0
    # just above the switch (x = 1e-3) the exact forms agree with the series
    # 1 +- 4 x^2 / 15 (truncation ~x^4; the exact forms lose ~1e-10 to the
    # cancellation in x - arctan x)
    x = 1.001e-3
    for sign in (1, -1):
        exact = alpha(sign * x**2, sigma)[0]
        assert exact == pytest.approx(1 + sign * 4 * x**2 / 15, abs=1e-9)
    x = 0.5
    expected = x**2 * math.atan(x) / (3 * (x - math.atan(x)))
    assert alpha(x**2, sigma)[0] == pytest.approx(expected, rel=1e-14)
    assert alpha(x**2, sigma)[0] > 1 > alpha(-(x**2), sigma)[0]
    with pytest.raises(ValueError, match="no real mode"):
        alpha(-4.0, sigma)


def _constants(groups=3, seed=1) -> FineConstants:
    """A random but physical fast-to-thermal medium (downscatter, fission in
    the last group, forward-peaked P1 scattering)."""
    rng = np.random.default_rng(seed)
    edges = tuple(np.geomspace(2e7, 1e-5, groups + 1))
    total = rng.uniform(0.3, 1.5, groups)
    s0 = np.zeros((groups, groups))
    for g in range(groups):
        s0[g, g:] = rng.uniform(0.05, 0.3, groups - g)
        s0[g] *= 0.9 * total[g] / s0[g].sum()  # 10% absorption
    s1 = 0.2 * s0
    nu = np.zeros(groups)
    nu[-1] = 0.15
    nu[0] = 0.005
    chi = np.zeros(groups)
    chi[0] = 1.0
    return FineConstants(
        "test", GroupStructure("test", edges), np.ones(groups), total,
        np.array([s0, s1]), nu, nu / 2.4, nu / 2.4 * 2e8, chi,
    )  # fmt: skip


def test_one_group_b1_closed_form():
    """One group: k(B^2) = nu Sigma_f / (Sigma_a + B^2 / (3 (alpha Sigma_t - Sigma_s1)))."""
    c = FineConstants(
        "one", GroupStructure("one", (2e7, 1e-5)), np.ones(1), np.array([0.8]),
        np.array([[[0.7]], [[0.1]]]), np.array([0.12]), np.array([0.05]),
        np.array([1e7]), np.array([1.0]),
    )  # fmt: skip
    for b2 in (0.0, 0.01, -0.01):
        a = alpha(b2, c.total)[0]
        expected = 0.12 / (0.1 + b2 / (3 * (a * 0.8 - 0.1)))
        assert mode(c, b2).k == pytest.approx(expected, rel=1e-13)


def test_critical_buckling_and_condensation():
    c = _constants(6)
    k_inf = c.k_infinity()
    assert mode(c, 0.0).k == pytest.approx(k_inf, rel=1e-13)
    crit = critical(c)
    assert k_inf > 1 and crit.b2 > 0
    assert crit.k == pytest.approx(1.0, abs=1e-11)
    edge = c.groups.boundaries[3]
    # the B1 spectrum: few-group diffusion is critical at the same buckling
    fg = condense(c, crit, (edge,))
    assert fg.k(crit.b2) == pytest.approx(1.0, abs=1e-11)
    # the infinite spectrum: few-group k-infinity = fine-group k-infinity
    fg0 = condense(c, mode(c, 0.0), (edge,))
    assert fg0.k(0.0) == pytest.approx(k_inf, rel=1e-12)
    # leakage hardens the spectrum: less weight in the thermal groups
    assert fg.flux[1] / fg.flux[0] < fg0.flux[1] / fg0.flux[0]
    with pytest.raises(ValueError, match="not a boundary"):
        condense(c, crit, (1.234,))


def test_diffusion_coefficient_p1_limit():
    """B -> 0 with a diagonal P1 matrix: D_g = 1 / (3 (Sigma_t - Sigma_s1,gg))."""
    c = _constants(4)
    c.scatter[1] = np.diag(np.diag(c.scatter[1]))
    d = mode(c, 1e-12).diffusion
    expected = 1 / (3 * (c.total - np.diag(c.scatter[1])))
    np.testing.assert_allclose(d, expected, rtol=1e-9)


@pytest.mark.skipif(not HAVE_DATA, reason="ENDF/B-VIII.1 neutron data not downloaded")
def test_fission_energy_release():
    from fission.fewgroup.kappa import COMPONENTS, recoverable_energy

    u235 = recoverable_energy("U235")  # LFC=1: thermal LIST + TAB1 components
    assert u235.thermal[COMPONENTS.index("ER")] == pytest.approx(1.940770e8)
    assert u235(0.0253) == pytest.approx(1.940770e8, rel=1e-4)
    t = u235.thermal  # ER = ET - ENU by definition (ENDF-102 sec. 1.5.1)
    assert t[7] == pytest.approx(t[8] - t[6], rel=1e-6)
    u234 = recoverable_energy("U234")  # LFC=0, NPLY=0: constant
    assert u234(1e6) == u234(0.0253)


@pytest.mark.skipif(not HAVE_DATA, reason="ENDF/B-VIII.1 neutron data not downloaded")
def test_homogenized_cell_keeps_the_lattice_eigenvalue():
    """Flux-volume weighting reproduces the MOC k-infinity of a small cell."""
    from fission.fewgroup import homogenize
    from fission.model import parse
    from fission.pincell.heterogeneous import Result, macroscopic
    from fission.processing.groupr import GroupLibrary
    from fission.transport import Geometry, cyclic_tracks
    from fission.transport.moc import solve, tracking, with_materials

    model = parse({
        "model": {"format": 1, "name": "pin"},
        "materials": {m: {"temperature": "moderator", "atom_densities": {"H1": 1.0}}
                      for m in ("fuel", "water")},
        "pins": {"p": {"radii": [0.45], "materials": ["fuel", "water"]}},
        "lattices": {"c": {"pitch": 1.26, "legend": {"P": "p"}, "map": "P"}},
    })  # fmt: skip
    g = Geometry.from_lattice(model, "c", rings={"fuel": 2}, sectors=4)
    p = tracking(g, cyclic_tracks(g.width, g.height, 16, 0.05), 2)
    groups = GroupStructure("two", (2e7, 0.625, 1e-5))

    def lib(name, total, s0, s1, fission=None):
        s = np.array([s0, s1])
        sigma = {1: np.array(total)}
        extra = {}
        if fission is not None:
            sigma[18] = np.array(fission)
            extra = {
                "nu_fission": 2.43 * np.array(fission),
                "chi": np.array([1.0, 0.0]),
            }
        return GroupLibrary(name, 0.0, groups, "", 1, np.ones(2), sigma, s,
                            np.zeros_like(s), **extra)  # fmt: skip

    libraries = {
        "fuel": {"U235": lib("U235", [10.0, 60.0], [[8.0, 0.2], [0.0, 10.0]],
                             [[1.0, 0.0], [0.0, 1.0]], [1.0, 40.0])},
        "water": {"H1": lib("H1", [20.0, 40.0], [[15.0, 4.9], [0.1, 39.5]],
                            [[8.0, 1.0], [0.0, 20.0]])},
    }  # fmt: skip
    densities = {"fuel": {"U235": 0.01}, "water": {"H1": 0.05}}
    xs = {m: macroscopic(m, [(densities[m][n], li) for n, li in libs.items()])
          for m, libs in libraries.items()}  # fmt: skip
    sol = solve(with_materials(p, xs), tol_k=1e-10, tol_source=1e-9)
    r = Result(None, "none", sol, {}, {}, densities, libraries)  # type: ignore[arg-type]
    c = homogenize(r)
    assert c.k_infinity() == pytest.approx(sol.k, rel=1e-7)
    np.testing.assert_allclose(mode(c, 0.0).flux / c.flux,
                               (mode(c, 0.0).flux / c.flux)[0], rtol=1e-6)  # fmt: skip
