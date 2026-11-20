"""Equivalence-theory self-shielding (step C1)."""

import math

import numpy as np
import pytest

from fission.data import EndfLibrary
from fission.model import parse
from fission.pincell.heterogeneous import macroscopic
from fission.pincell.shielding import (
    dancoff,
    escape_probability,
    potential_scattering,
    shielded,
    sigma0,
    surface,
)
from fission.processing.groupr import GroupLibrary
from fission.processing.groups import GroupStructure
from fission.processing.library import cached_groups
from fission.transport import Geometry, cyclic_tracks
from fission.transport.moc import tracking


def _cell(pitch=1.26):
    model = parse(
        {
            "model": {"format": 1, "name": "pin"},
            "materials": {
                m: {"temperature": "moderator", "atom_densities": {"H1": 1.0}}
                for m in ("fuel", "clad", "water")
            },
            "pins": {
                "p": {"radii": [0.41, 0.475], "materials": ["fuel", "clad", "water"]}
            },
            "lattices": {"c": {"pitch": pitch, "legend": {"P": "p"}, "map": "P"}},
        }
    )
    return Geometry.from_lattice(model, "c", rings={"fuel": 2}, sectors=8)


def _tracking(g, azim=32, spacing=0.02):
    return tracking(g, cyclic_tracks(g.width, g.height, azim, spacing), 3, faces=False)


def test_surfaces_of_rings():
    g = _cell()
    assert surface(g, "fuel") == pytest.approx(2 * math.pi * 0.41)
    assert surface(g, "clad") == pytest.approx(2 * math.pi * (0.41 + 0.475))
    with pytest.raises(ValueError, match="fills the cell"):
        surface(g, "water")


def test_dancoff_limits():
    g = _cell()
    p = _tracking(g)
    sigma = {"fuel": 0.4, "clad": 0.3, "water": 1.1}
    lattice = dancoff(p, sigma, "fuel")
    assert 0.1 < lattice.dancoff < 0.6
    assert lattice.chord == pytest.approx(2 * 0.41, rel=1e-12)  # 4 V / S of a cylinder
    assert lattice.sigma_e == pytest.approx((1 - lattice.dancoff) / lattice.chord)
    # a transparent moderator: every neutron leaving a pin enters another one
    clear = dancoff(p, {"fuel": 0.4, "clad": 1e-12, "water": 1e-12}, "fuel")
    assert clear.dancoff == pytest.approx(1.0, abs=1e-6)
    # a pin alone in thick water: no shadowing
    far = _cell(pitch=12.6)
    alone = dancoff(_tracking(far, 32, 0.05), sigma, "fuel")
    assert abs(alone.dancoff) < 1e-4  # the tracking error cancels in the ratio
    # denser water shadows less
    dense = dancoff(p, {"fuel": 0.4, "clad": 0.3, "water": 2.0}, "fuel")
    assert dense.dancoff < lattice.dancoff


def test_background_cross_section():
    data = pytest.importorskip("fission.data")
    if not (
        data.EndfLibrary().root / "archives" / "neutrons-version.VIII.1.tar.gz"
    ).exists():
        pytest.skip("ENDF/B-VIII.1 neutron data not downloaded")
    # H-1: 4 pi a^2 with the ENDF scattering radius is the free-atom value ~20.4 b
    assert potential_scattering("H1") == pytest.approx(20.44, abs=0.01)
    dens = {"U238": 0.02, "O16": 0.04}
    s0 = sigma0(dens, sigma_e=0.8)
    expected = (0.04 * potential_scattering("O16") + 0.8) / 0.02
    assert s0["U238"] == pytest.approx(expected, rel=1e-14)


@pytest.fixture(scope="module")
def zr91():
    if not (
        EndfLibrary().root / "archives" / "neutrons-version.VIII.1.tar.gz"
    ).exists():
        pytest.skip("ENDF/B-VIII.1 neutron data not downloaded")
    return cached_groups("Zr91", 293.6, thermal="free-gas")


def test_shielding_limits(zr91):
    huge = shielded("Zr91", 293.6, 1e12, "free-gas")
    for mt in (1, 2, 102):
        np.testing.assert_allclose(huge.sigma[mt], zr91.sigma[mt], rtol=1e-8)
    strong = shielded("Zr91", 293.6, 10.0, "free-gas")
    # a background of 10 b shields the large resonances: capture drops a lot
    # there. Where capture does not follow the total cross section it can
    # rise: in 454-677 eV Zr-91's total is almost all elastic (7.19 of 7.24 b)
    # and capture (0.05 b) sits where the total is low, which the weight
    # 1 / (sigma_t + sigma_0) favours (+3.5%).
    ratio = strong.sigma[102] / zr91.sigma[102]
    assert ratio.min() < 0.5
    assert ratio.max() < 1.05


def test_transport_correction_keeps_absorption():
    groups = GroupStructure("two", (2e7, 0.625, 1e-5))
    s = np.zeros((2, 2, 2))
    s[0] = [[0.5, 0.1], [0.0, 1.0]]
    s[1] = [[0.2, 0.02], [0.0, 0.1]]  # forward-peaked scattering
    lib = GroupLibrary(
        "X",
        0.0,
        groups,
        "",
        1,
        np.ones(2),
        {1: np.array([0.7, 1.3])},
        s,
        np.zeros_like(s),
    )
    xs = macroscopic("m", [(2.0, lib)])
    assert xs.total == pytest.approx(2 * (np.array([0.7, 1.3]) - [0.22, 0.1]))
    assert xs.absorption == pytest.approx(2 * (np.array([0.7, 1.3]) - s[0].sum(axis=1)))
    assert xs.scatter[0, 1] == pytest.approx(2 * 0.1)  # off-diagonal untouched


def test_escape_probability_limits_and_rational_form():
    g = _cell()
    p = _tracking(g)
    sigma = {"fuel": 0.4, "clad": 0.3, "water": 1.1}
    sig, esc = escape_probability(
        p, sigma, "fuel", np.array([1e-4, 0.1, 1.0, 10.0, 1e4])
    )
    assert esc[0] == pytest.approx(1.0, abs=1e-3)  # transparent fuel: all escape
    assert esc[-1] < 1e-3  # black fuel: none
    assert np.all(np.diff(esc) < 0)
    # for black fuel the escapes are surface events: P -> Sigma_e / Sigma with
    # Sigma_e = (1 - C) / chord, the same limit as the rational form
    e = dancoff(p, sigma, "fuel")
    assert esc[-1] * sig[-1] == pytest.approx(e.sigma_e, rel=0.02)
    # in between, the Wigner rational form underestimates escape
    # [Bell & Glasstone sec. 2.8b]
    rational = e.sigma_e / (sig + e.sigma_e)
    assert esc[2] > rational[2]


# ---------------------------------------------------------------- slowing down
def _grid(e_top=2e3, e_low=1.0, n=60000, history=0.3):
    e = np.geomspace(e_top * np.exp(history), e_low, n)
    u = np.log(e[0] / e)
    return e, u, int(np.searchsorted(u, history)) + 1


def _lorentz(e, e0, width, peak):
    return peak / (1 + ((e - e0) / (width / 2)) ** 2)


def test_march_conserves_neutrons_in_pure_scattering():
    from fission.pincell.shielding import march

    e, _, start = _grid()
    for alphas in ([0.983], [0.983, 0.779]):
        sc = [np.full_like(e, 0.3) for _ in alphas]
        psi = march(e, sum(sc), sc, alphas, np.zeros_like(e), np.ones(start), start)
        assert np.abs(psi - 1).max() < 1e-5


def test_march_absorption_decay():
    """Constant absorption in a heavy scatterer: psi ~ exp(-Sigma_a u / (xi Sigma_t))."""
    from fission.pincell.shielding import march

    e, u, start = _grid()
    a = 0.983
    xi = 1 + a * math.log(a) / (1 - a)
    sa, ss = 0.003, 0.3
    psi = march(e, np.full_like(e, sa + ss), [np.full_like(e, ss)], [a], np.zeros_like(e),
                np.ones(start), start)  # fmt: skip
    rate = -np.polyfit(u[5000:], np.log(psi[5000:]), 1)[0]
    assert rate == pytest.approx(sa / (xi * (sa + ss)), rel=0.01)


def test_march_resonance_limits():
    """A light moderator plus an absorber resonance.

    * absorption resonance: the flux dip follows NR (Sigma_p / Sigma_t);
    * a pure scattering resonance of a very heavy nucleus: a neutron scattered
      in it stays at the same energy, so the flux does not dip at all (the
      wide-resonance limit) while NR would predict a deep dip.
    """
    from fission.pincell.shielding import march

    e, _, start = _grid(history=0.8)  # the moderator's window is ln 2 wide
    moderator = np.full_like(e, 1.0)  # Sigma_s of a light moderator (alpha = 0.5)
    absorber_p = np.full_like(e, 0.1)  # absorber potential scattering
    res = _lorentz(e, 100.0, 0.5, 50.0)  # a 0.5 eV wide resonance at 100 eV
    peak = int(np.argmin(np.abs(e - 100.0)))
    # absorption resonance, absorber with A = 20 (narrow compared with its loss)
    total = moderator + absorber_p + res
    psi = march(e, total, [moderator, absorber_p], [0.5, ((19 / 21) ** 2)],
                np.zeros_like(e), np.ones(start), start)  # fmt: skip
    nr = (moderator + absorber_p) / total
    # NR to within 10%: NR takes the incoming source as exactly asymptotic,
    # while the resonance's own wings already lower it a little (the flux
    # just below the resonance is 0.85: its escape probability)
    assert psi[peak] == pytest.approx(nr[peak], rel=0.1)
    # scattering resonance of a very heavy nucleus (A = 1e6): a neutron makes
    # ~50 collisions in the resonance per moderator collision and drifts only
    # ~1e-4 in lethargy, far less than the resonance width (0.005): it stays
    alpha = ((1e6 - 1) / (1e6 + 1)) ** 2
    total = moderator + absorber_p + res
    psi = march(e, total, [moderator, absorber_p + res], [0.5, alpha],
                np.zeros_like(e), np.ones(start), start)  # fmt: skip
    assert psi[peak] > 0.95
    assert ((moderator + absorber_p) / total)[peak] < 0.05  # what NR would say
    # A = 1e4 drifts ~0.01 per moderator collision, twice the resonance
    # width: the flux dips again (neither NR nor WR)
    alpha = ((1e4 - 1) / (1e4 + 1)) ** 2
    psi = march(e, total, [moderator, absorber_p + res], [0.5, alpha],
                np.zeros_like(e), np.ones(start), start)  # fmt: skip
    assert 0.05 < psi[peak] < 0.9


def test_ring_escape_averages_to_the_region():
    """E_k by reciprocity, one sweep per Sigma: the volume average over the
    rings (equal areas here) is the escape probability of the whole region,
    computed directly with a source inside it; the rim escapes more."""
    from fission.pincell.shielding import ring_escape

    g = _cell()
    p = _tracking(g)
    sigma = {"fuel": 0.4, "clad": 0.3, "water": 1.1}
    sigmas = np.array([0.01, 0.5, 2.0, 20.0])
    _, whole = escape_probability(p, sigma, "fuel", sigmas)
    _, rings = ring_escape(p, sigma, "fuel", sigmas)
    assert rings.shape == (2, 4)
    np.testing.assert_allclose(rings.mean(axis=0), whole, rtol=2e-3, atol=1e-5)
    assert np.all(rings[1] > rings[0])  # rim (ring 1) above centre (ring 0)
    assert rings[1, -1] > 5 * rings[0, -1]  # black fuel: only the rim escapes


def test_split_rings_labels():
    g = _cell().split_rings("fuel")
    assert g.cells[0].materials == ("fuel/0", "fuel/1", "clad", "water")


def test_collision_matrix_escapes_match_reciprocity():
    """Two independent routes to E_k: forward sweeps from each ring
    (1 - sum_i P_k->i) and the single reciprocal sweep of ring_escape."""
    from fission.pincell.shielding import collision_matrix, ring_escape

    p = _tracking(_cell())
    sigma = {"fuel": 0.4, "clad": 0.3, "water": 1.1}
    sigmas = np.array([0.05, 1.0, 10.0])
    _, matrix, volumes = collision_matrix(p, sigma, "fuel", sigmas)
    _, rings = ring_escape(p, sigma, "fuel", sigmas)
    np.testing.assert_allclose(1 - matrix.sum(axis=2), rings.T, rtol=3e-3, atol=1e-5)
    np.testing.assert_allclose(volumes, volumes[0])  # equal-area rings
    # reciprocity V_j P_j->i = V_i P_i->j
    for m in matrix:
        np.testing.assert_allclose(m * volumes[:, None], (m * volumes[:, None]).T,
                                   rtol=5e-3)  # fmt: skip


def _sph_case(group_edges, resonance):
    """A fuel region with the NR flux of a pointwise Sigma_t(E) on a fine grid."""
    from fission.pincell.shielding import escape_probability

    p = _tracking(_cell())
    sigma_p = {"fuel": 0.4, "clad": 0.3, "water": 1.1}
    sigmas, esc = escape_probability(p, sigma_p, "fuel")
    e = np.geomspace(1.0, 100.0, 20001)
    total = 0.4 + (30.0 / (1 + ((e - 20.0) / 0.2) ** 2) if resonance else 0.0)
    total = np.broadcast_to(total, e.shape).astype(float)
    escape = np.interp(np.log(total), np.log(sigmas), esc)
    flux = escape + (1 - escape) * 0.4 / total
    details = {"energy": e, "total": total, "source": 0.4 / total, "level": np.ones_like(e),
                   "flux": flux}  # fmt: skip
    groups = GroupStructure("test", tuple(group_edges))
    u = np.log(e)

    def average(v, lo, hi):
        m = (e >= lo) & (e <= hi)
        return np.trapezoid(v[m], u[m]) / (u[m][-1] - u[m][0])

    sigma_t = np.array([average(total * flux, *groups.bounds(g))
                        / average(flux, *groups.bounds(g)) for g in range(groups.size)])  # fmt: skip
    return p, sigma_p, details, groups, sigma_t, (sigmas, esc), average


def _sph(p, sigma_p, details, groups, sigma_t):
    from fission.pincell.shielding import sph_factors

    references = {"fuel": (details["flux"], details["source"] * details["total"])}
    return sph_factors(p, sigma_p, np.array(p.geometry.materials()), details["energy"],
                       details["level"], references, {"fuel": sigma_t}, groups,
                       energy_range=(1.0, 100.0))  # fmt: skip


def test_sph_is_one_without_structure():
    p, sigma_p, details, groups, sigma_t, _, _ = _sph_case((100.0, 10.0, 1.0), False)
    mu = _sph(p, sigma_p, details, groups, sigma_t)
    np.testing.assert_allclose(mu["fuel"], 1.0, atol=1e-6)


def test_sph_preserves_the_collision_rate():
    """With one region the MOC iteration reduces to the scalar equation
    mu phi_MG(mu) = <phi>, phi_MG = P(mu Sigma) + (1 - P) Sigma_p / (mu Sigma);
    the resonance group gets mu < 1 (the flux-weighted Sigma_g alone
    overestimates the collision rate) and mu Sigma_g phi_MG = <Sigma_t phi>."""
    p, sigma_p, details, groups, sigma_t, (sigmas, esc), average = _sph_case(
        (100.0, 30.0, 10.0, 1.0), True
    )
    mu = _sph(p, sigma_p, details, groups, sigma_t)["fuel"]
    g = 1  # 10-30 eV holds the resonance
    assert mu[g] < 0.99
    np.testing.assert_allclose(mu[[0, 2]], 1.0, atol=2e-3)
    s = mu[g] * sigma_t[g]
    pg = np.interp(np.log(s), np.log(sigmas), esc)
    phi_mg = pg + (1 - pg) * 0.4 / s
    reference = average(details["total"] * details["flux"], *groups.bounds(g))
    assert s * phi_mg == pytest.approx(reference, rel=3e-3)  # P table vs MOC


def test_coupled_table_reciprocity_and_consistency():
    """Rings + clad: V_j Sigma_j P_j->i = V_i Sigma_i P_i->j with different
    cross sections in fuel and clad; with the clad at its outside value,
    collisions in the clad plus escape beyond it are the rings' escape."""
    from fission.pincell.shielding import coupled_collisions, ring_escape

    p = _tracking(_cell())
    sigma = {"fuel": 0.4, "clad": 0.3, "water": 1.1}
    t = coupled_collisions(p, sigma, np.array([0.1, 2.0]), np.array([0.3, 5.0]))
    assert t.labels == ["fuel/0", "fuel/1", "clad"]
    sig_r = np.array([2.0, 2.0, 5.0])
    m = t.volumes[:, None] * sig_r[:, None] * t.collide[1, 1]
    np.testing.assert_allclose(m, m.T, rtol=5e-3)
    _, rings = ring_escape(p, sigma, "fuel", np.array([0.1, 2.0]))
    for a in range(2):
        beyond = 1 - t.collide[a, 0, :2, :2].sum(axis=1)  # leaves the fuel
        np.testing.assert_allclose(beyond, rings[:, a], rtol=3e-3, atol=1e-5)
    # bilinear interpolation reproduces the nodes
    np.testing.assert_allclose(t.at(np.array([2.0]), np.array([0.3]))[0],
                               t.collide[1, 0], rtol=1e-12)  # fmt: skip
