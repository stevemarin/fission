"""Tests for secondary distributions: nu, MF4, MF5, MF6, kinematics."""

import math
import random

import pytest

from fission.data import EndfLibrary
from fission.endf import Evaluation
from fission.endf.distributions.angular import CM, angular_distribution
from fission.endf.distributions.energy import Panel, energy_distribution
from fission.endf.distributions.incident import (
    IncidentGrid,
    interpolate_panels,
    interpolate_parameter,
)
from fission.endf.distributions.kinematics import (
    cm_energy,
    elastic_alpha,
    kalbach_mann,
    kalbach_mann_mean_cosine,
    mean_lethargy_gain,
    two_body,
)
from fission.endf.distributions.legendre import (
    legendre,
    legendre_series,
    moments_of_linear,
)
from fission.endf.distributions.mf6 import Continuum, TwoBody, energy_angle
from fission.endf.distributions.nu import (
    delayed_groups,
    multiplicity_distribution,
    nu_bar,
)


def _integrate(f, a, b, n=4000):
    h = (b - a) / n
    return h * math.fsum(f(a + (k + 0.5) * h) for k in range(n))


# ---------------------------------------------------------------- Legendre
def test_legendre_recurrence_matches_explicit_polynomials():
    for mu in (-0.9, -0.3, 0.0, 0.4, 1.0):
        p = legendre(3, mu)
        assert p[2] == pytest.approx(0.5 * (3 * mu**2 - 1))
        assert p[3] == pytest.approx(0.5 * (5 * mu**3 - 3 * mu))


def test_legendre_series_is_normalized_and_moments_invert_it():
    a = [1.0, 0.31, -0.05, 0.02]
    assert _integrate(lambda m: legendre_series(a, m), -1, 1) == pytest.approx(
        1.0, rel=1e-7
    )
    mu = [-1 + 2 * k / 400 for k in range(401)]
    f = [legendre_series(a, m) for m in mu]  # tabulated finely: moments ~ a
    assert moments_of_linear(mu, f, 3) == pytest.approx(a, abs=2e-5)


def test_moments_of_a_linear_distribution_are_exact():
    # f = (1 + b mu)/2 has a_1 = b/3 and no higher moments
    b = 0.6
    a = moments_of_linear([-1.0, 1.0], [(1 - b) / 2, (1 + b) / 2], 4)
    assert a == pytest.approx([1.0, b / 3, 0.0, 0.0, 0.0], abs=1e-15)


# ------------------------------------------------- incident-energy interpolation
def test_bracket_laws_and_edges():
    grid = IncidentGrid.from_tab2([1.0, 2.0, 2.0, 4.0], [2, 4], [1, 22])
    assert grid.bracket(0.5) == grid.bracket(1.0)  # clamped below
    assert grid.bracket(1.5).t == 0.0  # histogram region
    b = grid.bracket(3.0)  # repeated energy: upper side, unit base
    assert (b.i, b.t, b.unit_base) == (2, 0.5, True)
    assert grid.bracket(9.0).t == 1.0


def test_unit_base_keeps_normalization_where_cartesian_does_not():
    p1 = Panel((0.0, 1.0, 2.0), (0.0, 1.0, 0.0))  # triangle on [0, 2]
    p2 = Panel((0.0, 2.0, 4.0), (0.0, 0.5, 0.0))  # triangle on [0, 4]
    for unit in (True, False):
        total = _integrate(
            lambda e, u=unit: interpolate_panels(p1, p1.range, p2, p2.range, 0.5, u, e),
            0,
            4,
        )
        if unit:
            assert total == pytest.approx(1.0, rel=1e-6)
    # unit base gives one triangle on [0, 3] peaked at 1.5; Cartesian two bumps
    peak = interpolate_panels(p1, p1.range, p2, p2.range, 0.5, True, 1.5)
    assert peak == pytest.approx(2 / 3)


def test_unit_base_equals_cartesian_for_equal_ranges():
    p1 = Panel((0.0, 1.0, 3.0), (0.0, 0.8, 0.2))
    p2 = Panel((0.0, 2.0, 3.0), (0.1, 0.5, 0.3))
    for e in (0.3, 1.2, 2.9):
        u = interpolate_panels(p1, p1.range, p2, p2.range, 0.3, True, e)
        c = interpolate_panels(p1, p1.range, p2, p2.range, 0.3, False, e)
        assert u == pytest.approx(c)


def test_parameters_are_interpolated_without_jacobian():
    r1 = Panel((0.0, 2.0), (0.4, 0.4))
    r2 = Panel((0.0, 4.0), (0.4, 0.4))
    assert interpolate_parameter(
        r1, r1.range, r2, r2.range, 0.5, True, 2.9
    ) == pytest.approx(0.4)


def test_panel_integral_and_mean_are_exact():
    p = Panel((1.0, 2.0, 4.0), (0.0, 2.0, 0.0))
    assert p.integral() == pytest.approx(3.0)
    assert p.mean() == pytest.approx((1 + 2 + 4) / 3)  # centroid of a triangle
    h = Panel((0.0, 1.0, 3.0), (0.5, 0.25, 0.0), law=1)
    assert h.integral() == pytest.approx(1.0)


# --------------------------------------------------------------- kinematics
@pytest.mark.parametrize("awr", [0.99917, 11.9, 236.0])
def test_isotropic_cm_scattering_has_lab_mean_cosine_two_thirds_a(awr):
    # for A ~ 1, mu_lab ~ sqrt(1 + mu_cm) near mu_cm = -1: the midpoint rule
    # converges slowly there, hence 2e-5
    mean = _integrate(lambda mu: 0.5 * two_body(1e3, mu, awr)[1], -1, 1, 20000)
    assert mean == pytest.approx(2 / (3 * awr), rel=2e-5)


def test_elastic_energy_range_and_lethargy_gain():
    awr, e = 15.858, 1e6
    assert two_body(e, -1.0, awr)[0] == pytest.approx(elastic_alpha(awr) * e)
    assert two_body(e, 1.0, awr)[0] == pytest.approx(e)
    mean_log = _integrate(lambda mu: 0.5 * math.log(e / two_body(e, mu, awr)[0]), -1, 1)
    assert mean_log == pytest.approx(mean_lethargy_gain(awr), rel=1e-5)
    assert mean_lethargy_gain(1.0) == 1.0


def test_inelastic_threshold():
    awr, q = 236.0, -44.9e3
    threshold = -q * (awr + 1) / awr
    assert cm_energy(threshold, awr, q) == pytest.approx(0.0, abs=1e-6)
    assert cm_energy(1.1 * threshold, awr, q) > 0


@pytest.mark.parametrize("a, r", [(0.0, 0.3), (0.5, 0.0), (2.0, 0.7), (6.0, 1.0)])
def test_kalbach_mann_shape(a, r):
    assert _integrate(lambda mu: kalbach_mann(mu, a, r), -1, 1) == pytest.approx(
        1.0, rel=1e-6
    )
    mean = _integrate(lambda mu: mu * kalbach_mann(mu, a, r), -1, 1)
    assert mean == pytest.approx(kalbach_mann_mean_cosine(a, r), abs=1e-6)


# ------------------------------------------------------------ real data
NUCLIDES = [
    "H1",
    "O16",
    "U234",
    "U235",
    "U238",
    "Zr90",
    "Zr91",
    "Zr92",
    "Zr94",
    "Zr96",
    "B10",
    "B11",
]
ACTINIDES = ("U234", "U235", "U238")


@pytest.fixture(scope="module")
def evals():
    lib = EndfLibrary()
    if not (lib.root / "archives" / "neutrons-version.VIII.1.tar.gz").exists():
        pytest.skip("ENDF/B-VIII.1 neutron data not downloaded")
    return {nuc: Evaluation.from_file(lib.neutron(nuc)) for nuc in NUCLIDES}


def test_every_distribution_of_the_pin_cell_nuclides_is_readable(evals):
    for ev in evals.values():
        for mf, mt in ev.sections:
            if mf == 4:
                angular_distribution(ev, mt)
            elif mf == 5:
                energy_distribution(ev, mt)
            elif mf == 6:
                energy_angle(ev, mt)


@pytest.mark.parametrize("nuclide", ACTINIDES)
def test_total_nu_is_prompt_plus_delayed(evals, nuclide):
    ev = evals[nuclide]
    nt, nd, npr = (nu_bar(ev, mt) for mt in (452, 455, 456))
    assert nt.table is not None
    for e in nt.table.x:  # nu is written with 6 significant digits
        assert nt(e) == pytest.approx(npr(e) + nd(e), rel=1e-5)


@pytest.mark.parametrize("nuclide", ACTINIDES)
def test_delayed_groups(evals, nuclide):
    groups = delayed_groups(evals[nuclide])
    assert len(groups.decay_constants) == 6
    assert list(groups.half_lives) == sorted(groups.half_lives, reverse=True)
    assert math.fsum(groups.abundances(0.0253)) == pytest.approx(1.0, abs=1e-6)
    # delayed neutrons are soft: a few hundred keV, against ~2 MeV for prompt
    assert all(1e5 < m < 1e6 for m in groups.mean_energies(0.0253))


@pytest.mark.parametrize("nuclide", ("U235", "U238"))
def test_multiplicity_distribution(evals, nuclide):
    """sum nu P~ = 1 holds tightly. The nu_p implied by P~ (1/sum P~) differs
    from MF1 MT456 by up to 1.2% (U-235) / 1.9% (U-238) at the tabulated
    energies: an inconsistency in ENDF/B-VIII.1, so only the normalized P(nu)
    should be used, with nu-bar from MF1."""
    ev = evals[nuclide]
    md = multiplicity_distribution(ev)
    assert md is not None
    nu_p = nu_bar(ev, 456)
    for e in md.tables[0].x:
        assert math.fsum(n * p for n, p in enumerate(md.scaled(e))) == pytest.approx(
            1.0, abs=1e-6
        )
        assert md.mean(e) == pytest.approx(nu_p(e), rel=0.02)
        assert math.fsum(md.probabilities(e)) == pytest.approx(1.0)


@pytest.mark.parametrize("nuclide", ACTINIDES)
def test_fission_spectra_are_normalized(evals, nuclide):
    chi = energy_distribution(evals[nuclide], 18)
    for part in chi.partials:  # U-238 at 20 MeV is off by 3.2e-4; the rest < 1e-6
        for panel in part.law.panels:  # type: ignore[union-attr]
            assert panel.integral() == pytest.approx(1.0, abs=5e-4)
    # mean energy ~2 MeV, rising slowly with incident energy
    assert 1.9e6 < chi.mean(0.0253) < 2.1e6
    assert chi.mean(1.4e7) > chi.mean(0.0253)


def test_hydrogen_elastic_is_isotropic_in_cm(evals):
    ev = evals["H1"]
    ad = angular_distribution(ev, 2)
    assert ad.frame == CM
    awr = ev.header.awr
    mean_lab = _integrate(lambda mu: ad(1.0, mu) * two_body(1.0, mu, awr)[1], -1, 1)
    assert mean_lab == pytest.approx(2 / (3 * awr), rel=1e-4)


@pytest.mark.parametrize("nuclide", ["O16", "U238", "Zr90"])
def test_elastic_distributions_are_positive(evals, nuclide):
    ad = angular_distribution(evals[nuclide], 2)
    assert ad.grid is not None
    for e in ad.grid.energies:
        assert min(ad(e, -1 + k / 20) for k in range(41)) > -1e-6


def test_continuum_spectra_are_normalized(evals):
    for ev in evals.values():
        for mf, mt in ev.sections:
            if mf != 6:
                continue
            for p in energy_angle(ev, mt).products:
                if isinstance(p.distribution, Continuum):
                    for panel in p.distribution.panels:
                        assert panel.component(0).integral() == pytest.approx(
                            1.0, abs=1e-4
                        )


def test_continuum_angular_distributions_integrate_to_the_energy_pdf(evals):
    rng = random.Random(2)
    for nuclide, mt in (("U238", 91), ("O16", 91), ("Zr90", 16)):
        (p,) = energy_angle(evals[nuclide], mt).products[:1]
        d = p.distribution
        assert isinstance(d, Continuum)
        for _ in range(5):
            e = rng.uniform(d.grid.energies[1], d.grid.energies[-1])
            e_out = rng.uniform(0, 0.5 * e)
            f0 = d.energy_pdf(e, e_out)
            assert _integrate(
                lambda mu, d=d, e=e, eo=e_out: d(e, eo, mu), -1, 1, 400
            ) == pytest.approx(f0, rel=1e-6, abs=1e-15)


def test_discrete_two_body_distributions(evals):
    ea = energy_angle(evals["U235"], 51)
    (p,) = ea.products
    assert isinstance(p.distribution, TwoBody)
    assert ea.frame == CM
    assert p.distribution.moments(1e6, 2)[0] == 1.0
