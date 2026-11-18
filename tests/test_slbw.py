"""Tests for channel quantities and Single-Level Breit-Wigner cross sections."""

import math

import pytest

from fission.data import EndfLibrary
from fission.endf import Evaluation
from fission.endf.channel import (
    K_CONSTANT,
    penetrability_shift_phase,
    spin_factor,
    wave_number,
)
from fission.endf.slbw import potential_scattering, single_resonance, slbw


# ----------------------------------------------------- channel quantities
def riccati_bessel(l, rho):
    """Regular/irregular Riccati-Bessel functions F_l, G_l for l = 0, 1, 2."""
    s, c = math.sin(rho), math.cos(rho)
    if l == 0:
        return s, c
    if l == 1:
        return s / rho - c, c / rho + s
    if l == 2:
        return (3 / rho**2 - 1) * s - 3 * c / rho, (3 / rho**2 - 1) * c + 3 * s / rho
    a, b = 15 / rho**3 - 6 / rho, 15 / rho**2 - 1
    return a * s - b * c, a * c + b * s


@pytest.mark.parametrize(
    "l, rho",
    [(l, r) for l in (0, 1, 2) for r in (0.01, 0.3, 1.0, 2.5)]
    + [
        (3, r) for r in (0.3, 1.0, 2.5, 5.0)
    ],  # l=3 Bessel forms cancel badly at tiny rho
)
def test_hard_sphere_functions_match_bessel_definitions(l, rho):
    """P = rho/(F^2+G^2), S = rho (F F' + G G')/(F^2+G^2), phi = atan(F/G).

    An independent route to the closed forms used by ENDF.
    """
    f, g = riccati_bessel(l, rho)
    h = 1e-6
    fp = (riccati_bessel(l, rho + h)[0] - riccati_bessel(l, rho - h)[0]) / (2 * h)
    gp = (riccati_bessel(l, rho + h)[1] - riccati_bessel(l, rho - h)[1]) / (2 * h)
    p, s, phi = penetrability_shift_phase(l, rho, rho)
    assert p == pytest.approx(rho / (f * f + g * g), rel=1e-10)
    assert s == pytest.approx(
        rho * (f * fp + g * gp) / (f * f + g * g), rel=1e-6, abs=1e-9
    )
    # phases may differ by a multiple of pi; compare via tan
    assert math.tan(phi) == pytest.approx(f / g, rel=1e-9)


def test_wave_number_constant():
    # ENDF-102 quotes 2.196771e-3 (older constants); CODATA 2018 gives 2.196808e-3
    assert pytest.approx(2.19680e-3, rel=2e-5) == K_CONSTANT
    assert wave_number(4.0, 1e12) == pytest.approx(
        2 * K_CONSTANT
    )  # heavy target: k ~ sqrt(E)


def test_spin_factor():
    assert spin_factor(0.5, 0.0) == 1.0  # even-even target, s-wave: always J = 1/2
    assert spin_factor(4.0, 3.5) == pytest.approx(9 / 16)  # U-235, J = 4
    assert spin_factor(3.0, 3.5) == pytest.approx(7 / 16)  # J = 3; the two add to 1


# ---------------------------------------------------- U-238, 6.67 eV level
@pytest.fixture(scope="module")
def u238():
    lib = EndfLibrary()
    if not (lib.root / "archives" / "neutrons-version.VIII.1.tar.gz").exists():
        pytest.skip("ENDF/B-VIII.1 neutron data not downloaded")
    (iso,) = Evaluation.from_file(lib.neutron("U238")).resonances
    rr = iso.resolved
    group = rr.groups[0]
    res = min(group.resonances, key=lambda r: abs(r.energy - 6.674))
    return rr, group, res


def test_peak_capture(u238):
    """At E = Er (s-wave, no shift): sigma_g = (4 pi / k^2) g Gn Gg / G^2."""
    rr, group, res = u238
    k = wave_number(res.energy, group.awri)
    gamma = res.gn + res.gg + abs(res.gfa) + abs(res.gfb)
    expected = 4 * math.pi / k**2 * res.gn * res.gg / gamma**2
    assert slbw(rr, res.energy, only=[res]).capture == pytest.approx(
        expected, rel=1e-12
    )
    assert expected == pytest.approx(22_765, rel=1e-3)


def test_subthreshold_fission_is_tiny(u238):
    """The evaluation gives U-238's low resonances a ~10 neV fission width
    (subthreshold fission): real, but ~1e-6 of capture."""
    rr, _, res = u238
    xs = slbw(rr, res.energy, only=[res])
    assert 0 < xs.fission < 1e-5 * xs.capture


def test_capture_area(u238):
    """Integrated capture over the resonance: (2 pi^2 / k^2) g Gn Gg / G, the
    shape-independent quantity measured in experiments. Integrate +-W and
    compare with the fraction of a Lorentzian inside the window."""
    rr, group, res = u238
    gamma = res.gn + res.gg + abs(res.gfa) + abs(res.gfb)
    w, n = 20 * gamma, 4000
    energies = [res.energy - w + 2 * w * i / n for i in range(n + 1)]
    values = [slbw(rr, e, only=[res]).capture for e in energies]
    area = (2 * w / n) * (sum(values) - 0.5 * (values[0] + values[-1]))
    k = wave_number(res.energy, group.awri)
    full = 2 * math.pi**2 / k**2 * res.gn * res.gg / gamma
    assert area == pytest.approx(
        full * 2 / math.pi * math.atan(2 * w / gamma), rel=5e-3
    )


def test_capture_is_1_over_v_far_below(u238):
    rr, _, res = u238
    ref = slbw(rr, 1e-5, only=[res]).capture * math.sqrt(1e-5)
    for energy in (1e-4, 1e-3, 1e-2):
        value = slbw(rr, energy, only=[res]).capture * math.sqrt(energy)
        assert value == pytest.approx(ref, rel=5e-3)


def test_interference_dip_below_and_bump_above(u238):
    rr, group, res = u238
    pot = potential_scattering(rr, group, res.energy)
    assert pot == pytest.approx(
        4 * math.pi * group.apl**2, rel=1e-3
    )  # 4 pi a^2 = 11.3 b
    below = min(
        slbw(rr, res.energy - 0.01 * i, only=[res]).elastic for i in range(1, 40)
    )
    above = slbw(rr, res.energy + 0.1, only=[res]).elastic
    assert below < 0.1 * pot < pot < above


def test_p_wave_potential_scattering_negligible_at_low_energy(u238):
    rr, _, _ = u238
    assert potential_scattering(rr, rr.groups[1], 6.674) < 1e-6


def test_single_resonance_ingredients(u238):
    rr, group, res = u238
    c = single_resonance(rr, group, res, res.energy)
    assert c.gn == pytest.approx(res.gn)  # Gn(E) = Gn(Er) at E = Er
    assert c.shifted_energy == res.energy  # s-wave: no shift
    assert c.g == 1.0


def test_slbw_elastic_can_go_negative(u238):
    """Known SLBW defect: summing independent levels' interference terms can
    drive the elastic cross section below zero. Multi-level formalisms
    (Reich-Moore, step 2c) do not have this problem."""
    rr, _, _ = u238
    assert slbw(rr, 6.55).elastic < 0
