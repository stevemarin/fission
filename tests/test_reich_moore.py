"""Tests for Reich-Moore resonance cross sections."""

import math

import pytest

from fission.data import EndfLibrary
from fission.endf import Evaluation
from fission.endf.reich_moore import (
    reich_moore,
    solve_spin_group,
    spin_groups,
    total_from_optical_theorem,
)
from fission.endf.resonances import LGroup, ResolvedRange, Resonance
from fission.endf.slbw import slbw


def one_level(resonance, spi=0.0, l=0):
    """A Reich-Moore range holding a single resonance."""
    group = LGroup(l=l, awri=236.0, apl=0.95, resonances=(resonance,))
    return ResolvedRange(
        el=1e-5, eh=1e4, lrf=3, naps=1, spi=spi, ap=0.95, lad=0, nlsc=1, groups=(group,)
    )


@pytest.mark.parametrize(
    "resonance, spi",
    [
        (Resonance(6.67, 0.5, 1.5e-3, 23e-3, 0.0, 0.0), 0.0),  # capture only
        (Resonance(1.14, 4.0, 1.5e-5, 35e-3, 0.12, -0.05), 3.5),  # two fission channels
        (Resonance(-2.0, 3.0, 4e-4, 36e-3, 0.11, 0.09), 3.5),  # bound level
    ],
)
@pytest.mark.parametrize("energy", [0.0253, 1.0, 1.14, 1.3, 6.6, 6.67, 50.0])
def test_single_level_equals_slbw(resonance, spi, energy):
    """With one resonance there is nothing to interfere with: RM == SLBW exactly."""
    rr = one_level(resonance, spi)
    rm, bw = reich_moore(rr, energy), slbw(rr, energy)
    assert rm.capture == pytest.approx(bw.capture, rel=1e-9)
    assert rm.fission == pytest.approx(bw.fission, rel=1e-9, abs=1e-300)
    assert rm.elastic == pytest.approx(bw.elastic, rel=1e-9)


def test_two_levels_interfere():
    """Two levels in the same spin group: RM differs from the SLBW sum."""
    a = Resonance(10.0, 0.5, 5e-3, 20e-3, 0.0, 0.0)
    b = Resonance(10.1, 0.5, 5e-3, 20e-3, 0.0, 0.0)
    group = LGroup(l=0, awri=236.0, apl=0.95, resonances=(a, b))
    rr = ResolvedRange(1e-5, 1e4, 3, 1, 0.0, 0.95, 0, 1, (group,))
    assert reich_moore(rr, 10.05).capture != pytest.approx(
        slbw(rr, 10.05).capture, rel=1e-2
    )


# --------------------------------------------------- real ENDF/B-VIII.1 data
def _resolved(nuclide):
    lib = EndfLibrary()
    if not (lib.root / "archives" / "neutrons-version.VIII.1.tar.gz").exists():
        pytest.skip("ENDF/B-VIII.1 neutron data not downloaded")
    (iso,) = Evaluation.from_file(lib.neutron(nuclide)).resonances
    return iso.resolved


@pytest.fixture(scope="module")
def u238():
    return _resolved("U238")


@pytest.fixture(scope="module")
def u235():
    return _resolved("U235")


THERMAL = 0.0253


def test_u238_thermal_capture(u238):
    # Mughabghab, Atlas of Neutron Resonances: 2.680 +- 0.019 b
    assert reich_moore(u238, THERMAL).capture == pytest.approx(2.68, rel=1e-2)


def test_u235_thermal_values(u235):
    # Mughabghab: fission 585.1 b, capture ~99 b
    xs = reich_moore(u235, THERMAL)
    assert xs.fission == pytest.approx(585.1, rel=5e-3)
    assert xs.capture == pytest.approx(99.0, rel=2e-2)


def test_slbw_is_badly_wrong_for_u235(u235):
    """U-235 levels overlap (large fission widths): ignoring interference costs ~10%."""
    rm, bw = reich_moore(u235, THERMAL), slbw(u235, THERMAL)
    assert abs(bw.fission / rm.fission - 1) > 0.05


def log_energies(lo, hi, n):
    return [lo * (hi / lo) ** (i / (n - 1)) for i in range(n)]


@pytest.mark.parametrize("name", ["u238", "u235"])
def test_unitarity_and_positivity(name, request):
    """|U_nn| <= 1 in every spin group; hence elastic >= 0. Capture >= 0 too."""
    rr = request.getfixturevalue(name)
    groups = spin_groups(rr)
    for energy in log_energies(1e-4, 0.999 * rr.eh, 150):
        for sg in groups:
            assert abs(solve_spin_group(sg, energy).u_nn) <= 1 + 1e-12
        xs = reich_moore(rr, energy)
        assert xs.elastic >= 0 and xs.capture >= -1e-9 * xs.total and xs.fission >= 0


def test_elastic_stays_positive_where_slbw_fails(u238):
    energies = [6.50 + 0.005 * i for i in range(30)]  # across the interference dip
    assert min(slbw(u238, e).elastic for e in energies) < 0
    assert min(reich_moore(u238, e).elastic for e in energies) > 0


def test_isolated_capture_matches_slbw(u238):
    """At the 6.67 eV peak there is no capture interference: RM ~ SLBW."""
    assert reich_moore(u238, 6.67428).capture == pytest.approx(
        slbw(u238, 6.67428).capture, rel=1e-4
    )


def test_potential_scattering_weight_bookkeeping(u235):
    """Every l up to NLSC-1 contributes (2l+1) worth of spin weight. At low
    energy only s-waves matter, so elastic far from resonances should be near
    4 pi AP^2 (U-235 AP = 9.6 fm -> ~11.6 b) modulated by resonance tails."""
    xs = reich_moore(u235, 1e-5)
    assert 5 < xs.elastic < 30
    assert math.isfinite(xs.total)


@pytest.mark.parametrize("name", ["u238", "u235"])
@pytest.mark.parametrize("energy", [0.0253, 1.14, 6.67428, 20.87, 102.5, 1500.0])
def test_optical_theorem(name, energy, request):
    """Total from 1 - Re U equals elastic + capture + fission (unitarity)."""
    rr = request.getfixturevalue(name)
    if energy > rr.eh:
        pytest.skip("outside the resolved range")
    xs = reich_moore(rr, energy)
    assert total_from_optical_theorem(rr, energy) == pytest.approx(xs.total, rel=1e-8)
