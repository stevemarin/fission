"""Tests for Multi-Level Breit-Wigner cross sections."""

import pytest

from fission.data import EndfLibrary
from fission.endf import Evaluation
from fission.endf.mlbw import mlbw
from fission.endf.reconstruct import resolved_range_xs
from fission.endf.resonances import LGroup, ResolvedRange, Resonance
from fission.endf.slbw import slbw


def one_level(res, spi=0.0, l=0):
    group = LGroup(l=l, awri=89.13, apl=0.0, resonances=(res,))
    return ResolvedRange(1e-5, 2e5, 2, 0, spi, 0.716, 0, 0, (group,))


@pytest.mark.parametrize(
    "res, spi, l",
    [
        (Resonance(3861.0, 0.5, 10.9, 0.15, 0.0, 0.0), 0.0, 0),
        (Resonance(292.0, 2.0, 0.6, 0.14, 0.0, 0.0), 2.5, 0),  # Zr-91-like, I = 5/2
        (Resonance(4008.0, 1.5, 0.34, 0.09, 0.0, 0.0), 0.0, 1),  # p-wave: shift factor
    ],
)
@pytest.mark.parametrize("energy", [0.0253, 100.0, 290.0, 3850.0, 3861.0, 4008.0, 1e5])
def test_single_level_equals_slbw(res, spi, l, energy):
    rr = one_level(res, spi, l)
    m, s = mlbw(rr, energy), slbw(rr, energy)
    assert m.capture == pytest.approx(s.capture, rel=1e-12)
    assert m.elastic == pytest.approx(s.elastic, rel=1e-9)


def test_two_levels_interfere_in_elastic_only():
    a = Resonance(100.0, 0.5, 2.0, 0.1, 0.0, 0.0)
    b = Resonance(104.0, 0.5, 2.0, 0.1, 0.0, 0.0)
    group = LGroup(l=0, awri=89.13, apl=0.0, resonances=(a, b))
    rr = ResolvedRange(1e-5, 2e5, 2, 0, 0.0, 0.716, 0, 0, (group,))
    m, s = mlbw(rr, 102.0), slbw(rr, 102.0)
    assert m.capture == pytest.approx(s.capture, rel=1e-12)
    assert m.elastic != pytest.approx(s.elastic, rel=1e-2)
    assert min(mlbw(rr, 99.0 + 0.05 * i).elastic for i in range(120)) >= 0


# --------------------------------------------------- real ENDF/B-VIII.1 data
@pytest.fixture(scope="module")
def lib():
    lib = EndfLibrary()
    if not (lib.root / "archives" / "neutrons-version.VIII.1.tar.gz").exists():
        pytest.skip("ENDF/B-VIII.1 neutron data not downloaded")
    return lib


ABUNDANCE = {
    "Zr90": 0.5145,
    "Zr91": 0.1122,
    "Zr92": 0.1715,
    "Zr94": 0.1738,
    "Zr96": 0.0280,
}


def test_natural_zirconium_thermal_capture(lib):
    """Mughabghab, Atlas of Neutron Resonances: natural Zr sigma_gamma = 0.185 b."""
    capture = 0.0
    for nuclide, abundance in ABUNDANCE.items():
        _, f = resolved_range_xs(Evaluation.from_file(lib.neutron(nuclide)))
        capture += abundance * f(0.0253)[1]
    assert capture == pytest.approx(0.185, rel=3e-2)


def test_u234_thermal_capture(lib):
    """Mughabghab: U-234 sigma_gamma = 100.1 +- 1.3 b."""
    _, f = resolved_range_xs(Evaluation.from_file(lib.neutron("U234")))
    assert f(0.0253)[1] == pytest.approx(100.1, rel=2e-2)


def test_zr91_elastic_positive(lib):
    rr = Evaluation.from_file(lib.neutron("Zr91")).resonances[0].resolved
    energies = [1e-5 * (rr.eh / 1e-5) ** (i / 400) for i in range(400)]
    assert min(mlbw(rr, e).elastic for e in energies) > 0
