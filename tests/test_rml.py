"""R-Matrix Limited (LRF=7) resonance parameters and cross sections."""

import pytest

from fission.data import EndfLibrary
from fission.endf import Evaluation
from fission.endf.reconstruct import resolved_range_xs
from fission.endf.reich_moore import reich_moore
from fission.endf.resonances import (
    Channel,
    LGroup,
    ParticlePair,
    ResolvedRange,
    Resonance,
    RMLRange,
    RMLResonance,
    RMLSpinGroup,
)
from fission.endf.rml import RMatrixLimited

AWRI, AP = 27.737, 0.45
GAMMA = ParticlePair(0.0, AWRI + 1, 0.0, 14.0, 1.0, 0.0, 0.0, 0, 0, 102)
NEUTRON = ParticlePair(1.0, AWRI, 0.0, 14.0, 0.5, 0.0, 0.0, 1, 0, 2)

# (l, J, [(E_r, G_n, G_g)]): two interfering s-wave levels, a bound level,
# and p-wave levels in both J; l = 1, J = 1.5 has no levels (potential only)
LEVELS = [
    (0, 0.5, [(-3000.0, 900.0, 1.2), (55000.0, 4000.0, 0.8), (61000.0, 2500.0, 1.5)]),
    (1, 0.5, [(20000.0, 30.0, 0.6)]),
    (1, 1.5, []),
]


def lrf3() -> ResolvedRange:
    groups = []
    for l in (0, 1):
        res = tuple(
            Resonance(e, j, gn, gg, 0.0, 0.0)
            for ll, j, levels in LEVELS
            if ll == l
            for e, gn, gg in levels
        )
        groups.append(LGroup(l=l, awri=AWRI, apl=AP, resonances=res))
    return ResolvedRange(1e-5, 1e6, 3, 1, 0.0, AP, 0, 2, tuple(groups))


def lrf7() -> RMLRange:
    groups = []
    for l, j, levels in LEVELS:
        channels = (Channel(0, 0, 0.0, 0.0, 0.0, 0.0), Channel(1, l, 0.5, 0.0, AP, AP))
        res = tuple(RMLResonance(e, (gg, gn)) for e, gn, gg in levels)
        groups.append(RMLSpinGroup(j, 1 if l % 2 == 0 else -1, channels, res))
    return RMLRange(1e-5, 1e6, 0, 3, 0, (GAMMA, NEUTRON), tuple(groups))


@pytest.mark.parametrize("energy", [0.0253, 1.0, 2e4, 5.4e4, 5.5e4, 5.8e4, 6.1e4, 3e5])
def test_lrf7_reproduces_lrf3(energy):
    """The same levels written as LRF=3 and as LRF=7 (KRM=3) are the same
    Reich-Moore cross sections; the potential scattering of the empty
    (l = 1, J = 1.5) group comes from NLSC in LRF=3, from a listed channel
    in LRF=7."""
    a, b = reich_moore(lrf3(), energy), RMatrixLimited(lrf7())(energy)
    assert b.elastic == pytest.approx(a.elastic, rel=1e-12)
    assert b.capture == pytest.approx(a.capture, rel=1e-12)
    assert b.fission == 0.0


@pytest.mark.parametrize("energy", [0.0253, 2e4, 5.5e4, 6.1e4])
def test_unitarity(energy):
    """elastic + capture = total from the optical theorem."""
    model = RMatrixLimited(lrf7())
    xs = model(energy)
    assert xs.elastic + xs.capture == pytest.approx(model.total(energy), rel=1e-9)


# --------------------------------------------------- real ENDF/B-VIII.1 data
def _evaluation(nuclide):
    lib = EndfLibrary()
    if not (lib.root / "archives" / "neutrons-version.VIII.1.tar.gz").exists():
        pytest.skip("ENDF/B-VIII.1 neutron data not downloaded")
    return Evaluation.from_file(lib.neutron(nuclide))


def test_si28_structure():
    rr = _evaluation("Si28").resonances[0].resolved
    assert isinstance(rr, RMLRange)
    assert (rr.krm, rr.ifg, rr.eh) == (3, 0, 1.75e6)
    assert [p.mt for p in rr.pairs] == [102, 2]
    assert [(g.j, g.parity, g.channels[1].l) for g in rr.groups] == [
        (0.5, 1, 0), (0.5, -1, 1), (1.5, -1, 1), (1.5, 1, 2),
        (2.5, 1, 2), (2.5, -1, 3), (3.5, -1, 3),
    ]  # fmt: skip


@pytest.mark.parametrize(
    "nuclide, elastic, capture",
    # NJOY 2016.78 RECONR on the same files (0.1% tolerance), 0.0253 eV
    [("Si28", 1.8945, 0.186175), ("Si29", 2.72193, 0.119011),
     ("Si30", 2.44922, 0.107474), ("Fe54", 2.1216, 2.25323)],
)  # fmt: skip
def test_thermal_values_match_njoy(nuclide, elastic, capture):
    """Resonance part + MF3 background at 0.0253 eV, against NJOY."""
    _, f = resolved_range_xs(_evaluation(nuclide))
    el, ca, _ = f(0.0253)
    assert el == pytest.approx(elastic, rel=5e-6)
    assert ca == pytest.approx(capture, rel=5e-6)


def test_inelastic_pairs_are_refused():
    """Fe-57 has an inelastic pair (14.4 keV level): not implemented."""
    rr = _evaluation("Fe57").resonances[0].resolved
    with pytest.raises(NotImplementedError, match="MT=51"):
        RMatrixLimited(rr)
