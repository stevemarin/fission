"""Tests for reading MF2 resonance parameters (no cross sections yet)."""

import pytest
from helpers import endf_line

from fission.data import EndfLibrary
from fission.endf import Evaluation, ResolvedRange, UnresolvedRange
from fission.endf.records import Cursor
from fission.endf.resonances import EmptyRange, read_mf2


def mf2(fields_per_line):
    return Cursor([endf_line(f, mat=9237, mf=2, mt=151) for f in fields_per_line])


def test_reich_moore_section():
    cursor = mf2(
        [
            ["9.223800+4", "2.360058+2", "0", "0", "1", "0"],  # HEAD: NIS = 1
            ["9.223800+4", "1.0", "0", "0", "1", "0"],  # isotope: NER = 1
            ["1.0-5", "2.0+4", "1", "3", "0", "0"],  # range: resolved, Reich-Moore
            ["0.0", "0.948", "0", "0", "2", "2"],  # SPI, AP, LAD, 0, NLS=2, NLSC
            ["2.360058+2", "0.948", "0", "0", "12", "2"],  # l = 0, two resonances
            ["6.674+0", "0.5", "1.492-3", "2.27-2", "0.0", "0.0"],
            ["2.087+1", "0.5", "1.008-2", "2.28-2", "0.0", "0.0"],
            ["2.360058+2", "0.948", "1", "0", "6", "1"],  # l = 1, one resonance
            ["4.4075+0", "1.5", "5.56-8", "2.25-2", "0.0", "0.0"],
        ]
    )
    (iso,) = read_mf2(cursor)
    rr = iso.resolved
    assert isinstance(rr, ResolvedRange) and rr.formalism == "Reich-Moore"
    assert (rr.el, rr.eh, rr.spi, rr.nlsc) == (1e-5, 2e4, 0.0, 2)
    assert [g.l for g in rr.groups] == [0, 1]
    first = rr.groups[0].resonances[0]
    assert (first.energy, first.j, first.gn, first.gg) == (
        6.674,
        0.5,
        1.492e-3,
        2.27e-2,
    )
    assert len(rr.resonances) == 3


def test_mlbw_column_order():
    """BW formalisms list GT (total width) first; it must not be read as GN."""
    cursor = mf2(
        [
            ["4.009000+4", "8.913000+1", "0", "0", "1", "0"],
            ["4.009000+4", "1.0", "0", "0", "1", "0"],
            ["1.0-5", "2.0+5", "1", "2", "0", "0"],  # MLBW
            ["0.0", "0.716", "0", "0", "1", "0"],
            ["8.913000+1", "0.0", "0", "0", "6", "1"],
            ["3.861+3", "0.5", "21.7", "10.9", "0.8", "0.0"],  # ER, AJ, GT, GN, GG, GF
        ]
    )
    res = read_mf2(cursor)[0].resolved.resonances[0]
    assert (res.gn, res.gg, res.gfa) == (10.9, 0.8, 0.0)


def test_no_resonances():
    cursor = mf2(
        [
            ["1.001000+3", "9.991673-1", "0", "0", "1", "0"],
            ["1.001000+3", "1.0", "0", "0", "1", "0"],
            ["1.0-5", "1.0+5", "0", "0", "0", "0"],  # LRU = 0
            ["0.5", "1.275246", "0", "0", "0", "0"],
        ]
    )
    (iso,) = read_mf2(cursor)
    assert isinstance(iso.ranges[0], EmptyRange) and iso.resolved is None


# --------------------------------------------------- real ENDF/B-VIII.1 data
def _isotope(nuclide):
    lib = EndfLibrary()
    if not (lib.root / "archives" / "neutrons-version.VIII.1.tar.gz").exists():
        pytest.skip("ENDF/B-VIII.1 neutron data not downloaded")
    (iso,) = Evaluation.from_file(lib.neutron(nuclide)).resonances
    return iso


@pytest.fixture(scope="module")
def u238():
    return _isotope("U238")


@pytest.fixture(scope="module")
def u235():
    return _isotope("U235")


def test_u238_ranges(u238):
    rr, ur = u238.resolved, u238.unresolved
    assert rr.formalism == "Reich-Moore" and (rr.el, rr.eh) == (1e-5, 2e4)
    assert isinstance(ur, UnresolvedRange) and (ur.el, ur.eh) == (2e4, 1.5e5)
    assert {g.l for g in rr.groups} == {0, 1}


def test_u238_6_67_eV_resonance(u238):
    """Mughabghab, Atlas of Neutron Resonances: E = 6.674 eV, J = 1/2,
    Gamma_n = 1.493 meV, Gamma_gamma = 23.0 meV."""
    res = min(u238.resolved.resonances, key=lambda r: abs(r.energy - 6.674))
    assert res.energy == pytest.approx(6.674, abs=2e-3)
    assert res.j == 0.5
    assert res.gn == pytest.approx(1.493e-3, rel=2e-2)
    assert res.gg == pytest.approx(23.0e-3, rel=3e-2)


def test_u238_lowest_s_wave_resonances(u238):
    s_wave = sorted(
        r.energy for r in u238.resolved.groups[0].resonances if r.energy > 0
    )
    assert s_wave[:3] == pytest.approx([6.67, 20.87, 36.68], abs=0.01)


@pytest.mark.parametrize("name", ["u235", "u238"])
def test_spins_obey_angular_momentum_coupling(name, request):
    """J must lie between |I - l - 1/2| and I + l + 1/2 (target spin I, neutron spin 1/2)."""
    rr = request.getfixturevalue(name).resolved
    for group in rr.groups:
        lowest = abs(abs(rr.spi - group.l) - 0.5)
        allowed = {lowest + k for k in range(int(rr.spi + group.l + 0.5 - lowest) + 1)}
        assert {abs(r.j) for r in group.resonances} <= allowed


def test_u235_spin_and_fission_channels(u235):
    rr = u235.resolved
    assert rr.formalism == "Reich-Moore" and rr.spi == 3.5 and rr.eh == 2250
    assert any(r.gfa != 0 and r.gfb != 0 for r in rr.resonances)  # two fission channels
    assert any(r.gfa < 0 or r.gfb < 0 for r in rr.resonances)  # signed amplitudes
