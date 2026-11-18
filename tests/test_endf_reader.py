"""Tests for the ENDF-6 reader.

Unit tests use small hand-written records. The physics tests read the real
ENDF/B-VIII.1 files and compare with independently published values; they are
skipped when the data has not been downloaded (python -m nodal.data ...).
"""

import math

import pytest
from helpers import endf_line

from nodal.data import EndfLibrary
from nodal.endf import Evaluation, Tabulated1D, mt_name
from nodal.endf.interpolation import HISTOGRAM, LIN_LIN, LIN_LOG, LOG_LIN, LOG_LOG
from nodal.endf.records import Cursor, parse_float


# ------------------------------------------------------------------ records
@pytest.mark.parametrize(
    "text, value",
    [
        (" 9.223800+4", 92238.0),
        ("-1.234567-5", -1.234567e-5),
        (" 2.360058+2", 236.0058),
        (" 1.0E+00   ", 1.0),
        ("1.0e-5", 1e-5),
        ("           ", 0.0),
        (" 0.5       ", 0.5),
        ("-2.5+10", -2.5e10),
    ],
)
def test_parse_float(text, value):
    assert parse_float(text) == pytest.approx(value, rel=1e-15)


def test_tab1_record_spanning_lines():
    lines = [
        endf_line(["0.0", "2.5+6", "0", "0", "2", "7"]),  # CONT: NR=2, NP=7
        endf_line(["3", "2", "7", "5"]),  # NBT/INT pairs
        endf_line(["1.0-5", "1.0", "2.0-5", "2.0", "3.0-5", "3.0"]),
        endf_line(["4.0-5", "4.0", "5.0-5", "5.0", "6.0-5", "6.0"]),
        endf_line(["7.0-5", "7.0"]),
    ]
    cont, tab = Cursor(lines).tab1()
    assert cont.c2 == 2.5e6 and cont.n2 == 7
    assert tab.breaks == [3, 7] and tab.laws == [2, 5]
    assert tab.x == pytest.approx([1e-5 * i for i in range(1, 8)])
    assert tab.y == [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0]


def test_list_record():
    lines = [
        endf_line(["1.0", "2.0", "0", "0", "8", "0"]),
        endf_line([f"{i}.0" for i in range(1, 7)]),
        endf_line(["7.0", "8.0"]),
    ]
    cont, values = Cursor(lines).list_record()
    assert cont.n1 == 8 and values == [float(i) for i in range(1, 9)]


# ------------------------------------------------------------ interpolation
@pytest.mark.parametrize(
    "law, f",
    [
        (LIN_LIN, lambda x: 3.0 + 2.0 * x),
        (LIN_LOG, lambda x: 3.0 + 2.0 * math.log(x)),
        (LOG_LIN, lambda x: 3.0 * math.exp(-0.7 * x)),
        (LOG_LOG, lambda x: 3.0 * x**-0.5),  # 1/v: sigma ~ E^-1/2
    ],
)
def test_each_law_is_exact_for_its_function_family(law, f):
    t = Tabulated1D([1.0, 4.0], [f(1.0), f(4.0)], [2], [law])
    for x in (1.3, 2.0, 3.7):
        assert t(x) == pytest.approx(f(x), rel=1e-13)


def test_histogram_and_regions():
    t = Tabulated1D(
        [1, 2, 3, 4], [10, 20, 30, 40], breaks=[2, 4], laws=[HISTOGRAM, LIN_LIN]
    )
    assert t(1.5) == 10  # region 1: constant
    assert t(2.5) == pytest.approx(25)  # region 2: linear
    assert t.law(0) == HISTOGRAM and t.law(1) == LIN_LIN and t.law(2) == LIN_LIN


def test_discontinuity_and_outside_range():
    t = Tabulated1D([1, 2, 2, 3], [1, 1, 5, 5])
    assert t(1.9) == pytest.approx(1) and t(2) == 5 and t(2.5) == 5
    assert t(0.5) == 0 and t(3.5) == 0


def test_bad_tables_rejected():
    with pytest.raises(ValueError):
        Tabulated1D([2, 1], [0, 0])
    with pytest.raises(ValueError):
        Tabulated1D([1, 2, 3], [0, 0, 0], breaks=[2], laws=[2])


def test_mt_names():
    assert mt_name(102) == "(n,gamma) capture"
    assert mt_name(51) == "(n,n') level 1"
    assert mt_name(649) == "(n,p) continuum"
    assert mt_name(801) == "(n,alpha) level 1"


def test_evaluation_from_lines():
    head = endf_line(["1.001000+3", "9.991673-1", "0", "0", "0", "0"], mf=3, mt=102)
    cont = endf_line(["2.224648+6", "2.224648+6", "0", "0", "1", "2"], mf=3, mt=102)
    interp = endf_line(["2", "5"], mf=3, mt=102)
    data = endf_line(["1.0-5", "4.0", "1.0+0", "1.264911-2"], mf=3, mt=102)  # 1/v
    send = endf_line([], mf=3, mt=0)
    ev = Evaluation.from_lines([head, cont, interp, data, send])
    xs = ev.cross_section(102)
    assert ev.mat == 125 and ev.reactions == [102]
    assert xs.qi == pytest.approx(2.224648e6)
    assert xs(1e-3) == pytest.approx(0.4, rel=1e-6)  # log-log is exact for 1/v


# --------------------------------------------------- real ENDF/B-VIII.1 data
def _evaluation(nuclide):
    lib = EndfLibrary()
    archive = lib.root / "archives" / "neutrons-version.VIII.1.tar.gz"
    if not archive.exists():
        pytest.skip("ENDF/B-VIII.1 neutron data not downloaded")
    return Evaluation.from_file(lib.neutron(nuclide))


@pytest.fixture(scope="module")
def h1():
    return _evaluation("H1")


@pytest.fixture(scope="module")
def o16():
    return _evaluation("O16")


@pytest.fixture(scope="module")
def u238():
    return _evaluation("U238")


THERMAL = 0.0253


def test_header_identifies_material(h1, u238):
    assert (h1.header.za, h1.header.lfi, h1.header.nsub, h1.header.nver) == (
        1001,
        0,
        10,
        8,
    )
    # AWR is the mass of the *atom* (proton + electron) over the neutron mass
    assert h1.header.awr == pytest.approx(1.007825 / 1.008665, rel=1e-5)
    assert u238.header.za == 92238 and u238.header.lfi == 1
    assert u238.header.zsymam == "92-U -238"


@pytest.mark.parametrize("name", ["h1", "u238"])
def test_directory_matches_file(name, request):
    """MF1/MT451 lists every section and its line count; check it against the file."""
    ev = request.getfixturevalue(name)
    listed = {(mf, mt): nc for mf, mt, nc, _ in ev.header.directory}
    assert listed == {key: len(lines) for key, lines in ev.sections.items()}


def test_hydrogen_thermal_values(h1):
    # Mughabghab, Atlas of Neutron Resonances: sigma_gamma(H) = 0.3326 b;
    # free-atom scattering = bound 82.02 b * (A/(A+1))^2 = 20.5 b
    assert h1.cross_section(102)(THERMAL) == pytest.approx(0.3326, rel=2e-3)
    assert h1.cross_section(2)(THERMAL) == pytest.approx(20.5, rel=1e-2)


def test_hydrogen_capture_follows_1_over_v(h1):
    capture = h1.cross_section(102)
    reference = capture(THERMAL) * math.sqrt(THERMAL)
    for energy in (1e-5, 1e-3, 0.1, 1.0):
        assert capture(energy) * math.sqrt(energy) == pytest.approx(reference, rel=1e-3)


def test_hydrogen_total_is_sum_of_partials(h1):
    total, elastic, capture = (h1.cross_section(mt) for mt in (1, 2, 102))
    for energy in total.table.x:
        assert total(energy) == pytest.approx(
            elastic(energy) + capture(energy), rel=1e-6
        )


def test_oxygen_thermal_scattering(o16):
    # free-atom value from bound 4.232 b * (A/(A+1))^2 = 3.75 b
    assert o16.cross_section(2)(THERMAL) == pytest.approx(3.75, rel=2e-2)


def test_uranium_238_resonance_flag_and_background(u238):
    """LRP=1: below the resonance range MF3 is only a background; the 2.68 b
    thermal capture comes entirely from the MF2 resonance parameters."""
    assert u238.header.lrp == 1
    assert u238.cross_section(102)(THERMAL) < 1e-3
    assert 190e6 < u238.cross_section(18).qi < 210e6  # ~200 MeV per fission
    assert u238.cross_section(1)(1e6) == pytest.approx(
        7.0, rel=3e-2
    )  # above resonances
