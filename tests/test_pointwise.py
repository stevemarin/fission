"""Tests for full-energy-range pointwise tables."""

import math
import random

import pytest

from nodal.data import EndfLibrary
from nodal.endf import Evaluation, fast
from nodal.endf.pointwise import (
    SHADE,
    evaluated_total_inconsistency,
    partial_mts,
    pointwise,
    shaded_points,
)
from nodal.endf.reconstruct import resolved_range_xs


def test_shading_splits_repeated_energies():
    points, jumps = shaded_points([1.0, 2.0, 2.0, 3.0])
    assert jumps == {2.0}
    assert points == {1.0, 2.0 * (1 - SHADE), 2.0 * (1 + SHADE), 3.0}


# --------------------------------------------------- real ENDF/B-VIII.1 data
@pytest.fixture(scope="module")
def lib():
    lib = EndfLibrary()
    if not (lib.root / "archives" / "neutrons-version.VIII.1.tar.gz").exists():
        pytest.skip("ENDF/B-VIII.1 neutron data not downloaded")
    return lib


def _ev(lib, nuclide):
    return Evaluation.from_file(lib.neutron(nuclide))


def test_partials_drop_redundant_sums(lib):
    mts = partial_mts(_ev(lib, "U238"))
    assert {2, 4, 16, 18, 102} <= set(mts)
    assert not {1, 3} & set(mts)
    assert not set(range(51, 92)) & set(mts)  # covered by MT4
    assert not {19, 20, 21, 38} & set(mts)  # covered by MT18


@pytest.mark.parametrize("nuclide", ["H1", "O16", "U235", "U238", "Zr90"])
def test_partials_reproduce_evaluated_total(lib, nuclide):
    """At the evaluation's own breakpoints (outside the resolved range, where
    MF3 is complete), MT1 = sum of our partials to ENDF's 7-digit precision.
    A double-counted or missing reaction would show up here."""
    ev = _ev(lib, nuclide)
    rr = ev.resonances[0].resolved if ev.header.lrp == 1 else None
    grid = sorted(
        {x for mt in [1, *partial_mts(ev)] for x in ev.cross_section(mt).table.x}
    )
    if rr is not None:
        grid = [e for e in grid if e > rr.eh * 1.0001]
    assert evaluated_total_inconsistency(ev, grid) < 1e-4


@pytest.fixture(scope="module")
def h1(lib):
    ev = _ev(lib, "H1")
    return ev, pointwise(ev)


def test_linearization_of_log_log_capture(h1):
    """H-1 capture is stored log-log (1/v); the table must be lin-lin accurate
    to the bisection criterion tol*sigma + atol (atol = 1e-5 b dominates where
    capture falls to a few mb)."""
    ev, table = h1
    capture = ev.cross_section(102)
    assert len(table) > len(capture.table)  # points were added
    rng = random.Random(3)
    for _ in range(500):
        e = 10 ** rng.uniform(-5, math.log10(2e7))
        assert table(102, e) == pytest.approx(capture(e), rel=3e-3, abs=2e-5)


def test_total_is_sum_of_partials(h1):
    _, table = h1
    for i in range(len(table)):
        partials = math.fsum(table.xs[mt][i] for mt in table.mts if mt != 1)
        assert table.xs[1][i] == partials


def test_jumps_are_shaded_not_on_grid(lib):
    ev = _ev(lib, "Zr90")
    table = pointwise(ev, hi=3e5)
    rr = ev.resonances[0].resolved
    assert rr.eh not in table.energies
    assert (
        rr.eh * (1 - SHADE) in table.energies and rr.eh * (1 + SHADE) in table.energies
    )


def test_resolved_region_matches_reconstruction(lib):
    """Inside the resolved range the table is resonance + background."""
    ev = _ev(lib, "Zr90")
    table = pointwise(ev, lo=3000.0, hi=5000.0)
    _, f = resolved_range_xs(ev)
    for i in range(0, len(table), 20):
        el, ca, _ = f(table.energies[i])
        assert table.xs[2][i] == pytest.approx(el, rel=1e-12)
        assert table.xs[102][i] == pytest.approx(ca, rel=1e-12)


def test_u238_window_agrees_with_rust_reconstruction(lib):
    ev = _ev(lib, "U238")
    table = pointwise(ev, lo=6.3, hi=7.0)
    rec = fast.reconstruct(ev, 6.3, 7.0)
    for e in rec.energies[::10]:
        assert table(102, e) == pytest.approx(rec("capture", e), rel=2e-3)


def test_boron_10_thermal_absorption(lib):
    """Soluble boron works through B-10(n,alpha), MT107: 3842 +- 8 b at
    0.0253 eV (Mughabghab). It must survive the choice of partials."""
    table = pointwise(_ev(lib, "B10"))
    assert table(107, 0.0253) == pytest.approx(3842, rel=5e-3)
    assert table(107, 0.0253) / table(1, 0.0253) > 0.99  # almost pure absorber


@pytest.mark.slow
def test_u238_full_table(lib):
    table = pointwise(_ev(lib, "U238"))
    assert table(102, 0.0253) == pytest.approx(2.68, rel=1e-2)
    # capture RI 0.5 eV - 20 MeV: Mughabghab 275 +- 3 b (resolved part ~274.4 b)
    assert table.integral_over_e(102, 0.5, 2e7) == pytest.approx(275.0, rel=1.5e-2)
