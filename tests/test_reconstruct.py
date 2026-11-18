"""Tests for adaptive pointwise reconstruction of resonance cross sections."""

import math
import random

import pytest

from fission.data import EndfLibrary
from fission.endf import Evaluation
from fission.endf.channel import wave_number
from fission.endf.reconstruct import (
    PointwiseXS,
    adaptive_grid,
    reconstruct,
    resolved_range_xs,
)
from fission.endf.reich_moore import reich_moore
from fission.endf.resonances import LGroup, ResolvedRange, Resonance


def lorentzian(e):
    peak = 1000.0 / (1 + ((e - 5.0) / 0.01) ** 2)
    return (10.0 + peak, 0.1 * peak, 0.0)


def test_adaptive_grid_meets_tolerance_between_points():
    tol = 1e-3
    energies, values = adaptive_grid(lorentzian, [1.0, 5.0, 10.0], tol=tol, atol=1e-6)
    p = PointwiseXS(energies, *map(list, zip(*values, strict=True)))
    rng = random.Random(0)
    for _ in range(2000):
        e = rng.uniform(1.0, 10.0)
        exact = lorentzian(e)[0]
        assert (
            abs(p("elastic", e) - exact) <= 3 * tol * exact
        )  # midpoint test is not a bound
    # points crowd where sigma curves: much denser within 10 widths of the peak
    # than far away (Lorentzian tails fall like 1/x^2 and still need points)
    density_near = sum(1 for e in energies if abs(e - 5.0) < 0.1) / 0.2
    density_far = sum(1 for e in energies if 1.0 <= e <= 4.0) / 3.0
    assert density_near > 20 * density_far


def test_linear_function_needs_no_refinement():
    seeds = [1.0, 2.0, 7.0]
    energies, _ = adaptive_grid(lambda e: (2 * e + 1, 3.0, 0.0), seeds)
    assert energies == seeds


def test_integral_over_e_is_exact_for_lin_lin():
    p = PointwiseXS([1.0, 4.0], [0.0, 0.0], [3.0, 9.0], [0.0, 0.0])  # sigma = 1 + 2E
    assert p.integral_over_e("capture", 1.0, 4.0) == pytest.approx(math.log(4) + 2 * 3)
    assert p.integral_over_e("capture", 2.0, 3.0) == pytest.approx(math.log(1.5) + 2)


def test_single_level_capture_area():
    """Reconstructed capture of one isolated level integrates to (2 pi^2/k^2) g Gn Gg / G."""
    res = Resonance(20.0, 0.5, 2e-3, 23e-3, 0.0, 0.0)
    group = LGroup(l=0, awri=236.0, apl=0.95, resonances=(res,))
    rr = ResolvedRange(1e-5, 1e4, 3, 1, 0.0, 0.95, 0, 1, (group,))
    gamma = res.gn + res.gg
    lo, hi = 20.0 - 400 * gamma, 20.0 + 400 * gamma

    def f(e):
        x = reich_moore(rr, e)
        return x.elastic, x.capture, x.fission

    energies, values = adaptive_grid(
        f, [lo, 20.0 - gamma / 2, 20.0, 20.0 + gamma / 2, hi]
    )
    p = PointwiseXS(energies, *map(list, zip(*values, strict=True)))
    area = sum(
        0.5 * (p.capture[i] + p.capture[i + 1]) * (energies[i + 1] - energies[i])
        for i in range(len(energies) - 1)
    )
    k = wave_number(20.0, 236.0)
    expected = 2 * math.pi**2 / k**2 * res.gn * res.gg / gamma
    # the window holds 2/pi atan(800) of a Lorentzian; 1/k^2 and Gn(E) skew slightly
    assert area == pytest.approx(expected * 2 / math.pi * math.atan(800), rel=5e-3)


# --------------------------------------------------- real ENDF/B-VIII.1 data
@pytest.fixture(scope="module")
def u238():
    lib = EndfLibrary()
    if not (lib.root / "archives" / "neutrons-version.VIII.1.tar.gz").exists():
        pytest.skip("ENDF/B-VIII.1 neutron data not downloaded")
    return Evaluation.from_file(lib.neutron("U238"))


@pytest.fixture(scope="module")
def u238_window(u238):
    return reconstruct(u238, 6.3, 7.0)


def test_grid_points_are_exact(u238, u238_window):
    _, f = resolved_range_xs(u238)
    for i in range(0, len(u238_window), 25):
        e = u238_window.energies[i]
        assert (
            u238_window.elastic[i],
            u238_window.capture[i],
            u238_window.fission[i],
        ) == f(e)


def test_peak_is_on_the_grid(u238_window):
    assert 6.67428 in u238_window.energies
    i = u238_window.energies.index(6.67428)
    assert u238_window.capture[i] == max(u238_window.capture)


def test_interpolation_accuracy_in_window(u238, u238_window):
    _, f = resolved_range_xs(u238)
    rng = random.Random(1)
    for _ in range(200):
        e = rng.uniform(6.3, 7.0)
        exact = f(e)
        for j, name in enumerate(("elastic", "capture")):
            assert abs(u238_window(name, e) - exact[j]) <= 3e-3 * abs(exact[j]) + 1e-4
