"""Tests for Doppler broadening (free-gas kernel on linear-linear tables)."""

import math
import random

import numpy as np
import pytest

from fission import _rust
from fission.data import EndfLibrary
from fission.endf import Evaluation
from fission.endf.pointwise import PointwiseTable, pointwise
from fission.processing.doppler import (
    GAUSS_RULES,
    alpha,
    broaden,
    broaden_table,
    broadened_at,
    broadening_weights,
    cullen_weisbin,
    doppler_width,
    thin,
)


def _logspace(lo, hi, n):
    return [lo * (hi / lo) ** (i / (n - 1)) for i in range(n)]


def test_alpha_and_doppler_width():
    # U-238 (AWR 236.0058) at 293.6 K: k_B T = 0.025300 eV
    assert alpha(236.0058, 293.6) == pytest.approx(236.0058 / 0.0253006, rel=1e-5)
    # Gamma_D = 2 sqrt(E kT / A) = 4 sqrt(E / alpha) / 2
    a = alpha(236.0058, 900.0)
    assert doppler_width(6.674, 236.0058, 900.0) == pytest.approx(
        2 * math.sqrt(6.674 / a)
    )


def test_constant_cross_section_matches_closed_form():
    """A constant sigma_0 broadens to sigma_0 [(1 + 1/2y^2) erf(y) + exp(-y^2)/(y sqrt(pi))].

    At low energy the moving targets add relative speed, so the reaction rate
    per neutron stays finite and sigma rises like 1/v (the free-gas
    'thermal rise' of scattering cross sections).
    """
    e = _logspace(1e-12, 1e3, 40)  # coarse is fine: the table is exact
    table = PointwiseTable(e, {2: [20.0] * len(e)})
    awr, temp = 0.9992, 293.6
    a = alpha(awr, temp)
    for energy in [1e-4, 1e-3, 0.0253, 0.1, 1.0, 10.0]:
        y = math.sqrt(a * energy)
        exact = (1 + 0.5 / y**2) * math.erf(y) + math.exp(-(y**2)) / (
            y * math.sqrt(math.pi)
        )
        (got,) = broadened_at(table, awr, temp, [energy])[2]
        assert got == pytest.approx(20.0 * exact, rel=1e-12)


def test_one_over_v_is_unchanged():
    """A 1/v cross section is temperature independent: v sigma is constant,
    so the averaged reaction rate does not depend on the target motion."""
    e = _logspace(
        1e-5, 1e2, 7 * 2000 + 1
    )  # dense: lin-lin represents 1/sqrt(E) to 1e-7
    table = PointwiseTable(e, {102: [0.3326 * math.sqrt(0.0253 / x) for x in e]})
    probe = [1e-5, 3e-4, 0.0253, 0.7, 12.0]
    for temp in [293.6, 1200.0]:
        got = broadened_at(table, 0.9992, temp, probe)[102]
        for energy, value in zip(probe, got, strict=True):
            assert value == pytest.approx(0.3326 * math.sqrt(0.0253 / energy), rel=1e-6)


def test_quadrature_matches_closed_form_kernel():
    """Gauss-Legendre panels agree with Cullen & Weisbin's analytic integration."""
    rng = random.Random(7)
    e = sorted(rng.uniform(1e-3, 3.0) for _ in range(30))
    sigma = [rng.uniform(0.0, 50.0) for _ in e]
    awr, temp = 10.0, 600.0
    a = alpha(awr, temp)
    x = [math.sqrt(a * v) for v in e]
    for energy in [1e-3, 0.05, 0.5, 1.7, 2.9]:
        y = math.sqrt(a * energy)
        w = broadening_weights(x, y)
        quad = math.fsum(wi * sigma[i] for i, wi in w.items())
        assert quad == pytest.approx(cullen_weisbin(x, sigma, y), rel=1e-11)


def test_weights_are_shared_by_all_reactions_and_keep_the_sum():
    rng = random.Random(1)
    e = sorted(rng.uniform(1.0, 2.0) for _ in range(200))
    parts = {mt: [rng.uniform(0, 10) for _ in e] for mt in (2, 18, 102)}
    total = [math.fsum(v) for v in zip(*parts.values(), strict=True)]
    hot = broaden(PointwiseTable(e, {1: total, **parts}), 236.0, 900.0)
    for j in range(len(e)):
        assert hot.xs[1][j] == pytest.approx(
            math.fsum(hot.xs[mt][j] for mt in parts), rel=1e-13
        )


def _lorentzian_table(e0=10.0, gamma=0.05, peak=1e4, lo=6.0, hi=14.0, n=801):
    e = [lo + (hi - lo) * i / (n - 1) for i in range(n)]
    s = [peak / (1 + (2 * (x - e0) / gamma) ** 2) + 5.0 for x in e]
    return PointwiseTable(e, {102: s})


def test_successive_broadening_adds_temperatures():
    """Maxwellians compose: broadening to T1, then by T2 - T1, equals broadening to T2.

    The intermediate (300 K) table is linearly interpolated on a 4 meV grid,
    which limits the agreement to ~1e-4 here; a wrong temperature dependence
    of the kernel would show up at the percent level.
    """
    table, awr = _lorentzian_table(), 50.0
    probe = [9.9, 9.97, 10.0, 10.02, 10.1]
    direct = broadened_at(table, awr, 900.0, probe)[102]
    grid = [8.5 + 0.004 * i for i in range(751)]
    warm = PointwiseTable(grid, broadened_at(table, awr, 300.0, grid))
    twice = broadened_at(warm, awr, 600.0, probe)[102]
    for d, t in zip(direct, twice, strict=True):
        assert t == pytest.approx(d, rel=3e-4)


def test_emax_keeps_cold_values_above():
    table = _lorentzian_table(lo=9.0, hi=11.0, n=201)
    hot = broaden(table, 50.0, 900.0, emax=10.0)
    above = [j for j, e in enumerate(table.energies) if e > 10.0]
    assert all(hot.xs[102][j] == table.xs[102][j] for j in above)
    assert hot.xs[102][above[0] - 1] != table.xs[102][above[0] - 1]
    assert broaden(table, 50.0, 0.0).xs == table.xs


def test_gauss_legendre_constants():
    for n, (nodes, weights) in GAUSS_RULES.items():
        ref_nodes, ref_weights = np.polynomial.legendre.leggauss(n)
        assert np.allclose(nodes, ref_nodes, rtol=0, atol=1e-15)
        assert np.allclose(weights, ref_weights, rtol=0, atol=1e-15)


def test_thin_drops_only_redundant_points():
    e = [float(i) for i in range(11)]
    line = [(2.0 * x + 1.0,) for x in e]
    assert thin(e, line, 1e-3, 0.0) == [0, 10]
    kinked = [(abs(x - 4.0),) for x in e]  # a corner at x = 4 must stay
    assert thin(e, kinked, 1e-3, 0.0) == [0, 4, 10]


def test_broadened_table_is_accurate_between_its_points():
    """The thinned table reproduces the broadened function within ~tol."""
    table, awr, temp, tol = _lorentzian_table(n=401), 50.0, 300.0, 1e-3
    hot = broaden_table(table, awr, temp, tol=tol)
    assert len(hot) < len(table)  # smoother than the 0 K curve
    rng = random.Random(3)
    probe = [rng.uniform(8.0, 12.0) for _ in range(200)]
    exact = broadened_at(table, awr, temp, probe)[102]
    for e, ref in zip(probe, exact, strict=True):
        assert hot(102, e) == pytest.approx(ref, rel=2 * tol, abs=1e-5)


def test_sparse_kinks_are_resolved():
    """Where 0 K points are far apart compared with the Doppler width, the
    broadened curve is straight between them except for rounded corners at
    each 0 K kink. The corners on both ends of an interval can cancel at its
    midpoint, so the midpoint test alone would accept a wrong chord; extra
    seeds at a few Doppler widths around each 0 K point prevent that."""
    e = [90.0, 100.0, 110.0, 120.0]  # a ramp: one concave and one convex corner
    table = PointwiseTable(e, {18: [50.0, 50.0, 10.0, 10.0]})
    awr, temp, tol = 233.0, 293.6, 1e-3  # Gamma_D(100 eV) = 0.21 eV
    for use_rust in (False, True):
        hot = broaden_table(table, awr, temp, tol=tol, use_rust=use_rust)
        probe = [95.0 + 0.01 * i for i in range(2001)]  # inside the table
        exact = broadened_at(table, awr, temp, probe)[18]
        for energy, ref in zip(probe, exact, strict=True):
            assert hot(18, energy) == pytest.approx(ref, rel=2 * tol, abs=1e-5), energy


# ------------------------------------------------------------ Rust parity
def _random_table(seed, n=120):
    rng = random.Random(seed)
    e = sorted({rng.uniform(1e-3, 5.0) for _ in range(n)})
    cols = [[rng.uniform(0.0, 30.0) for _ in e] for _ in range(3)]
    return e, cols


@pytest.mark.parametrize("seed", [0, 1, 2])
def test_rust_weights_match_python(seed):
    e, _ = _random_table(seed)
    a = alpha(20.0, 600.0)
    x = [math.sqrt(a * v) for v in e]
    for y in [0.05, 1.0, x[0], x[len(x) // 2], x[-1] + 2.0, x[-1] + 10.0]:
        py = broadening_weights(x, y)
        start, w = _rust.doppler_weights(x, y)
        dense = np.zeros(len(x))
        for i, wi in py.items():
            dense[i] = wi
        rust = np.zeros(len(x))
        rust[start : start + len(w)] = w
        np.testing.assert_allclose(rust, dense, rtol=1e-13, atol=1e-15 * max(dense))


def test_rust_evaluate_matches_python():
    e, cols = _random_table(4)
    table = PointwiseTable(e, {2: cols[0], 18: cols[1], 102: cols[2]})
    at = [1e-4, 0.01, 0.3, 2.0, 4.99, 7.0]
    py = broadened_at(table, 20.0, 900.0, at)
    rust = _rust.doppler_evaluate(e, cols, 20.0, 900.0, at)
    for mt, col in zip((2, 18, 102), rust, strict=True):
        np.testing.assert_allclose(col, py[mt], rtol=1e-13)


def test_rust_thin_matches_python():
    rng = random.Random(5)
    e = [i * 0.01 for i in range(1, 2000)]
    cols = [[math.sin(v) + 1e-4 * rng.random() for v in e], [math.exp(-v) for v in e]]
    rows = list(zip(*cols, strict=True))
    assert _rust.doppler_thin(e, cols, 1e-3, 1e-5) == thin(e, rows, 1e-3, 1e-5)


def test_rust_broaden_table_matches_python():
    """Same grid point for point, same values."""
    cold = _lorentzian_table(lo=8.0, hi=12.0, n=201)
    cold.xs[2] = [10.0 + 0.1 * v for v in cold.energies]
    py = broaden_table(cold, 50.0, 300.0, emax=11.5, use_rust=False)
    rust = broaden_table(cold, 50.0, 300.0, emax=11.5, use_rust=True)
    assert rust.energies == py.energies
    for mt in cold.xs:
        np.testing.assert_allclose(rust.xs[mt], py.xs[mt], rtol=1e-12)


# --------------------------------------------------- real ENDF/B-VIII.1 data
@pytest.fixture(scope="module")
def u238():
    lib = EndfLibrary()
    if not (lib.root / "archives" / "neutrons-version.VIII.1.tar.gz").exists():
        pytest.skip("ENDF/B-VIII.1 neutron data not downloaded")
    ev = Evaluation.from_file(lib.neutron("U238"))
    (level,) = [
        r
        for g in ev.resonances[0].resolved.groups
        for r in g.resonances
        if 6 < r.energy < 7
    ]
    return ev.header.awr, level, pointwise(ev, mts=[102], lo=4.0, hi=9.5)


def test_u238_667_resonance_at_900k(u238):
    """Peak and area of the 6.67 eV capture resonance against the
    single-level psi function (Bethe-Placzek) estimate."""
    awr, level, cold = u238
    grid = [e for e in cold.energies if 6.0 <= e <= 7.4]
    hot = broadened_at(cold, awr, 900.0, grid)[102]

    def area(v):
        return math.fsum(
            0.5 * (v[i] + v[i + 1]) * (grid[i + 1] - grid[i])
            for i in range(len(grid) - 1)
        )

    cold_vals = [cold(102, e) for e in grid]
    assert area(hot) == pytest.approx(area(cold_vals), rel=1e-3)  # area is conserved

    # psi(0, xi) = sqrt(pi)/2 xi exp(xi^2/4) erfc(xi/2), xi = Gamma / Gamma_D
    gamma = (
        level.gn + level.gg + abs(level.gfa) + abs(level.gfb)
    )  # total width, 0.0242 eV
    xi = gamma / doppler_width(level.energy, awr, 900.0)
    psi = 0.5 * math.sqrt(math.pi) * xi * math.exp(xi**2 / 4) * math.erfc(xi / 2)
    assert max(hot) == pytest.approx(max(cold_vals) * psi, rel=5e-3)  # agrees to 0.2%


@pytest.mark.slow
def test_u238_resonance_integral_vs_temperature(u238):
    """Full cached tables: the dilute capture RI does not depend on T (area is
    conserved), but with a flux dip (narrow-resonance weighting with
    sigma_0 = 50 b) capture grows with T: the Doppler effect."""
    import numpy as np

    from fission.endf.cache import cached_pointwise
    from fission.processing.cache import cached_doppler

    def shielded(t, s0=50.0):
        e, g, tot = (np.asarray(v) for v in (t.energies, t.xs[102], t.xs[1]))
        k = (e >= 0.5) & (e <= 2e4)
        return np.trapezoid(g[k] * s0 / ((tot[k] + s0) * e[k]), e[k])

    cold = cached_pointwise("U238")
    hot = {temp: cached_doppler("U238", temp) for temp in (293.6, 900.0)}
    ri0 = cold.integral_over_e(102, 0.5, 2e4)
    for t in hot.values():
        assert t.integral_over_e(102, 0.5, 2e4) == pytest.approx(ri0, rel=5e-4)
    assert shielded(cold) < shielded(hot[293.6]) < shielded(hot[900.0])
    assert shielded(hot[900.0]) / shielded(hot[293.6]) == pytest.approx(1.09, abs=0.01)
