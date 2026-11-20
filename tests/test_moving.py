"""Free-gas elastic kernel with the 0 K cross section, and the march with it."""

import math

import numpy as np
import pytest

from fission.pincell.moving import erfcx, free_gas_kernel
from fission.processing.thermal import free_gas_total

A = 236.0058  # U-238
ALPHA = ((A - 1) / (A + 1)) ** 2
CONSTANT = (np.array([1e-5, 2e7]), np.array([10.0, 10.0]))  # sigma_0 = 10 b


def test_erfcx_matches_the_library_erfc():
    y = np.linspace(0.0, 25.9, 5001)
    ref = np.exp(y * y) * np.array([math.erfc(v) for v in y])
    assert np.max(np.abs(erfcx(y) / ref - 1)) < 2e-13


def test_constant_cross_section_gives_the_free_gas_total():
    """Integrated over E_out, the kernel is the free-gas sigma_s(E) of a
    constant sigma_0, (1 + 1/2y^2) erf(y) + ... (chapter 8). The mesh
    integration error goes as dx^2."""
    errors = []
    for dx in (1 / 8, 1 / 16):
        k = free_gas_kernel("X", *CONSTANT, A, 600.0, 9.0, 11.0, dx=dx, dz=1 / 128)
        e = k.energies
        rows = np.nonzero(~np.isnan(k.norm))[0]
        closed = np.array([10.0 * free_gas_total(A, 600.0, e[m]) for m in rows])
        errors.append(np.max(np.abs(k.norm[rows] / closed - 1)))
    assert errors[1] < 5e-4
    assert errors[0] / errors[1] == pytest.approx(4.0, rel=0.1)


def test_cold_gas_is_the_target_at_rest_box():
    """T -> 0: uniform 1 / ((1 - alpha) E) on [alpha E, E], zero outside."""
    k = free_gas_kernel("X", *CONSTANT, A, 0.25, 9.9, 10.0, dx=0.5, dz=1 / 32)
    e = k.energies
    m = k.top
    row = k.dense()[m] * (1 - ALPHA) * e[m]
    inside = (e > ALPHA * e[m] + 0.02) & (e < e[m] - 0.02)
    outside = (e < ALPHA * e[m] - 0.02) | (e > e[m] + 0.02)
    assert np.max(np.abs(row[inside] - 1)) < 5e-3  # edge normalization, ~sqrt(T)
    assert np.max(row[outside]) < 1e-12


@pytest.fixture(scope="module")
def u238():
    from fission.endf import cache
    from fission.processing.cache import COLD

    t = cache.cached_pointwise("U238", **COLD)
    return np.asarray(t.energies), np.asarray(t.xs[2])


def test_resonance_kernel_integrates_to_the_broadened_cross_section(u238):
    """Across U-238's 6.67 eV resonance the kernel's integral is the
    Doppler-broadened sigma_s of doppler.py (an independent calculation of
    the same free-gas average) to the dx^2 mesh error; and it satisfies
    detailed balance, E e^{-E/kT} sigma(E -> E') symmetric, to rounding."""
    from fission.endf.pointwise import PointwiseTable
    from fission.processing.doppler import BOLTZMANN_EV_PER_K, broadened_at

    e0, s0 = u238
    k = free_gas_kernel("U238", e0, s0, A, 600.0, 6.0, 7.5, dx=1 / 16, dz=1 / 128)
    e = k.energies
    rows = np.nonzero(~np.isnan(k.norm))[0]
    sel = rows[:: max(len(rows) // 8, 1)]
    lo, hi = np.searchsorted(e0, [1.0, 20.0])
    sub = PointwiseTable(list(e0[lo:hi]), {2: list(s0[lo:hi])})
    ref = np.array(broadened_at(sub, A, 600.0, list(e[sel]))[2])
    assert np.max(np.abs(k.norm[sel] / ref - 1)) < 1e-3
    kt = BOLTZMANN_EV_PER_K * 600.0
    full = k.dense() * np.nan_to_num(k.norm)[:, None]
    worst = 0.0
    for m in rows[40:-40:9]:
        for n in rows[40:-40:9]:
            if full[m, n] > 1e-6 * full[m].max() and full[n, m] > 0:
                left = e[m] * math.exp(-(e[m] - e[n]) / kt) * full[m, n]
                worst = max(worst, abs(left / (e[n] * full[n, m]) - 1))
    assert worst < 1e-12


def _resonance_problem(temperature, u238):
    """U-238 (fuel density) with a constant moderator, 4-40 eV, homogeneous."""
    from fission.pincell.shielding import march
    from fission.processing.cache import cached_doppler

    e0, s0 = u238
    k = free_gas_kernel("U238", e0, s0, A, temperature, 4.0, 12.0)
    tb = cached_doppler("U238", temperature)
    eb = np.asarray(tb.energies)
    grid = np.union1d(eb[(eb >= 4.0) & (eb <= 40.0)], np.geomspace(4.0, 40.0, 4000))
    mesh = k.energies
    grid = np.union1d(grid, mesh[(mesh >= 4.0) & (mesh <= 40.0)])[::-1]
    n_u, moderator = 0.0223, 1.1
    total = n_u * np.interp(grid, eb, np.asarray(tb.xs[1])) + moderator
    scatter = n_u * np.interp(grid, eb, np.asarray(tb.xs[2]))
    start = int(np.searchsorted(-grid, -(40.0 / 2.2), side="right"))
    history = (moderator + scatter[:start]) / total[:start]
    u = np.log(grid[0] / grid)
    window = (grid > 5.5) & (grid < 8.0)

    def absorption(moving):
        psi = march(grid, total, [np.full_like(grid, moderator), scatter], [0.5, ALPHA],
                    np.zeros_like(grid), history, start, moving=moving)  # fmt: skip
        rate = (total - moderator - scatter) * psi
        return np.trapezoid(rate[window], u[window]), psi

    return absorption, k, grid


def test_moving_targets_raise_resonance_capture_more_when_hot(u238):
    """Thermal motion lets neutrons scatter into the 6.67 eV resonance from
    below as well as from above: the march captures more there than with
    targets at rest, and more at 1200 K than at 600 K. Far from the
    resonance the flux is the same."""
    effect = []
    for temperature in (600.0, 1200.0):
        absorption, k, grid = _resonance_problem(temperature, u238)
        rest, psi_rest = absorption(None)
        hot, psi_hot = absorption([None, k])
        effect.append(hot / rest - 1)
        far = (grid > 13.0) & (grid < 16.0)  # above the kernel's range and tail
        assert np.max(np.abs(psi_hot[far] / psi_rest[far] - 1)) < 1e-6
    assert 0 < effect[0] < effect[1] < 0.01


def test_pure_scattering_keeps_the_slowing_down_flux():
    """Constant sigma_0, no absorption: the moving march keeps psi = 1 to the
    small physical heating correction (U-238 carries ~1% of the slowing down
    here and gains ~2kT/E per collision: ~5e-4 at 1200 K), including at the
    bottom of the range, where neutrons from below scatter up into it. Near
    the top of the kernel's range the flux has a physical seam: below it
    targets move, above they are at rest (OpenMC's resonance scattering has
    the same one at its energy_max)."""
    from fission.pincell.shielding import march

    k = free_gas_kernel("X", *CONSTANT, A, 1200.0, 4.0, 12.0)
    grid = np.geomspace(4.0, 40.0, 6000)
    mesh = k.energies
    grid = np.union1d(grid, mesh[(mesh >= 4.0) & (mesh <= 40.0)])[::-1]
    u238, moderator = np.full_like(grid, 0.223), np.full_like(grid, 0.5)
    start = int(np.searchsorted(-grid, -(40.0 / 2.2), side="right"))
    psi = march(grid, u238 + moderator, [moderator, u238], [0.5, ALPHA],
                np.zeros_like(grid), np.ones(start), start, moving=[None, k])  # fmt: skip
    away = grid < 11.0
    assert np.max(np.abs(psi[away] - 1.0)) < 1e-3
    # below the seam's echo through the moderator (alpha = 0.5: 11.8 -> 5.9 eV)
    assert np.all(psi[grid < 5.5] > 1.0)  # heating, not a dip, at the bottom
    assert np.max(np.abs(psi[(grid > 11.7) & (grid < 12.1)] - 1)) > 0.05  # the seam
