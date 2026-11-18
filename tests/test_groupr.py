"""Tests for group averaging (fission.processing.groupr)."""

import math
import random

import numpy as np
import pytest

from fission.data import EndfLibrary
from fission.endf import Evaluation
from fission.endf.distributions.angular import CM, AngularDistribution
from fission.endf.distributions.mf6 import Continuum, energy_angle
from fission.endf.pointwise import PointwiseTable
from fission.processing.groupr import (
    Kernel,
    continuum_kernel,
    group_integrals,
    group_library,
    hat_integrals,
    two_body_kernel,
)
from fission.processing.groups import XMAS_172, GroupStructure, equal_lethargy
from fission.processing.weights import Constant, InverseE, ThermalSlowingDownFission

ISOTROPIC = AngularDistribution(2, CM, None, ())


# ------------------------------------------------------------ structures, weights
def test_xmas_172_structure():
    assert XMAS_172.size == 172
    assert XMAS_172.boundaries[0] == pytest.approx(1.9640e7, rel=1e-4)
    assert XMAS_172.boundaries[-1] == pytest.approx(1e-5, rel=1e-4)
    g = XMAS_172.find(0.0253)
    lo, hi = XMAS_172.bounds(g)
    assert lo <= 0.0253 <= hi
    assert XMAS_172.find(1e8) is None


def test_weight_is_continuous_with_continuous_slope():
    w = ThermalSlowingDownFission.at(293.6)
    for e in (w.eb, w.ec):
        below, above = w(e * (1 - 1e-9)), w(e * (1 + 1e-9))
        assert below == pytest.approx(above, rel=1e-7)
        h = 1e-4
        slope_below = math.log(w(e) / w(e * (1 - h))) / -math.log(1 - h)
        slope_above = math.log(w(e * (1 + h)) / w(e)) / math.log(1 + h)
        assert slope_below == pytest.approx(slope_above, abs=2e-3)
        assert slope_below == pytest.approx(-1.0, abs=2e-3)


# ------------------------------------------------------------ hat integrals
def test_flux_of_inverse_e_is_the_lethargy_width():
    groups = equal_lethargy(5, 1.0, 1e5)
    flux = group_integrals(
        [1.0, 1e5], [[1.0, 1.0]], InverseE(), groups, Kernel.constant(1.0, 1e5), False
    )
    assert flux[:, 0, 0] == pytest.approx([math.log(10.0)] * 5, rel=1e-12)


def test_linear_sigma_and_linear_kernel_are_integrated_exactly():
    groups = GroupStructure("two", (10.0, 4.0, 1.0))
    x, s = [1.0, 3.0, 10.0], [2.0, 5.0, 1.0]
    z, k = (1.0, 6.0, 10.0), np.array([[1.0], [0.5], [2.0]])
    out = group_integrals(x, [s], Constant(), groups, Kernel(z, k), False)[:, 0, 0]

    def brute(a, b, n=200000):
        e = np.linspace(a, b, n + 1)
        f = np.interp(e, x, s) * np.interp(e, z, k[:, 0])
        return np.trapezoid(f, e)

    assert out == pytest.approx([brute(4.0, 10.0), brute(1.0, 4.0)], rel=1e-8)


@pytest.mark.parametrize(
    "weight", [Constant(), InverseE(), ThermalSlowingDownFission(0.1, 1e3)]
)
def test_rust_hat_integrals_match_python(weight):
    from fission.processing.groupr import _hat_integrals

    rng = random.Random(4)
    x = sorted({10 ** rng.uniform(-3, 4) for _ in range(300)})
    cols = [[rng.uniform(0, 5) for _ in x] for _ in range(3)]
    z = sorted({10 ** rng.uniform(-3, 4) for _ in range(40)})
    groups = equal_lethargy(12, 1e-3, 1e4)
    py = hat_integrals(x, cols, weight, groups, z)
    rust = _hat_integrals(x, cols, weight, groups, z, use_rust=True)
    assert py.keys() == rust.keys()
    for key, v in py.items():
        assert rust[key] == pytest.approx(v, rel=1e-12, abs=1e-300)


# ------------------------------------------------------------ transfer kernels
def test_hydrogen_like_isotropic_scattering_is_uniform_in_energy():
    """For A = 1 and isotropic CM scattering, E' is uniform on [0, E]."""
    groups = equal_lethargy(6, 1e-2, 1e4)
    e = 7000.0
    k = two_body_kernel(ISOTROPIC, 1.0, 0.0, groups, 1, e)
    for g in range(groups.size):
        lo, hi = groups.bounds(g)
        assert k[0, g] == pytest.approx(max(0.0, min(hi, e) - lo) / e, abs=1e-12)


@pytest.mark.parametrize("awr", [1.0, 11.9, 236.0])
def test_isotropic_kernel_conserves_and_has_mean_cosine_two_thirds_a(awr):
    groups = equal_lethargy(40, 1.0, 1e7)
    e = 2.5e6
    k = two_body_kernel(ISOTROPIC, awr, 0.0, groups, 1, e)
    # for A = 1 a fraction floor/E is scattered below the structure (E' uniform on [0, E])
    lost = 1.0 / e if awr == 1.0 else 0.0
    assert k[0].sum() == pytest.approx(1.0 - lost, rel=1e-12)
    assert k[1].sum() == pytest.approx(2 / (3 * awr), rel=1e-7)


def test_inelastic_kernel_is_zero_below_threshold():
    groups = equal_lethargy(10, 1.0, 1e7)
    q, awr = -45e3, 236.0
    assert not two_body_kernel(ISOTROPIC, awr, q, groups, 0, 44e3).any()
    assert two_body_kernel(ISOTROPIC, awr, q, groups, 0, 1e6)[0].sum() == pytest.approx(
        1.0
    )


# ------------------------------------------------------------ real data
@pytest.fixture(scope="module")
def lib():
    lib = EndfLibrary()
    if not (lib.root / "archives" / "neutrons-version.VIII.1.tar.gz").exists():
        pytest.skip("ENDF/B-VIII.1 neutron data not downloaded")
    return lib


COARSE = equal_lethargy(10, 1e-5, 2e7, "10 groups")


def _small_table(ev, mts, lo=1e-5, hi=2e7, n=4000):
    """A coarse pointwise table straight from MF3 (no resonances needed here)."""
    e = [lo * (hi / lo) ** (k / (n - 1)) for k in range(n)]
    return PointwiseTable(e, {mt: [ev.cross_section(mt)(x) for x in e] for mt in mts})


def test_continuum_kernels_conserve_neutrons(lib):
    for nuclide, mt in (("U238", 91), ("O16", 91), ("Zr90", 16)):
        ev = Evaluation.from_file(lib.neutron(nuclide))
        ea = energy_angle(ev, mt)
        cont = ea.products[0].distribution
        assert isinstance(cont, Continuum)
        e = cont.grid.energies[len(cont.grid.energies) // 2]
        k = continuum_kernel(cont, ea.frame, ev.header.awr, COARSE, 1, e)
        assert k[0].sum() == pytest.approx(
            cont.panels[0].component(0).integral(), rel=2e-3
        )
        assert k[0].sum() == pytest.approx(1.0, rel=1e-4)


@pytest.fixture(scope="module")
def u235(lib):
    ev = Evaluation.from_file(lib.neutron("U235"))
    from fission.endf.cache import cached_pointwise

    table = cached_pointwise("U235")
    return group_library(ev, table, COARSE, ThermalSlowingDownFission.at(293.6), lmax=1)


def test_library_conservation_and_fission(u235):
    gl = u235
    awr = 233.0248
    alpha = ((awr - 1) / (awr + 1)) ** 2
    for g in range(COARSE.size):
        if COARSE.bounds(g)[0] * alpha < COARSE.boundaries[-1]:
            continue  # scattering can leave the bottom of the structure
        assert gl.scatter_elastic[0][g].sum() == pytest.approx(gl.sigma[2][g], rel=1e-6)
    # no upscattering for a target at rest
    assert np.allclose(np.tril(gl.scatter_elastic[0], -1), 0.0)
    assert gl.chi.sum() == pytest.approx(1.0)
    assert gl.chi_delayed.sum(axis=1) == pytest.approx([1.0] * 6)
    ratio = gl.nu_fission / gl.sigma[18]
    assert np.all((ratio > 2.3) & (ratio < 5.5))  # nu-bar
    assert gl.nu_prompt[-1] == pytest.approx(2.414, rel=2e-3)  # thermal nu_p


def test_inelastic_matrix_matches_total_inelastic_cross_section(u235):
    """Below the (n,2n) threshold the level-by-level and continuum matrices
    add up to the evaluated total inelastic cross section MT4."""
    gl = u235
    for g in range(COARSE.size):
        lo, hi = COARSE.bounds(g)
        if 1e5 < lo and hi < 5e6:
            assert gl.scatter_inelastic[0][g].sum() == pytest.approx(
                gl.sigma[4][g], rel=2e-3
            )


def test_one_over_e_group_average_is_exact(lib):
    """With w = 1/E, sigma_g = Int sigma dE/E / ln(E_hi/E_lo) exactly."""
    from fission.endf.cache import cached_pointwise

    table = cached_pointwise("U238")
    groups = GroupStructure("res", (2e4, 100.0, 20.0, 5.0, 0.5))
    out = group_integrals(
        table.energies, [table.xs[102]], InverseE(), groups, Kernel.constant(0.5, 2e4)
    )[:, 0, 0]
    for g in range(groups.size):
        lo, hi = groups.bounds(g)
        assert out[g] == pytest.approx(table.integral_over_e(102, lo, hi), rel=1e-9)
