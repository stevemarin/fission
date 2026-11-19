"""Parity of the Rust kernels of steps A3-A5 with their Python references."""

import numpy as np
import pytest

from fission import _rust
from fission.data import EndfLibrary
from fission.endf import Evaluation
from fission.endf.distributions.angular import (
    CM,
    AngularDistribution,
    LegendrePanel,
    TabulatedPanel,
    angular_distribution,
)
from fission.endf.distributions.incident import IncidentGrid
from fission.endf.distributions.mf6 import energy_angle
from fission.endf.interpolation import Tabulated1D
from fission.pincell.infinite import Mixture, solve
from fission.pincell.vera import Component
from fission.processing.groupr import (
    GroupLibrary,
    _angular_for_rust,
    group_library,
    two_body_kernel,
    two_body_kernel_table,
)
from fission.processing.groups import XMAS_172, GroupStructure, equal_lethargy
from fission.processing.thermal import FreeGas, Tabulated, thermal_kernels
from fission.processing.weights import ThermalSlowingDownFission

GROUPS = equal_lethargy(30, 1e-5, 2e7)


def _same(rust, python, rel=1e-11):
    scale = max(np.max(np.abs(python)), 1e-300)
    assert np.max(np.abs(np.asarray(rust) - np.asarray(python))) <= rel * scale


def _kernel_both(dist, awr, q, energy, groups=GROUPS, lmax=3):
    py = two_body_kernel(dist, awr, q, groups, lmax, energy)
    rs = _rust.two_body_kernel(
        *_angular_for_rust(dist), awr, q, list(groups.boundaries), lmax, energy
    ).reshape(lmax + 1, groups.size)
    return rs, py


def _synthetic():
    grid = IncidentGrid.from_tab2([1e3, 1e6, 1e6, 2e7], [4], [2])
    panels = (
        LegendrePanel((0.1, 0.02)),
        LegendrePanel((0.4, 0.2, 0.05, 0.01)),
        TabulatedPanel(Tabulated1D([-1.0, 0.0, 0.9, 1.0], [0.1, 0.3, 1.0, 2.8])),
        TabulatedPanel(Tabulated1D([-1.0, 0.5, 1.0], [0.05, 0.2, 3.0])),
    )
    return AngularDistribution(2, CM, grid, panels)


@pytest.mark.parametrize("n", [1, 2, 5, 16, 33, 64, 128])
def test_gauss_rules_match_numpy(n):
    x, w = _rust.gauss_legendre(n)
    xr, wr = np.polynomial.legendre.leggauss(n)
    assert np.max(np.abs(x - xr)) < 1e-14
    assert np.max(np.abs(w - wr)) < 1e-14


@pytest.mark.parametrize("energy", [10.0, 5e4, 9.9e5, 1e6, 3e6, 1.5e7])
@pytest.mark.parametrize("q", [0.0, -2e5])
def test_two_body_kernel_synthetic(energy, q):
    rs, py = _kernel_both(_synthetic(), 15.9, q, energy)
    _same(rs, py)
    iso = AngularDistribution(2, CM, None, ())
    _same(*_kernel_both(iso, 1.0, 0.0, energy))


@pytest.fixture(scope="module")
def lib():
    lib = EndfLibrary()
    if not (lib.root / "archives" / "neutrons-version.VIII.1.tar.gz").exists():
        pytest.skip("ENDF/B-VIII.1 neutron data not downloaded")
    return lib


def test_two_body_kernel_real_data(lib):
    u238 = Evaluation.from_file(lib.neutron("U238"))
    ad = angular_distribution(u238, 2)  # Legendre below, tabulated above (LTT = 3)
    for e in (1.0, 2e4, 1e6, 8e6, 1.9e7):
        _same(*_kernel_both(ad, u238.header.awr, 0.0, e, XMAS_172))
    u235 = Evaluation.from_file(lib.neutron("U235"))
    level = energy_angle(u235, 51).products[0].distribution  # MF6 LAW = 2
    q = u235.cross_section(51).qi
    for e in (1e5, 1e6, 1e7):
        _same(*_kernel_both(level, u235.header.awr, q, e, XMAS_172))


def test_kernel_tables(lib):
    zr = Evaluation.from_file(lib.neutron("Zr90"))
    cs = zr.cross_section(51)
    ad = angular_distribution(zr, 51)
    seeds = [3e6, 5e6, 8e6, 1.2e7, 2e7]
    k_py = two_body_kernel_table(
        ad, zr.header.awr, cs.qi, GROUPS, 1, seeds, cs.table, False, False
    )
    k_rs = two_body_kernel_table(
        ad, zr.header.awr, cs.qi, GROUPS, 1, seeds, cs.table, False, True
    )
    assert k_rs.z == k_py.z
    _same(k_rs.values, k_py.values)
    el = angular_distribution(zr, 2)
    a_py = two_body_kernel_table(
        el,
        zr.header.awr,
        0.0,
        GROUPS,
        1,
        [1.0, 1e4, 1e6, 2e7],
        adaptive=True,
        use_rust=False,
    )
    a_rs = two_body_kernel_table(
        el,
        zr.header.awr,
        0.0,
        GROUPS,
        1,
        [1.0, 1e4, 1e6, 2e7],
        adaptive=True,
        use_rust=True,
    )
    assert a_rs.z == a_py.z  # same refinement decisions
    _same(a_rs.values, a_py.values)


def test_thermal_kernels():
    from fission.endf.distributions.sab import scattering_law

    lib = EndfLibrary()
    models = [FreeGas(0.99917, 565.0), FreeGas(236.0, 565.0)]
    if (lib.root / "archives" / "thermal_scatt-version.VIII.1.tar.gz").exists():
        models.append(
            Tabulated(
                scattering_law(Evaluation.from_file(lib.thermal("HinH2O")), 565.0)
            )
        )
    energies = [1e-4, 0.0253, 0.5, 3.9]
    for model in models:
        py = thermal_kernels(model, XMAS_172, 3, energies, use_rust=False)
        rs = thermal_kernels(model, XMAS_172, 3, energies, use_rust=True)
        for (kp, tp), (kr, tr) in zip(py, rs, strict=True):
            _same(kr, kp, 1e-12)
            assert tr == pytest.approx(tp, rel=1e-12)


def test_infinite_medium_power_iteration_matches_eigen_solve():
    rng = np.random.default_rng(3)
    g = 12
    groups = equal_lethargy(g, 1e-5, 2e7)
    comps, libs = [], []
    for i in range(3):
        scatter = np.triu(rng.uniform(0, 0.3, (g, g)))  # downscatter
        scatter += np.tril(rng.uniform(0, 0.01, (g, g)), -1)  # some upscatter
        total = scatter.sum(axis=1) + rng.uniform(0.05, 0.4, g)
        chi = np.zeros(g)
        chi[:4] = rng.uniform(0, 1, 4)
        libs.append(GroupLibrary(
            "X", 0.0, groups, "", 0, np.ones(g), {1: total}, scatter[None],
            np.zeros((1, g, g)), nu_fission=rng.uniform(0, 0.5, g) * (i < 2),
            chi=chi / chi.sum(),
        ))  # fmt: skip
        comps.append(Component(f"X{i}", "fuel", 0.0, rng.uniform(0.5, 1.5), "free-gas"))
    mix = Mixture(comps, libs)
    eig, power = solve(mix, use_rust=False), solve(mix, use_rust=True)
    assert power.k == pytest.approx(eig.k, rel=1e-11)
    assert power.flux == pytest.approx(eig.flux, rel=1e-8)


def test_group_library_rust_matches_python(lib):
    ev = Evaluation.from_file(lib.neutron("Zr90"))
    from fission.endf.cache import cached_pointwise

    table = cached_pointwise("Zr90")
    groups = GroupStructure.from_ascending(
        "coarse", [1e-5, 1.0, 1e3, 1e5, 1e6, 5e6, 2e7]
    )
    w = ThermalSlowingDownFission.at(293.6)
    py = group_library(ev, table, groups, w, lmax=1, use_rust=False)
    rs = group_library(ev, table, groups, w, lmax=1, use_rust=True)
    for name in ("scatter_elastic", "scatter_inelastic"):
        _same(getattr(rs, name), getattr(py, name), 1e-9)
