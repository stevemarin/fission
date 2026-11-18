"""Tests for thermal scattering: MF7 reader, kernels, thermal group matrices."""

import numpy as np
import pytest

from fission.data import EndfLibrary
from fission.endf import Evaluation
from fission.endf.distributions.sab import available_temperatures, scattering_law
from fission.processing.groupr import GroupLibrary, Kernel, group_integrals
from fission.processing.groups import XMAS_172, equal_lethargy
from fission.processing.thermal import (
    FreeGas,
    Tabulated,
    ddxs,
    free_gas_total,
    thermal_kernel,
    with_thermal,
)
from fission.processing.weights import ThermalSlowingDownFission


@pytest.mark.parametrize("awr", [0.99917, 15.858, 236.0])
@pytest.mark.parametrize("energy", [0.01, 0.0253, 0.5, 4.0])
def test_free_gas_kernel_integrates_to_the_closed_form(awr, energy):
    """The free-gas kernel, integrated over E' and mu, is the Doppler-broadened
    cross section of a constant sigma (chapter 8)."""
    model = FreeGas(awr, 293.6)
    _, total = thermal_kernel(model, XMAS_172, 0, energy)
    expected = (
        free_gas_total(awr, 293.6, energy) * model.sigma_bound * (awr / (awr + 1)) ** 2
    )
    assert total == pytest.approx(expected, rel=1e-4)


def test_free_gas_mean_cosine_tends_to_two_thirds_a():
    model = FreeGas(15.858, 293.6)
    k, _ = thermal_kernel(model, XMAS_172, 1, 4.0)
    assert k[1].sum() / k[0].sum() == pytest.approx(2 / (3 * 15.858), rel=1e-3)


@pytest.fixture(scope="module")
def water():
    lib = EndfLibrary()
    if not (lib.root / "archives" / "thermal_scatt-version.VIII.1.tar.gz").exists():
        pytest.skip("ENDF/B-VIII.1 thermal scattering data not downloaded")
    return Evaluation.from_file(lib.thermal("HinH2O"))


@pytest.fixture(scope="module")
def law(water):
    return scattering_law(water, 293.6)


def test_reader(water, law):
    temps = available_temperatures(water)
    assert 293.6 in temps and 600.0 in temps
    assert law.count == 2  # two H per H2O
    assert law.sigma_bound == pytest.approx(4 * law.sigma_free, rel=2e-3)  # A ~ 1
    with pytest.raises(KeyError):
        scattering_law(water, 1200.0)


def test_detailed_balance_pointwise(law):
    rng = np.random.default_rng(1)
    for model in (FreeGas(1.0, 293.6), Tabulated(law)):
        kt = model.kt
        for _ in range(20):
            e, ep, mu = (
                rng.uniform(1e-3, 1.0),
                rng.uniform(1e-3, 1.0),
                rng.uniform(-1, 1),
            )
            down = ddxs(model, e, np.array([ep]), np.array([mu]))[0]
            up = ddxs(model, ep, np.array([e]), np.array([mu]))[0]
            assert e * np.exp(-e / kt) * down == pytest.approx(
                ep * np.exp(-ep / kt) * up, rel=1e-12
            )


def test_bound_hydrogen(law):
    """Water scatters more than free hydrogen at thermal energies (binding),
    and tends to the free-atom value at a few eV (impulse approximation)."""
    model = Tabulated(law)
    free = model.law.sigma_free
    _, high = thermal_kernel(model, XMAS_172, 0, 4.0)
    assert high == pytest.approx(free_gas_total(0.99917, 293.6, 4.0) * free, rel=0.01)
    k, thermal = thermal_kernel(model, XMAS_172, 1, 0.0253)
    assert thermal > 1.5 * free_gas_total(0.99917, 293.6, 0.0253) * free
    assert 0 < k[1].sum() / k[0].sum() < 2 / 3


def test_group_matrices_obey_detailed_balance(law):
    """With a Maxwellian weight, phi_g S(g -> g') = phi_g' S(g' -> g)."""
    groups = equal_lethargy(8, 1e-4, 0.045)  # below 2 kT: weight is Maxwellian
    weight = ThermalSlowingDownFission.at(293.6)
    flux = group_integrals(
        [1e-4, 0.045], [[1.0, 1.0]], weight, groups, Kernel.constant(1e-4, 0.045), False
    )[:, 0, 0]
    zeros = np.zeros(groups.size)
    lib = GroupLibrary(
        "H1", 293.6, groups, "", 0, flux, {1: zeros.copy(), 2: zeros.copy()},
        np.zeros((1, 8, 8)), np.zeros((1, 8, 8)),
    )  # fmt: skip
    model = Tabulated(law)
    out = with_thermal(lib, model, weight, cutoff=0.045, use_rust=False)
    rate = flux[:, None] * out.scatter_elastic[0]
    big = rate > 1e-3 * rate.max()
    # second-order convergence in the kernel spacing: ~0.1% for du = 0.02
    assert rate[big] == pytest.approx(rate.T[big], rel=2e-3)
    assert np.tril(out.scatter_elastic[0], -1).sum() > 0  # upscattering exists


def test_library_with_free_gas_conserves_the_elastic_cross_section():
    lib = EndfLibrary()
    if not (lib.root / "archives" / "neutrons-version.VIII.1.tar.gz").exists():
        pytest.skip("ENDF/B-VIII.1 neutron data not downloaded")
    from fission.endf.cache import cached_pointwise
    from fission.processing.groupr import group_library

    ev = Evaluation.from_file(lib.neutron("O16"))
    table = cached_pointwise("O16")
    groups = equal_lethargy(12, 1e-5, 2e7)
    weight = ThermalSlowingDownFission.at(293.6)
    gl = group_library(ev, table, groups, weight, lmax=1)
    hot = with_thermal(
        gl, FreeGas(ev.header.awr, 293.6), weight, (table.energies, table.xs[2])
    )
    rows = [
        g
        for g in range(groups.size)
        if groups.bounds(g)[1] <= 4.0 and groups.bounds(g)[0] > 1e-3
    ]
    for g in rows:
        assert hot.scatter_elastic[0][g].sum() == pytest.approx(
            gl.sigma[2][g], rel=1e-4
        )
    assert np.tril(hot.scatter_elastic[0], -1).sum() > 0
    assert (
        np.tril(gl.scatter_elastic[0], -1).sum() < 1e-12
    )  # target at rest: rounding only
