"""Tests for the VERA pin-cell data and the infinite-medium solver."""

import math

import numpy as np
import pytest

from fission.pincell.infinite import Mixture, solve
from fission.pincell.vera import CASES, CLAD, FUEL, PITCH, homogenized, region_areas
from fission.processing.groupr import GroupLibrary
from fission.processing.groups import GroupStructure


def test_region_areas_fill_the_cell():
    areas = region_areas()
    assert sum(areas.values()) == pytest.approx(PITCH**2, rel=1e-14)
    assert areas["fuel"] == pytest.approx(math.pi * 0.4096**2)


def test_homogenization_conserves_atoms():
    case = CASES["1A"]
    comps = homogenized(case, include_unsupported=True)
    areas = region_areas()
    u238 = sum(c.density for c in comps if c.nuclide == "U238") * PITCH**2
    assert u238 == pytest.approx(FUEL["U238"] * areas["fuel"])
    o16 = {c.region: c for c in comps if c.nuclide == "O16"}
    assert set(o16) == {"fuel", "moderator"}  # two components, two temperatures
    assert o16["moderator"].thermal == "free-gas"
    (h,) = [c for c in comps if c.nuclide == "H1"]
    assert h.thermal == "HinH2O"
    assert len(homogenized(case)) == len(comps) - 4  # Fe-54, -57, -58, Cr-53 left out
    left_out = sum(CLAD[n] for n in ("Fe54", "Fe57", "Fe58", "Cr53"))
    assert left_out < 2e-5
    assert left_out / sum(CLAD.values()) < 5e-4


def _two_group(sigma_t, scatter, nu_fission, chi):
    groups = GroupStructure("two", (2e7, 0.625, 1e-5))
    fission = np.asarray(nu_fission) / 2.5
    removal = np.asarray(sigma_t) - np.asarray(scatter).sum(axis=1)
    sigma = {1: np.asarray(sigma_t), 18: fission, 102: removal - fission}
    lib = GroupLibrary(
        "X", 0.0, groups, "", 0, np.ones(2), sigma,
        np.asarray(scatter)[None], np.zeros((1, 2, 2)),
        nu_fission=np.asarray(nu_fission), chi=np.asarray(chi),
    )  # fmt: skip
    from fission.pincell.vera import Component

    return Mixture([Component("X", "fuel", 0.0, 1.0, "free-gas")], [lib])


def test_two_group_k_infinity_matches_the_formula():
    st1, st2, s11, s12, s22 = 0.52, 1.30, 0.48, 0.02, 1.20
    nf1, nf2 = 0.008, 0.16
    mix = _two_group([st1, st2], [[s11, s12], [0.0, s22]], [nf1, nf2], [1.0, 0.0])
    sol = solve(mix)
    expected = (nf1 + nf2 * s12 / (st2 - s22)) / (st1 - s11)
    assert sol.k == pytest.approx(expected, rel=1e-12)
    # neutron balance: production / removal = k
    assert (mix.nu_fission @ sol.flux) / (mix.removal @ sol.flux) == pytest.approx(
        sol.k
    )


def test_four_factors_multiply_to_k():
    mix = _two_group(
        [0.52, 1.30], [[0.48, 0.02], [0.001, 1.20]], [0.008, 0.16], [1.0, 0.0]
    )
    sol = solve(mix)
    ff = sol.four_factors()
    assert ff["product"] == pytest.approx(sol.k, rel=1e-12)
    assert 0 < ff["p"] < 1 and 0 < ff["f"] <= 1 and ff["epsilon"] >= 1
