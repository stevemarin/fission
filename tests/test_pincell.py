"""Tests for the infinite-medium solver (the VERA model is tested in test_model)."""

import numpy as np
import pytest

from fission.model import Component
from fission.pincell.infinite import Mixture, solve
from fission.processing.groupr import GroupLibrary
from fission.processing.groups import GroupStructure


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
