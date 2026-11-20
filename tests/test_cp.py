"""Collision probabilities of an annular cell (step C4)."""

import math

import numpy as np
import pytest

from fission.transport.cp import bickley, collision_probabilities, eigenvalue


def test_bickley_functions():
    assert bickley(1, np.array(0.0)) == pytest.approx(math.pi / 2, rel=1e-12)
    assert bickley(3, np.array(0.0)) == pytest.approx(math.pi / 4, rel=1e-12)
    assert bickley(2, np.array(0.0)) == pytest.approx(1.0, rel=1e-12)
    # dKi_n/dx = -Ki_(n-1)
    for x in (0.05, 0.5, 3.0):
        h = 1e-5
        slope = (bickley(3, np.array(x + h)) - bickley(3, np.array(x - h))) / (2 * h)
        assert slope == pytest.approx(-bickley(2, np.array(x)), rel=1e-7)


RADII = np.array([0.24, 0.34, 0.4096, 0.418, 0.475, 0.6, 1.26 / math.sqrt(math.pi)])
SIGMA = np.array(
    [[0.6, 0.6, 0.6, 1e-4, 0.3, 1.2, 1.2], [9.0, 9.0, 9.0, 1e-4, 0.3, 2.5, 2.5]]
)


def test_conservation_reciprocity_and_volumes():
    cp = collision_probabilities(RADII, SIGMA)
    np.testing.assert_allclose(cp.collide.sum(axis=2) + cp.escape, 1.0, atol=1e-12)
    exact = math.pi * np.diff(np.concatenate([[0.0], RADII**2]))
    np.testing.assert_allclose(cp.volumes, exact, rtol=2e-6)
    for g in range(2):  # V_i Sigma_i P_ij is symmetric
        m = cp.volumes[:, None] * SIGMA[g][:, None] * cp.collide[g]
        np.testing.assert_allclose(m, m.T, rtol=1e-6, atol=1e-14)
    white = cp.white()
    np.testing.assert_allclose(white.sum(axis=2), 1.0, atol=1e-10)


def test_escape_from_a_bare_cylinder():
    """Limits: transparent (P_esc -> 1) and black (P_esc -> S / (4 V Sigma))."""
    r = np.array([0.4])
    thin = collision_probabilities(r, np.array([[1e-6]]))
    assert thin.escape[0, 0] == pytest.approx(1.0, abs=1e-5)
    thick = collision_probabilities(r, np.array([[200.0]]))
    surface, volume = 2 * math.pi * 0.4, math.pi * 0.16
    assert thick.escape[0, 0] * 4 * volume * 200 / surface == pytest.approx(
        1.0, rel=0.02
    )


def test_homogeneous_cell_gives_k_infinity():
    """Every ring the same material, white boundary: k = k_inf of the medium."""
    n = len(RADII)
    total = np.array([0.5, 1.2])
    scatter = np.tile(np.array([[0.4, 0.08], [0.0, 1.1]]), (n, 1, 1))
    nu = np.tile(np.array([0.005, 0.15]), (n, 1))
    chi = np.tile(np.array([1.0, 0.0]), (n, 1))
    cp = collision_probabilities(RADII, np.tile(total[:, None], (1, n)))
    k, phi = eigenvalue(cp, scatter, nu, chi)
    a = np.diag(total) - scatter[0].T
    k_inf = nu[0] @ np.linalg.solve(a, chi[0])
    assert k == pytest.approx(k_inf, rel=1e-9)
    np.testing.assert_allclose(phi / phi[:, :1], 1.0, rtol=1e-8)  # flat in space
