import sys
from pathlib import Path

import numpy as np
import pytest

from nodal import Core, Material, NodalSolver, nem

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "benchmarks"))
import iaea  # type: ignore

FUEL = Material.two_group(1.5, 0.4, 0.01, 0.08, 0.02, 0.0, 0.135)


def k_bare(mat, B2):
    D1, D2 = mat.D
    sr1, sa2 = mat.sigma_r[0], mat.sigma_a[1]
    return mat.nu_sigma_f[1] * mat.sigma_s[0, 1] / ((D1 * B2 + sr1) * (D2 * B2 + sa2))


# ---------------------------------------------------------------- NEM kernels
def _cos_moments(D, sr, h, B, x0):
    """Exact node average and source moments for phi = cos(B x) on [x0-h/2, x0+h/2]."""
    xi = np.linspace(-0.5, 0.5, 40001)
    w = np.full_like(xi, 1 / (xi.size - 1))
    w[[0, -1]] *= 0.5
    phi = np.cos(B * (x0 + h * xi))
    S = (D * B**2 + sr) * phi
    return (w * phi).sum(), (w * xi * S).sum(), (w * (3 * xi**2 - 0.25) * S).sum()


def test_one_node_reproduces_cosine():
    D, sr, h, B, x0, f = 1.3, 0.05, 20.0, 0.03, 4.0, 1.0
    phibar, S1, S2 = _cos_moments(D, sr, h, B, x0)
    xm, xp = x0 - h / 2, x0 + h / 2
    phi = lambda x: np.cos(B * x)
    J = lambda x: D * B * np.sin(B * x)
    jin_m = phi(xm) / 4 + J(xm) / 2   # J+ entering the - face
    jin_p = phi(xp) / 4 - J(xp) / 2   # J- entering the + face
    o = np.ones(1)
    _, jo_m, jo_p = nem.one_node(D * o, h * o, sr * o, f * o, phibar * o, S1 * o, S2 * o,
                                 jin_m * o, jin_p * o)
    assert jo_m[0] == pytest.approx(phi(xm) / 4 - J(xm) / 2, rel=1e-5)
    assert jo_p[0] == pytest.approx(phi(xp) / 4 + J(xp) / 2, rel=1e-5)


def test_leakage_fit_is_exact_for_quadratics():
    # L(x) = 1 + 2x + 3x^2 on nodes of widths 4 (centre), 2 (left), 6 (right)
    h0, hm, hp = 4.0, 2.0, 6.0
    avg = lambda a, b: (b - a + (b**2 - a**2) + (b**3 - a**3)) / (b - a)
    L0, Lm, Lp = avg(-2, 2), avg(-4, -2), avg(2, 8)
    L1, L2 = nem.leakage_moments(L0, Lm, Lp, -0.5 - hm / h0, -0.5, 0.5, 0.5 + hp / h0)
    xi = np.linspace(-0.5, 0.5, 20001)
    w = np.full_like(xi, 1 / (xi.size - 1))
    w[[0, -1]] *= 0.5
    L = 1 + 2 * h0 * xi + 3 * (h0 * xi) ** 2
    assert L1 == pytest.approx((w * xi * L).sum(), rel=1e-6)
    assert L2 == pytest.approx((w * (3 * xi**2 - 0.25) * L).sum(), rel=1e-6)


# ------------------------------------------------------------ analytic cores
def test_infinite_medium():
    core = Core([[10.0] * 3, [10.0] * 3, [1.0]], np.zeros((3, 3, 1), int), [FUEL], {})
    r = NodalSolver(core).solve()
    assert r.k == pytest.approx(k_bare(FUEL, 0.0), abs=1e-7)


@pytest.mark.parametrize("n", [4, 10])
def test_bare_slab(n):
    a = 200.0
    core = Core([np.full(n, a / n), [1.0], [1.0]], np.zeros((n, 1, 1), int), [FUEL],
                {"x-": "zero_flux", "x+": "zero_flux"})
    kref = k_bare(FUEL, (np.pi / a) ** 2)
    nem_k = NodalSolver(core).solve().k
    fdm_k = NodalSolver(core, "fdm").solve().k
    assert abs(nem_k - kref) < 1e-6           # < 0.1 pcm, even at 50 cm nodes
    assert abs(fdm_k - kref) > 10 * abs(nem_k - kref)


def test_bare_square_exercises_transverse_leakage():
    a, n = 160.0, 8
    core = Core([np.full(n, a / n), np.full(n, a / n), [1.0]], np.zeros((n, n, 1), int),
                [FUEL], {s: "zero_flux" for s in ("x-", "x+", "y-", "y+")})
    kref = k_bare(FUEL, 2 * (np.pi / a) ** 2)
    assert NodalSolver(core).solve().k == pytest.approx(kref, abs=5e-6)  # error is 0.35 pcm at 20 cm


def test_reflective_boundary_equals_mirrored_core():
    m2 = Material.two_group(2.0, 0.3, 0.0, 0.01, 0.04, 0.0, 0.0)  # reflector
    half = [0, 0, 0, 1, 1]
    full = half[::-1] + half
    bc = {"x+": "vacuum"}
    core_h = Core([[20.0] * 5, [1.0], [1.0]], np.array(half)[:, None, None], [FUEL, m2], bc)
    core_f = Core([[20.0] * 10, [1.0], [1.0]], np.array(full)[:, None, None], [FUEL, m2],
                  {"x-": "vacuum", "x+": "vacuum"})
    tight = {"tol_k": 1e-10, "tol_phi": 1e-9}
    k_h = NodalSolver(core_h).solve(**tight).k
    k_f = NodalSolver(core_f).solve(**tight).k
    assert k_h == pytest.approx(k_f, abs=1e-9)


def test_uniform_adf_has_no_effect():
    # interior continuity of f*phi is unchanged by a uniform f; zero-flux edges
    # keep the boundary condition independent of f too
    a, n = 200.0, 5
    mk = lambda adf: Core([np.full(n, a / n), [1.0], [1.0]], np.zeros((n, 1, 1), int),
                          [Material.two_group(1.5, 0.4, 0.01, 0.08, 0.02, 0.0, 0.135, adf=adf)],
                          {"x-": "zero_flux", "x+": "zero_flux"})
    k1 = NodalSolver(mk(None)).solve().k
    k2 = NodalSolver(mk([1.1, 0.9])).solve().k
    assert k1 == pytest.approx(k2, abs=1e-7)


# ------------------------------------------------------------- IAEA benchmark
def test_iaea_2d():
    core, ox, oy, _ = iaea.build(2, 20.0)
    r = NodalSolver(core).solve()
    assert r.converged
    assert abs(r.k - iaea.K_REF[2]) < 15e-5
    assert abs(NodalSolver(core, "fdm").solve().k - iaea.K_REF[2]) > 200e-5
    p = iaea.assembly_power(core, r, ox, oy)
    assert np.allclose(p, p.T, equal_nan=True, atol=1e-6)  # octant symmetry


def test_iaea_3d():
    core, *_ = iaea.build(3, 20.0)
    r = NodalSolver(core).solve()
    assert r.converged
    assert abs(r.k - iaea.K_REF[3]) < 15e-5
