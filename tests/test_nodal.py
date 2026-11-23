"""Nodal diffusion (M4): mesh, finite difference, NEM."""

import numpy as np
import pytest

from fission.model import load_model
from fission.nodal import diffusion, nem
from fission.nodal.mesh import core_mesh

D, SIGMA_A, NU_SIGMA_F, C = 1.2, 0.08, 0.1, 0.5

XS = """
[library]
format = 1
groups = 1
[materials.f]
diffusion = [1.2]
absorption = [0.08]
total = [0.08]
scatter = [[0.0]]
nu_fission = [0.1]
chi = [1.0]
"""

CUBE = """
[model]
format = 1
name = "cube"
multigroup = "xs.toml"
[materials.f]
macroscopic = "f"
[lattices.f]
material = "f"
[assemblies.a]
axial = [{ lattice = "f", height = 60.0 }]
[core]
assembly_pitch = 30.0
boundary = { left = "reflective", top = "reflective", right = "vacuum", bottom = "vacuum" }
legend = { A = "a" }
map = '''
A A
A A
'''
"""


@pytest.fixture
def cube(tmp_path):
    (tmp_path / "xs.toml").write_text(XS)
    (tmp_path / "cube.toml").write_text(CUBE)
    return load_model(tmp_path / "cube.toml")


def buckling(half_width: float) -> float:
    """B with B tan(B a) = c / D: cos(B x) between a reflective plane at 0
    and a vacuum face at a with d phi / dn = -(c / D) phi."""
    lo, hi = 1e-12, np.pi / (2 * half_width) - 1e-12
    for _ in range(200):
        mid = (lo + hi) / 2
        lo, hi = (mid, hi) if mid * np.tan(mid * half_width) < C / D else (lo, mid)
    return mid


def exact_k(*bucklings: float) -> float:
    return NU_SIGMA_F / (SIGMA_A + D * sum(b * b for b in bucklings))


def test_finite_difference_converges_quadratically(cube):
    """60 x 60 cm (reflective on two sides) x 60 cm high (vacuum both ends)."""
    k = exact_k(buckling(60.0), buckling(60.0), buckling(30.0))
    errors = [diffusion.solve(core_mesh(cube, refine=r), tol_k=1e-10,
                              tol_source=1e-9, inner=10).k - k for r in (2, 4)]  # fmt: skip
    assert errors[0] > 0 and errors[1] > 0
    assert errors[0] / errors[1] == pytest.approx(4.0, rel=0.05)


def test_two_node_problem_reproduces_a_cosine():
    k = 1.0
    b = np.sqrt((NU_SIGMA_F / k - SIGMA_A) / D)
    h, x0 = 20.0, 7.0

    def f(x):
        return np.cos(b * (x - x0))

    def average(a, c):
        return (np.sin(b * (c - x0)) - np.sin(b * (a - x0))) / (b * (c - a))

    def current(x):  # J = -D phi'
        return D * b * np.sin(b * (x - x0))

    def node(a, c, j_far):
        one = np.ones((1, 1))
        return {"d": D * one, "removal": SIGMA_A * one, "c": NU_SIGMA_F / k * np.ones((1, 1, 1)),
                "h": np.array([h]), "adf": one, "phi": average(a, c) * one,
                "j_far": j_far * one, "l1": 0 * one, "l2": 0 * one}  # fmt: skip

    # 20 cm nodes: 2.6 rad of the cosine per node, the current is still within 0.2%
    _, j = nem.two_node(node(-20, 0, current(-20)), node(0, 20, current(20)), 1)
    assert j[0, 0] == pytest.approx(current(0.0), rel=2e-3)
    h = 10.0  # and the error falls fast with the node size
    phi_s, j = nem.two_node(node(-10, 0, current(-10)), node(0, 10, current(10)), 1)
    assert phi_s[0, 0] == pytest.approx(f(0.0), rel=1e-4)
    assert j[0, 0] == pytest.approx(current(0.0), rel=1e-4)


def test_leakage_shape_is_exact_for_a_quadratic(cube):
    mesh = core_mesh(cube, refine=3)  # 10 cm nodes in x
    x = np.cumsum(mesh.widths[0]) - mesh.widths[0] / 2

    def q(t):
        return 2.0 + 0.3 * t - 0.004 * t * t

    def avg(lo, hi):  # average of q over [lo, hi]
        def integral(t):
            return 2.0 * t + 0.15 * t * t - 0.004 * t**3 / 3

        return (integral(hi) - integral(lo)) / (hi - lo)

    lbar = np.zeros((1, *mesh.shape))
    lbar[...] = np.array([avg(c - 5, c + 5) for c in x])
    l1, l2 = nem.leakage_shape(mesh, lbar, 0)
    i = 2  # an interior node: q = L + l1 P1 + l2 P2 in node units (h = 10)
    assert l1[0, 0, 0, i] == pytest.approx(10 * (0.3 - 0.008 * x[i]), rel=1e-10)
    assert l2[0, 0, 0, i] == pytest.approx(-0.004 * 100 / 3, rel=1e-10)


def test_nem_is_exact_in_one_dimension(cube):
    mesh = core_mesh(cube, refine=1)  # two 30 cm nodes along x
    mesh.boundary = dict.fromkeys(mesh.boundary, "reflective")
    mesh.boundary["east"] = "vacuum"
    k = exact_k(buckling(60.0))
    fd = diffusion.solve(mesh, tol_k=1e-10, tol_source=1e-9).k
    r = nem.solve(mesh, tol_k=1e-9, tol_flux=1e-8)
    assert abs(r.k - k) < 1e-6 < abs(fd - k) / 100


def test_nem_converges_in_two_dimensions(cube):
    k = exact_k(buckling(60.0), buckling(60.0))
    errors = []
    for refine in (1, 2):
        mesh = core_mesh(cube, refine=refine)
        mesh.boundary["bottom"] = mesh.boundary["top"] = "reflective"
        errors.append(abs(nem.solve(mesh).k - k))
    assert errors[1] < 1e-4 and errors[1] < errors[0] / 5


def test_iaea_mesh():
    m = load_model("iaea-3d")
    assert m.multigroup["fuel1"].diffusion_coefficient.tolist() == [1.5, 0.4]
    mesh = core_mesh(m, refine=1)
    assert mesh.shape == (19, 9, 9)  # 20 cm axial nodes; 9 x 9 quarter positions
    assert mesh.widths[0][0] == 10.0 and mesh.widths[0][1] == 20.0
    assert mesh.boundary["west"] == mesh.boundary["north"] == "reflective"
    names = np.array(mesh.names)
    # the centre column: reflector at the bottom, rodded fuel, rodded reflector on top
    column = names[mesh.material[:, 0, 0]]
    assert (
        column[0] == "reflector"
        and column[1] == "fuel2_rod"
        and column[-1] == "reflector_rod"
    )
    # assembly 18 (map row 2, column 2): rod only from 280 to 360 cm
    partial = names[mesh.material[:, 2, 2]]
    assert partial[13] == "fuel2" and partial[14] == "fuel2_rod"
    assert (mesh.material[:, 8, 4:] == -1).all()  # outside the core


def test_nem_on_the_iaea_3d_benchmark():
    """One node per assembly: within 15 pcm of NEM on 5 cm nodes (1.0290742,
    docs ch. 21), against 250 pcm for finite difference on the same mesh."""
    mesh = core_mesh(load_model("iaea-3d"), refine=1)
    r = nem.solve(mesh)
    assert abs(r.k - 1.0290742) < 1.5e-4
    fd = diffusion.solve(mesh, tol_k=1e-8, tol_source=1e-7).k
    assert abs(fd - 1.0290742) > 2e-3
    from fission.nodal.__main__ import assembly_powers

    p = assembly_powers(mesh, r.flux)
    assert np.nanmax(p) == pytest.approx(1.438, abs=2e-3)
    assert p[0, 0] < p[1, 1]  # the rodded centre


def test_nem_warm_start_gives_the_same_answer():
    mesh = core_mesh(load_model("iaea-3d"), refine=1)
    cold = nem.solve(mesh)
    warm = nem.solve(mesh, start=cold)
    assert warm.k == pytest.approx(cold.k, abs=2e-7)
    assert warm.iterations < cold.iterations


def test_control_rod_banks_cut_the_stacks():
    """wbn1-rods: 57 RCCAs in the full core; a bank inserted splits its
    assemblies at the absorber zones and uses the rodded lattices."""
    from fission.nodal.mesh import stacks

    m = load_model("wbn1-rods")
    banks = m.core.banks
    full = sum(
        1 if (r, c) == (0, 0) else 2 if 0 in (r, c) else 4
        for cells in banks.positions.values()
        for r, c in cells
    )
    assert full == 57
    out = stacks(m)  # all rods out: the stacks of the assemblies
    assert out[0, 0] == m.assemblies["a"].segments
    rod = m.control_rods["rcca"]
    for steps in (0, 100):
        segs = stacks(m, {"D": steps})[0, 0]
        assert sum(s.height for s in segs) == pytest.approx(m.assemblies["a"].height)
        names = [s.lattice for s in segs]
        assert names == [
            "reflector",
            "fuel_211",
            "fuel_211_aic",
            "fuel_211_b4c",
            "reflector",
        ]
        assert segs[1].height == pytest.approx(rod.bottom - 21.5 + steps * rod.step)
        assert segs[2].height == pytest.approx(101.6)
    # fully inserted, the B4C ends at the top of the fuel
    assert stacks(m, {"D": 0})[0, 0][3].height == pytest.approx(259.08)
    # fully withdrawn: the tip is in the upper reflector, nothing rodded
    assert stacks(m, {"D": rod.steps})[0, 0] == m.assemblies["a"].segments
    # other banks' positions are untouched
    assert stacks(m, {"D": 0})[0, 2] == m.assemblies["a"].segments
    with pytest.raises(ValueError, match="no bank"):
        stacks(m, {"X": 0})
