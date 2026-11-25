"""Thermal-hydraulic feedback (M4b): water properties, fuel rod, coupled core."""

import math

import numpy as np
import pytest

from fission.feedback import fuel, water
from fission.model import BRANCH_AXES, load_model
from fission.tables import XSTable


@pytest.mark.parametrize(
    ("t", "p", "v", "h", "cp"),
    [  # IAPWS R7-97(2012), Table 5
        (300.0, 3e6, 0.100215168e-2, 0.115331273e3, 0.417301218e1),
        (300.0, 80e6, 0.971180894e-3, 0.184142828e3, 0.401008987e1),
        (500.0, 3e6, 0.120241800e-2, 0.975542239e3, 0.465580682e1),
    ],
)
def test_if97_region1_verification_values(t, p, v, h, cp):
    assert water.specific_volume(p, t) == pytest.approx(v, rel=1e-8)
    assert water.enthalpy(p, t) / 1e3 == pytest.approx(h, rel=1e-8)
    assert water.heat_capacity(p, t) / 1e3 == pytest.approx(cp, rel=1e-8)
    assert water.temperature(p, water.enthalpy(p, t)) == pytest.approx(t, abs=1e-8)


def test_fuel_rod():
    rp, ri, ro = 0.4096e-2, 0.418e-2, 0.475e-2
    cold = fuel.rod(0.0, 580.0, rp, ri, ro)
    assert cold.centre == pytest.approx(580.0) and cold.doppler == pytest.approx(580.0)
    q = 2e4
    r = fuel.rod(q, 580.0, rp, ri, ro)
    # conduction in the pellet, gap and clad drops, Doppler weighting
    integral = fuel.k_fuel_integral(r.centre) - fuel.k_fuel_integral(r.surface)
    assert integral == pytest.approx(q / (4 * math.pi), rel=1e-9)
    assert r.surface - r.clad_inner == pytest.approx(
        q / (2 * math.pi * rp * 1e4), rel=1e-12
    )
    k = fuel.k_clad(0.5 * (r.clad_inner + r.clad_outer))
    assert r.clad_inner - r.clad_outer == pytest.approx(q * math.log(ro / ri) / (2 * math.pi * k),
                                                        rel=1e-6)  # fmt: skip
    assert r.doppler == pytest.approx(0.3 * r.centre + 0.7 * r.surface)
    assert r.centre > fuel.rod(q / 2, 580.0, rp, ri, ro).centre


def _table(path, fuel_node: bool):
    """Two-group constants that depend on the state like a PWR: Doppler
    raises fast absorption, less water lowers moderation, boron absorbs."""
    axes: dict[str, tuple[float, ...]] = {
        "fuel_temperature": (500.0, 1500.0),
        "moderator_temperature": (550.0, 620.0),
        "moderator_density": (0.6, 0.8),
        "boron_ppm": (0.0, 2000.0),
    }
    shape = (2, 2, 2, 2)
    data = {k: np.zeros(shape + (2,)) for k in ("D", "total", "sigma_a", "nu_sigma_f",
                                                "sigma_f", "kappa_sigma_f", "chi", "adf")}  # fmt: skip
    data["sigma_s"] = np.zeros(shape + (2, 2))
    for index in np.ndindex(shape):
        tf, _, rho, ppm = (axes[a][i] for a, i in zip(BRANCH_AXES, index, strict=True))
        s12 = 0.02 * rho / 0.7
        a1 = (
            0.010 * (1 + 0.15 * (math.sqrt(tf) - math.sqrt(900)) / 30)
            if fuel_node
            else 0.0
        )
        a2 = (0.08 if fuel_node else 0.01) * rho / 0.7 + 1e-5 * ppm * rho / 0.7
        nsf = np.array([0.006, 0.12]) if fuel_node else np.zeros(2)
        data["D"][index] = [1.4, 0.4] if fuel_node else [1.3, 0.2]
        data["sigma_a"][index] = [a1, a2]
        data["sigma_s"][index] = [[0.0, s12], [0.0, 0.0]]
        data["total"][index] = [a1 + s12, a2]
        data["nu_sigma_f"][index] = nsf
        data["sigma_f"][index] = nsf / 2.4
        data["kappa_sigma_f"][index] = nsf / 2.4 * 2e8
        data["chi"][index] = [1.0, 0.0]
        data["adf"][index] = [1.0, 1.0]
    XSTable(axes, data, {}).write(path)


CORE = """
[model]
format = 1
name = "small"
[materials.fuel]
temperature = "fuel"
table = "{fuel}"
[materials.reflector]
table = "{reflector}"
[lattices.fuel]
material = "fuel"
[lattices.reflector]
material = "reflector"
[assemblies.fuel]
axial = [{{ lattice = "reflector", height = 20.0 }}, {{ lattice = "fuel", height = 200.0 }},
         {{ lattice = "reflector", height = 20.0 }}]
[assemblies.reflector]
axial = [{{ lattice = "reflector", height = 240.0 }}]
[core]
assembly_pitch = 20.0
symmetry = "quarter"
state = "hfp"
legend = {{ F = "fuel", R = "reflector" }}
map = '''
F F R
F F R
R R .
'''
[core.thermal]
assembly_power = 1.2e7
assembly_flow = 80.0
inlet_temperature = 565.0
pressure = 15.5e6
fuel_heat_fraction = 0.974
rods = 264
pellet_radius = 0.4096
clad_inner = 0.418
clad_outer = 0.475
gap_conductance = 1.0e4
[states.hfp]
fuel_temperature = 900.0
moderator_temperature = 565.0
moderator_density = 0.74
boron_ppm = 1000
"""


def test_core_with_feedback(tmp_path):
    from fission.feedback.core import solve

    _table(tmp_path / "fuel.json", True)
    _table(tmp_path / "reflector.json", False)
    (tmp_path / "core.toml").write_text(
        CORE.format(fuel=tmp_path / "fuel.json", reflector=tmp_path / "reflector.json")
    )
    model = load_model(tmp_path / "core.toml")
    r = solve(model, axial=20.0)
    assert r.history[-1][1] < 0.05 and r.history[-1][2] < 0.05  # converged
    th = model.core.thermal
    fuel_volume = (
        2.25 * 20.0 * 20.0 * 200.0
    )  # 1 + 2 x 1/2 + 1/4 assemblies in the quarter
    assert r.power.sum() == pytest.approx(
        th.assembly_power / (20 * 20 * 200) * fuel_volume
    )
    # energy balance of a full assembly channel (map position 1, 1)
    column = r.power[:, 1, 1].sum()
    h_out = (
        water.enthalpy(th.pressure, th.inlet_temperature) + column / th.assembly_flow
    )
    top = r.states["moderator_temperature"][-1, 1, 1] - 0  # reflector above the fuel
    assert water.enthalpy(th.pressure, top) == pytest.approx(h_out, rel=1e-9)
    # the hotter the rod, the hotter its fuel
    tf = r.states["fuel_temperature"][1:-1, 0, 0]
    p = r.power[1:-1, 0, 0]
    assert np.argmax(tf) == np.argmax(p)
    assert (
        r.states["moderator_density"][-1, 0, 0] < r.states["moderator_density"][0, 0, 0]
    )


def test_reflector_slab_water_cells():
    from fission.fewgroup.reflector import water_cells
    from fission.transport import Geometry

    m = load_model("vera-p1")
    g = Geometry.from_lattice(m, "reflector", {"fuel": 3}, 8, 0.6)
    assert (g.nx, g.ny) == (34, 1)
    water = water_cells(g)
    cell = g.cell_of_regions()
    assert set(cell[water]) == set(range(17, 34))
    assert set(cell[~water]) == set(range(17))


@pytest.mark.parametrize(
    ("nu_sigma_f", "k", "current"),
    [((0.0, 0.0), np.inf, (1.0, 0.1)),  # baffle-like reflector, current entering
     ((0.0065, 0.13), 1.0, (-1.0, -0.3))],  # fuel at the slab's k, current leaving
)  # fmt: skip
def test_homogeneous_surface_flux_matches_finite_differences(nu_sigma_f, k, current):
    """The modal solution of a homogeneous medium against a fine-mesh
    finite-difference one driven by a current at x = 0 (with a little
    up-scattering; in the fuel the fundamental mode is a cosine)."""
    from fission.fewgroup.condense import FewGroup
    from fission.fewgroup.reflector import homogeneous_surface_flux

    two = np.ones(2)
    s = np.array([[0.5956, 0.0270], [0.00079, 1.641]])
    fg = FewGroup("r", (2e7, 0.625, 1e-5), two, np.array([1.338, 0.307]),
                  two, np.array([0.00226, 0.0369]), s, np.array(nu_sigma_f), 0 * two,
                  0 * two, np.array([1.0, 0.0]), 0.0)  # fmt: skip
    width, n = 21.42, 1500
    h = width / n
    a = np.diag(fg.sigma_a + np.array([s[0, 1], s[1, 0]])) - np.array(
        [[0.0, s[1, 0]], [s[0, 1], 0.0]]
    )
    a = a - np.outer(fg.chi, fg.nu_sigma_f) / k
    current = np.array(current)
    # unknowns (cell, group); D phi'' by cell-centred differences
    m = np.zeros((2 * n, 2 * n))
    rhs = np.zeros(2 * n)
    for i in range(n):
        for g in range(2):
            row = 2 * i + g
            d = fg.D[g] / h**2
            m[row, 2 * i : 2 * i + 2] += a[g]
            if i > 0:
                m[row, row] += d
                m[row, row - 2] -= d
            if i < n - 1:
                m[row, row] += d
                m[row, row + 2] -= d
    rhs[0:2] = current / h
    phi = np.linalg.solve(m, rhs).reshape(n, 2)
    average = phi.mean(axis=0)
    surface = phi[0] + current * h / (2 * fg.D)  # half a cell back to x = 0
    got = homogeneous_surface_flux(fg, average, width, k)
    np.testing.assert_allclose(got, surface, rtol=1e-4)
    # the balance that makes the average equivalent to the current
    np.testing.assert_allclose(a @ average * width, current, rtol=1e-8)


@pytest.mark.parametrize(
    ("t", "rho", "mu"),  # IAPWS R12-08, Table 4 (micro Pa s)
    [(298.15, 998.0, 889.735100), (298.15, 1200.0, 1437.649467),
     (373.15, 1000.0, 307.883622), (433.15, 1000.0, 217.685358),
     (873.15, 600.0, 77.430195), (1173.15, 400.0, 64.154608)],
)  # fmt: skip
def test_viscosity_verification_values(t, rho, mu):
    assert water.viscosity(t, rho) * 1e6 == pytest.approx(mu, rel=1e-8)


@pytest.mark.parametrize(
    ("p", "rho", "drho", "cp", "cv", "mu", "lam"),  # IAPWS R15-11, Table 7, 620 K
    [(20e6, 613.227777, 5.20937820, 7.63433705, 3.03793441, 70.9051068, 481.485195),
     (50e6, 699.226043, 1.84869007, 5.32047725, 2.91692653, 84.1527945, 545.038940)],
)  # fmt: skip
def test_conductivity_industrial_verification_values(p, rho, drho, cp, cv, mu, lam):
    t = 620.0
    assert water.density(p, t) == pytest.approx(rho, rel=1e-8)
    assert water.drho_dp(p, t) * 1e6 == pytest.approx(drho, rel=1e-8)
    assert water.heat_capacity(p, t) / 1e3 == pytest.approx(cp, rel=1e-8)
    assert water.heat_capacity_v(p, t) / 1e3 == pytest.approx(cv, rel=1e-8)
    assert water.viscosity(t, rho) * 1e6 == pytest.approx(mu, rel=1e-7)
    assert water.conductivity(p, t) * 1e3 == pytest.approx(lam, rel=1e-8)


def test_film_heat_transfer():
    from fission.feedback import convection

    p, t, g, d = 15.5e6, 580.0, 3500.0, 0.0112
    re, pr = convection.reynolds_prandtl(p, t, g, d)
    k = water.conductivity(p, t)
    assert convection.film_coefficient(p, t, g, d) == pytest.approx(
        0.023 * re**0.8 * pr**0.4 * k / d, rel=1e-12
    )
    assert convection.presser(1.326) == pytest.approx(
        0.9217 + 0.1478 * 1.326 - 0.1130 * math.exp(-7 * 0.326), rel=1e-12
    )
    rp, ri, ro, q, h = 0.4096e-2, 0.418e-2, 0.475e-2, 2e4, 4e4
    bare, filmed = (
        fuel.rod(q, 580.0, rp, ri, ro),
        fuel.rod(q, 580.0, rp, ri, ro, film=h),
    )
    assert filmed.clad_outer - 580.0 == pytest.approx(
        q / (2 * math.pi * ro * h), rel=1e-12
    )
    assert filmed.centre > bare.centre


def test_isothermal_is_the_zero_power_limit(tmp_path):
    """At (almost) zero power the coupled solve settles with fuel and water
    at the inlet temperature: the isothermal k at that temperature."""
    from fission.feedback.core import isothermal, solve

    _table(tmp_path / "fuel.json", True)
    _table(tmp_path / "reflector.json", False)
    (tmp_path / "core.toml").write_text(
        CORE.format(fuel=tmp_path / "fuel.json", reflector=tmp_path / "reflector.json")
    )
    model = load_model(tmp_path / "core.toml")
    k_cold, k_inlet, k_hot = isothermal(model, (555.0, 565.0, 575.0), axial=20.0)
    assert solve(model, axial=20.0, power_level=1e-7).k == pytest.approx(
        k_inlet, rel=1e-6
    )
    # the synthetic fuel is overmoderated (its thermal absorption, water and
    # boron, scales with the density): the isothermal coefficient is positive
    assert k_cold < k_inlet < k_hot


def test_critical_boron_secant():
    """rho linear in boron: the secant lands on it in one step; a little
    curvature costs a few more."""
    from fission.feedback.core import critical_boron

    def linear(ppm):
        return 1 / (1 - (0.05 - 1e-4 * ppm))  # critical at 500 ppm

    b, k, n = critical_boron(linear)
    assert b == pytest.approx(500.0, abs=1e-9) and k == pytest.approx(1.0)
    assert n == 4  # two starting points, the exact step, the confirming one

    def curved(ppm):
        return 1 / (1 - (0.05 - 1e-4 * ppm + 2e-8 * ppm**2))

    b, k, _ = critical_boron(curved, tol=1e-6)
    assert 1 - 1 / k == pytest.approx(0.0, abs=1e-12)


def test_critical_boron_of_the_core(tmp_path):
    """The synthetic core made critical by its boron, at zero power."""
    from fission.feedback.core import critical_boron, solve

    _table(tmp_path / "fuel.json", True)
    _table(tmp_path / "reflector.json", False)
    (tmp_path / "core.toml").write_text(
        CORE.format(fuel=tmp_path / "fuel.json", reflector=tmp_path / "reflector.json")
    )
    model = load_model(tmp_path / "core.toml")

    def k_of(ppm):
        return solve(model, axial=20.0, power_level=1e-7, boron=ppm).k

    b, k, _ = critical_boron(k_of)
    assert k == pytest.approx(1.0, abs=1e-5)
    # more boron, lower k: the answer brackets as it should
    assert k_of(b + 50) < 1 < k_of(b - 50)


def test_zppt_cases_reproduce_the_keno_results():
    """The KENO eigenvalues of Table P5-6, differenced as VERA does, give
    the KENO column of Table P5-7 (a check of the transcribed cases)."""
    from fission.feedback.zppt import CASES, derived

    d = derived({n: c.keno for n, c in CASES.items()})
    table_p5_7 = {"A": 898, "B": 875, "C": 984, "D": 1386, "SA": 447, "SB": 1066,
                  "SC": 499, "SD": 499}  # fmt: skip
    for bank, worth in table_p5_7.items():
        assert d[bank] == pytest.approx(worth, abs=1.0)
    assert d["DBW"] == pytest.approx(-10.21, abs=0.01)
    assert d["ITC"] == pytest.approx(-3.18, abs=0.02)
    assert len(CASES) == 32


def test_zppt_runs_the_cases_on_a_core(tmp_path):
    """Each case is the isothermal core at its own boron and temperature."""
    from fission.feedback.core import isothermal
    from fission.feedback.zppt import CASES, run

    _table(tmp_path / "fuel.json", True)
    _table(tmp_path / "reflector.json", False)
    (tmp_path / "core.toml").write_text(
        CORE.format(fuel=tmp_path / "fuel.json", reflector=tmp_path / "reflector.json")
    )
    model = load_model(tmp_path / "core.toml")
    k = run(model, {n: CASES[n] for n in (2, 20)}, axial=20.0)
    assert k[2] == isothermal(model, [565.0], axial=20.0, boron=1291.0)[0]
    assert k[20] == isothermal(model, [560.0], axial=20.0, boron=1291.0)[0]
