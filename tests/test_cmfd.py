"""Face currents and CMFD acceleration (step B3)."""

import numpy as np
import pytest

from fission.model import MODELS_DIR, parse
from fission.transport import Geometry, cmfd, cyclic_tracks
from fission.transport.moc import FOUR_PI, _rust_sweep, build, solve, sweep
from fission.transport.trace import face_counts

MIXED = {
    "left": "reflective",
    "bottom": "reflective",
    "right": "vacuum",
    "top": "vacuum",
}


def _lattice(boundary="reflective", azim=16, spacing=0.05, gap=0.0):
    """3 x 3 pins: UO2 and MOX around a water cell (in a water gap)."""
    model = parse(
        {
            "model": {"format": 1, "name": "3x3", "multigroup": "c5g7-xs.toml"},
            "materials": {
                "uo2": {"macroscopic": "UO2"},
                "mox": {"macroscopic": "MOX87"},
                "water": {"macroscopic": "moderator"},
            },
            "pins": {
                "u": {"radii": [0.54], "materials": ["uo2", "water"]},
                "m": {"radii": [0.54], "materials": ["mox", "water"]},
                "w": {"radii": [], "materials": ["water"]},
            },
            "lattices": {
                "c": {
                    "pitch": 1.26,
                    "legend": {"U": "u", "M": "m", "W": "w"},
                    "map": "U M U\nM W M\nU M U",
                    "gap": gap,
                }
            },
        },
        MODELS_DIR,
    )
    g = Geometry.from_lattice(model, "c", rings={"uo2": 2, "mox": 2}, sectors=4)
    tracks = cyclic_tracks(g.width, g.height, azim, spacing, boundary)
    return build(model, g, tracks, n_polar=2)


def _random_sweep_inputs(problem, seed=1):
    rng = np.random.default_rng(seed)
    sigma_t = problem.per_region("total")
    source = rng.uniform(0.1, 1.0, sigma_t.shape)
    psi = rng.uniform(
        0, 1, (len(problem.tracks), 2, len(problem.polar.sin), problem.groups)
    )
    psi = sweep(problem, sigma_t, source, psi)[1]  # zero where entering from vacuum
    return sigma_t, source, psi


@pytest.mark.parametrize("boundary", ["reflective", MIXED])
def test_currents_python_rust_and_balance(boundary):
    p = _lattice(boundary)
    sigma_t, source, psi = _random_sweep_inputs(p)
    tally, _, leakage, cur = sweep(p, sigma_t, source, psi)
    rt, _, rl, rc = _rust_sweep(p, sigma_t, source, psi)
    assert cur is not None and rc is not None
    np.testing.assert_allclose(rc, cur, rtol=1e-12, atol=1e-12 * np.abs(cur).max())
    np.testing.assert_allclose(rt, tally, rtol=1e-12, atol=1e-12 * np.abs(tally).max())
    assert rl == pytest.approx(leakage, rel=1e-12)
    g = p.geometry
    nx, ny = g.nx, g.ny
    n_xf, _ = face_counts(g)
    jx = cur[:n_xf].reshape(ny, nx + 1, -1)
    jy = cur[n_xf:].reshape(ny + 1, nx, -1)
    # per cell: sum of W delta psi over its regions = -(net outflow)
    per_cell = np.zeros((len(g.cells), p.groups))
    np.add.at(per_cell, g.cell_of_regions(), tally)
    outflow = (jx[:, 1:] - jx[:, :-1] + jy[1:, :] - jy[:-1, :]).reshape(-1, p.groups)
    np.testing.assert_allclose(per_cell, -outflow, atol=1e-11 * np.abs(outflow).max())
    # leakage = outward current through the vacuum sides
    vacuum = 0.0
    if p.tracks.boundary["right"] == "vacuum":
        vacuum += jx[:, -1].sum()
    if p.tracks.boundary["top"] == "vacuum":
        vacuum += jy[-1, :].sum()
    assert vacuum == pytest.approx(leakage, rel=1e-12, abs=1e-12)


def test_coarse_equations_reproduce_converged_moc():
    """At the MOC solution, the corrected coarse problem is satisfied exactly."""
    p = _lattice(MIXED)
    sol = solve(p, cmfd=False, tol_k=1e-11, tol_source=1e-10)
    sigma_t = p.per_region("total")
    source = (np.einsum("rgh,rg->rh", p.per_region("scatter"), sol.flux)
              + p.per_region("chi") * ((p.per_region("nu_fission") * sol.flux).sum(1)
                                       / sol.k)[:, None]) / FOUR_PI  # fmt: skip
    # one more sweep with converged boundary fluxes gives the currents
    psi = np.zeros((len(p.tracks), 2, len(p.polar.sin), p.groups))
    for _ in range(300):
        tally, psi, _, cur = _rust_sweep(p, sigma_t, source, psi)
    flux = FOUR_PI * source / sigma_t + tally / (sigma_t * p.areas[:, None])
    coarse = cmfd.homogenize(p, flux)
    op = cmfd.operator(coarse, *cmfd.couplings(p, coarse, cur))
    assert cmfd.residual(op, coarse.flux, sol.k) < 1e-8
    k, _, _ = cmfd.solve_coarse(op, coarse.flux, sol.k, tol=1e-12)
    assert k == pytest.approx(sol.k, rel=1e-8)


def test_coarse_solver_rust_matches_python():
    p = _lattice(MIXED)
    sigma_t, source, psi = _random_sweep_inputs(p, seed=4)
    tally, psi, _, cur = _rust_sweep(p, sigma_t, source, psi)
    flux = FOUR_PI * source / sigma_t + tally / (sigma_t * p.areas[:, None])
    coarse = cmfd.homogenize(p, flux)
    op = cmfd.operator(coarse, *cmfd.couplings(p, coarse, cur))
    a = cmfd.solve_coarse(op, coarse.flux, 1.0, max_outer=40)
    b = cmfd.solve_coarse_rust(op, coarse.flux, 1.0, max_outer=40)
    assert b[0] == pytest.approx(a[0], rel=1e-13)
    np.testing.assert_allclose(b[1], a[1], rtol=1e-12)
    assert a[2] == b[2]


@pytest.mark.parametrize(("boundary", "gap"), [("reflective", 0.0), (MIXED, 0.0),
                                             ("reflective", 0.1), (MIXED, 0.1)])  # fmt: skip
def test_cmfd_gives_the_same_answer_much_faster(boundary, gap):
    p = _lattice(boundary, gap=gap)
    plain = solve(p, cmfd=False, tol_k=1e-9, tol_source=1e-8)
    fast = solve(p, cmfd=True, tol_k=1e-9, tol_source=1e-8)
    assert fast.k == pytest.approx(plain.k, rel=1e-8)

    # eigenvectors have an arbitrary scale: compare at equal total production
    def shape(sol):
        return sol.flux / sol.rates("nu_fission").sum()

    np.testing.assert_allclose(shape(fast), shape(plain), rtol=1e-6)
    assert fast.iterations * 3 < plain.iterations  # 144 -> 36 on the leaky lattice
    assert fast.balance() == pytest.approx(fast.k, rel=1e-8)


def test_surface_flux_of_a_uniform_medium():
    """One material everywhere, reflective: the flux is flat, so the surface
    flux equals the volume flux (discontinuity factors of 1)."""
    from fission.transport.moc import surface_flux

    p = _lattice()
    uniform = {m: p.xs[0] for m in set(p.geometry.materials())}
    from fission.transport.moc import with_materials

    p = with_materials(p, uniform)
    sol = solve(p, tol_k=1e-10, tol_source=1e-9)
    average = (p.areas[:, None] * sol.flux).sum(axis=0) / p.areas.sum()
    for side, phi in surface_flux(p, sol.psi).items():
        np.testing.assert_allclose(phi, average, rtol=2e-3, err_msg=side)


def test_face_flux_of_a_uniform_medium():
    """Flat flux: a line across the domain sees the volume flux. (Single
    faces do not: at grazing angles one crossing stands for a whole face,
    so how a line's crossings split between its faces is noise, here
    -0.57/+1.15/-0.57% over the three rows.)"""
    from fission.transport.moc import face_flux, with_materials

    p = _lattice()
    p = with_materials(p, {m: p.xs[0] for m in set(p.geometry.materials())})
    sol = solve(p, tol_k=1e-10, tol_source=1e-9)
    average = (p.areas[:, None] * sol.flux).sum(axis=0) / p.areas.sum()
    lines = [[1, 5, 9], [2, 6, 10], [15, 16, 17], [18, 19, 20]]  # x = p, 2p; y = p, 2p
    for line in lines:
        np.testing.assert_allclose(
            face_flux(sol, line).mean(axis=0), average, rtol=1e-6
        )


def test_face_flux_on_the_boundary_is_the_surface_flux():
    """On a reflective side the face tally (both directions) and
    surface_flux (entering only, times 2) count the same flux."""
    from fission.transport.moc import face_flux, surface_flux

    p = _lattice()
    sol = solve(p, tol_k=1e-10, tol_source=1e-9)
    left = face_flux(sol, [0, 4, 8]).mean(axis=0)  # x = 0 in the three rows
    bottom = face_flux(sol, [12, 13, 14]).mean(axis=0)  # y = 0
    side = surface_flux(p, sol.psi)
    np.testing.assert_allclose(left, side["left"], rtol=1e-6)
    np.testing.assert_allclose(bottom, side["bottom"], rtol=1e-6)
