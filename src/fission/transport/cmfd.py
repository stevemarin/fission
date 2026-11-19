"""Coarse-mesh finite-difference (CMFD) acceleration of MOC.

Unaccelerated MOC needs ~1000 sweeps on C5G7 because each sweep moves
neutrons only one track length. CMFD adds, after every sweep, a small
problem on a coarse mesh (here: one cell per pin cell) that couples the whole
core at once and is cheap to solve to convergence.

**Coarse balance.** Integrating the transport equation over coarse cell c
(volume V_c) gives, exactly,

    sum_faces J_out + Sigma_t,c V_c phi_c = V_c [ sum_g' Sigma_s,g'->g,c phi_c,g'
                                                  + chi_c / k sum_g' nu Sigma_f,c,g' phi_c,g' ]

with flux-volume-weighted cross sections (Sigma_x,c = sum_i A_i Sigma_x,i phi_i
/ sum_i A_i phi_i) and J the net currents through the faces, tallied by the
sweep (sum W psi over the tracks crossing a face, signed by direction).

**Closing it.** To solve for new fluxes, the currents must be expressed in
the coarse fluxes. Finite-difference diffusion gives, across the face
between cells L and R (with phi = 0 outside the core),

    J = -D~ (phi_R - phi_L) - D^ (phi_R + phi_L),

D~ = 2 D_L D_R / (D_L + D_R), D = 1 / (3 Sigma_t) (per unit face length and
cell width, which cancel on a uniform mesh). On a vacuum side D~ = 2 D h /
(h + 4 D) (the Marshak condition), on a reflective side D~ = D^ = 0. The
**correction** D^ is chosen so that the formula reproduces the MOC current
exactly at the current fluxes:

    D^ = -(J_MOC + D~ (phi_R - phi_L)) / (phi_R + phi_L).

The coarse problem therefore has the same solution as MOC once MOC has
converged; diffusion theory only shapes how fast it gets there.

**Solving it.** Power iteration on k; each outer iteration sweeps the groups
(Gauss-Seidel in energy, so downscatter is used at once) and each group's
5-point system gets a few red-black successive over-relaxation (SOR) passes.

**Stability.** With optically thick coarse cells, plain CMFD can diverge;
the correction is therefore under-relaxed from one iteration to the next
(RELAX).

**Prolongation.** The coarse solution rescales the MOC solution: every FSR
flux in cell c is multiplied by phi_c(new) / phi_c(MOC), and so are the
incoming angular fluxes of the tracks that start in c; k is taken from the
coarse problem.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .trace import face_counts

# Under-relaxation of the current correction D^ between iterations. Plain
# CMFD (RELAX = 1) diverges on C5G7: pin cells of water are ~3 mean free
# paths thick in the thermal groups, and coarse cells that thick make the
# correction overshoot. Mixing in the previous D^ damps this; 0.5-0.7
# converged C5G7 in ~20 iterations in our tests (docs ch. 15).
RELAX = 0.7


@dataclass
class Coarse:
    """Homogenized coarse-mesh data; arrays are (groups, ny, nx), row 0 at the bottom."""

    volume: np.ndarray  # (ny, nx)
    flux: np.ndarray
    total: np.ndarray
    nu_fission: np.ndarray
    chi: np.ndarray
    scatter: np.ndarray  # (g, g', ny, nx): from g to g'


def homogenize(problem, flux: np.ndarray) -> Coarse:
    """Flux-volume-weighted cross sections of every lattice cell."""
    g = problem.geometry
    ny, nx, groups = g.ny, g.nx, problem.groups
    cell = g.cell_of_regions()
    n = len(g.cells)
    areas = problem.areas

    def cell_sum(values: np.ndarray) -> np.ndarray:
        """Sum per cell over FSRs; values (regions, ...) -> (cells, ...)."""
        out = np.zeros((n, *values.shape[1:]))
        np.add.at(out, cell, values)
        return out

    volume = cell_sum(areas)
    af = areas[:, None] * flux  # (regions, groups)
    cf = cell_sum(af)  # (cells, groups)

    def weighted(name: str) -> np.ndarray:
        return cell_sum(af * problem.per_region(name)) / cf

    production = (problem.per_region("nu_fission") * flux).sum(axis=1)
    fission_source = cell_sum(areas * production)
    chi = cell_sum(areas[:, None] * production[:, None] * problem.per_region("chi"))
    has_fission = fission_source[:, None] > 0
    chi = np.divide(
        chi, fission_source[:, None], out=np.zeros_like(chi), where=has_fission
    )
    scatter = cell_sum(af[:, :, None] * problem.per_region("scatter")) / cf[:, :, None]

    def grid(a: np.ndarray) -> np.ndarray:  # (cells, groups) -> (groups, ny, nx)
        return np.ascontiguousarray(a.T.reshape(groups, ny, nx))

    return Coarse(
        volume.reshape(ny, nx),
        grid(cf / volume[:, None]),
        grid(weighted("total")),
        grid(weighted("nu_fission")),
        grid(chi),
        np.ascontiguousarray(
            np.transpose(scatter, (1, 2, 0)).reshape(groups, groups, ny, nx)
        ),
    )


def couplings(
    problem, coarse: Coarse, currents: np.ndarray, previous=None, relax: float = 1.0
):
    """(D~, D^) on x-faces (groups, ny, nx+1) and y-faces (groups, ny+1, nx).

    previous, relax: D^ <- relax D^ + (1 - relax) D^(previous) (see RELAX).
    """
    g = problem.geometry
    ny, nx, h = g.ny, g.nx, g.pitch
    groups = problem.groups
    n_xf, _ = face_counts(g)
    jx = currents[:n_xf].T.reshape(groups, ny, nx + 1)
    jy = currents[n_xf:].T.reshape(groups, ny + 1, nx)
    bc = problem.tracks.boundary
    d = 1.0 / (3.0 * coarse.total)
    phi = coarse.flux

    def vacuum(dc: np.ndarray) -> np.ndarray:
        return 2 * dc * h / (h + 4 * dc)

    dtx = np.zeros((groups, ny, nx + 1))
    dtx[:, :, 1:-1] = 2 * d[:, :, :-1] * d[:, :, 1:] / (d[:, :, :-1] + d[:, :, 1:])
    if bc["left"] == "vacuum":
        dtx[:, :, 0] = vacuum(d[:, :, 0])
    if bc["right"] == "vacuum":
        dtx[:, :, -1] = vacuum(d[:, :, -1])
    dty = np.zeros((groups, ny + 1, nx))
    dty[:, 1:-1, :] = 2 * d[:, :-1, :] * d[:, 1:, :] / (d[:, :-1, :] + d[:, 1:, :])
    if bc["bottom"] == "vacuum":
        dty[:, 0, :] = vacuum(d[:, 0, :])
    if bc["top"] == "vacuum":
        dty[:, -1, :] = vacuum(d[:, -1, :])

    zero_x = np.zeros((groups, ny, 1))
    left = np.concatenate([zero_x, phi], axis=2)  # phi_L of every x-face
    right = np.concatenate([phi, zero_x], axis=2)  # phi_R
    dhx = -(jx + dtx * (right - left)) / (right + left)
    zero_y = np.zeros((groups, 1, nx))
    below = np.concatenate([zero_y, phi], axis=1)
    above = np.concatenate([phi, zero_y], axis=1)
    dhy = -(jy + dty * (above - below)) / (above + below)
    # reflective sides carry no net current
    if bc["left"] == "reflective":
        dhx[:, :, 0] = 0.0
    if bc["right"] == "reflective":
        dhx[:, :, -1] = 0.0
    if bc["bottom"] == "reflective":
        dhy[:, 0, :] = 0.0
    if bc["top"] == "reflective":
        dhy[:, -1, :] = 0.0
    if previous is not None:
        dhx = relax * dhx + (1 - relax) * previous[0]
        dhy = relax * dhy + (1 - relax) * previous[1]
    return dtx, dhx, dty, dhy


@dataclass
class Operator:
    """The coarse 5-point operator per group: diag phi_c + sum a_n phi_n."""

    diag: np.ndarray  # (groups, ny, nx), includes removal (Sigma_t - Sigma_s,gg) V
    east: np.ndarray
    west: np.ndarray
    north: np.ndarray
    south: np.ndarray
    inscatter: np.ndarray  # (g', g, ny, nx): V Sigma_s,g'->g, zero on the diagonal
    fission: np.ndarray  # (groups, ny, nx): V nu Sigma_f
    chi: np.ndarray


def operator(coarse: Coarse, dtx, dhx, dty, dhy) -> Operator:
    """Cell balance: outflow through the faces + removal = in-scatter + fission."""
    groups = coarse.total.shape[0]
    v = coarse.volume[None]
    self_scatter = np.stack([coarse.scatter[g, g] for g in range(groups)])
    diag = (coarse.total - self_scatter) * v
    # J(east face) = -dt (phi_E - phi_c) - dh (phi_E + phi_c), outflow +J
    diag = diag + dtx[:, :, 1:] - dhx[:, :, 1:] + dtx[:, :, :-1] + dhx[:, :, :-1]
    diag = diag + dty[:, 1:, :] - dhy[:, 1:, :] + dty[:, :-1, :] + dhy[:, :-1, :]
    inscatter = coarse.scatter * v[None]
    for g in range(groups):
        inscatter[g, g] = 0.0
    return Operator(
        diag,
        -dtx[:, :, 1:] - dhx[:, :, 1:],
        -dtx[:, :, :-1] + dhx[:, :, :-1],
        -dty[:, 1:, :] - dhy[:, 1:, :],
        -dty[:, :-1, :] + dhy[:, :-1, :],
        inscatter,
        coarse.nu_fission * v,
        coarse.chi,
    )


def _neighbours(phi: np.ndarray, op: Operator, g: int) -> np.ndarray:
    """sum_n a_n phi_n for group g (phi = 0 outside)."""
    p = np.pad(phi, 1)
    return (
        op.east[g] * p[1:-1, 2:]
        + op.west[g] * p[1:-1, :-2]
        + op.north[g] * p[2:, 1:-1]
        + op.south[g] * p[:-2, 1:-1]
    )


def solve_coarse(
    op: Operator,
    phi: np.ndarray,
    k: float,
    tol: float = 1e-10,
    max_outer: int = 5000,
    inner: int = 4,
    omega: float = 1.5,
) -> tuple[float, np.ndarray, int]:
    """Power iteration on the coarse problem (Python reference).

    Returns (k, phi, outer iterations). phi keeps the total fission
    production of the starting flux.
    """
    groups, ny, nx = phi.shape
    phi = phi.copy()
    red = (np.add.outer(np.arange(ny), np.arange(nx)) % 2) == 0
    production = float((op.fission * phi).sum())
    target = production
    for outer in range(1, max_outer + 1):
        source_f = op.fission * phi
        fission = source_f.sum(axis=0) / k  # (ny, nx)
        old = phi.copy()
        for g in range(groups):
            b = op.chi[g] * fission + np.einsum("hyx,hyx->yx", op.inscatter[:, g], phi)
            for _ in range(inner):
                for colour in (red, ~red):
                    update = (b - _neighbours(phi[g], op, g)) / op.diag[g]
                    phi[g] = np.where(
                        colour, (1 - omega) * phi[g] + omega * update, phi[g]
                    )
        new_production = float((op.fission * phi).sum())
        k_new = k * new_production / production
        phi *= target / new_production
        production = target
        change = float(np.max(np.abs(phi - old)) / np.max(np.abs(phi)))
        converged = abs(k_new - k) < tol and change < 100 * tol
        k = k_new
        if converged:
            return k, phi, outer
    return k, phi, max_outer


def residual(op: Operator, phi: np.ndarray, k: float) -> float:
    """Largest relative imbalance of the coarse equations at (phi, k)."""
    groups = phi.shape[0]
    fission = (op.fission * phi).sum(axis=0) / k
    worst = 0.0
    for g in range(groups):
        lhs = op.diag[g] * phi[g] + _neighbours(phi[g], op, g)
        rhs = op.chi[g] * fission + np.einsum("hyx,hyx->yx", op.inscatter[:, g], phi)
        worst = max(worst, float(np.max(np.abs(lhs - rhs)) / np.max(np.abs(rhs))))
    return worst


def solve_coarse_rust(
    op: Operator,
    phi: np.ndarray,
    k: float,
    tol: float = 1e-10,
    max_outer: int = 5000,
    inner: int = 4,
    omega: float = 1.5,
) -> tuple[float, np.ndarray, int]:
    """Rust port of ``solve_coarse`` (same algorithm, same arguments)."""
    from fission import _rust

    groups, ny, nx = phi.shape

    def flat(a: np.ndarray) -> np.ndarray:
        return np.ascontiguousarray(a, dtype=float).ravel()

    k_new, out, outer = _rust.cmfd_solve(
        groups, ny, nx, flat(op.diag), flat(op.east), flat(op.west), flat(op.north),
        flat(op.south), flat(op.inscatter), flat(op.fission), flat(op.chi), flat(phi),
        k, tol, max_outer, inner, omega,
    )  # fmt: skip
    return k_new, out.reshape(phi.shape), outer
