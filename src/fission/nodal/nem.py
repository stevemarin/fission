"""The nodal expansion method (NEM) with coarse-mesh finite difference.

Finite difference assumes the flux linear inside each node, which needs
nodes of a few centimetres. NEM keeps one node per assembly (20 cm) and
describes the flux inside it better. Integrating the diffusion equation
over the two transverse directions gives, along x in a node of width h
(local coordinate xi in [-1/2, 1/2]),

    -(D / h^2) phi''(xi) + Sigma_r phi(xi) = Q(xi) - L(xi),

with Q the fission and in-scattering source and L the transverse leakage
(what leaks out through the y and z faces, per unit volume). The flux is
expanded in polynomials of zero mean,

    phi(xi) = phi_bar + a1 P1 + a2 P2 + a3 P3 + a4 P4,
    P1 = xi, P2 = 3 xi^2 - 1/4, P3 = xi (xi^2 - 1/4),
    P4 = (xi^2 - 1/20) (xi^2 - 1/4).

P1 and P2 have the values -+1/2 and 1/2 at the faces, P3 and P4 vanish
there, so a1 = phi(1/2) - phi(-1/2) and a2 = phi(1/2) + phi(-1/2) - 2
phi_bar follow from the face fluxes; a3 and a4 from the first two moments
of the equation (weighted with P1 and P2). The source Q is expanded in the
same polynomials (it is a combination of the group fluxes), and L is a
quadratic that keeps the average transverse leakage of the node and of its
two neighbours.

**Two-node problems.** Across every interface, the unknowns are, per group,
a1..a4 of both nodes; the equations are the four moment equations, the
continuity of the flux and of the current at the interface, and the
currents on the two far faces, taken from the coarse problem (they satisfy
every node's balance, which keeps the iteration consistent). The current
found sets D^ of the coarse finite-difference problem (``diffusion``) so
that it reproduces it:

    J = -D~ (phi_hi - phi_lo) - D^ (phi_hi + phi_lo),

and the coarse problem gives the new eigenvalue and node averages. On an
external face the same is done with one node and the boundary condition.
This alternation (Smith's non-linear iteration) converges to the NEM
solution. Discontinuity factors f enter as phi_heterogeneous = f
phi_homogeneous at each face (1 for the IAEA benchmark).
"""

from __future__ import annotations

import time
from dataclasses import dataclass, field
from typing import Any

import numpy as np
from numpy.polynomial import Polynomial as P

from . import diffusion
from .diffusion import AXES, Couplings
from .mesh import CoreMesh

# Under-relaxation of D^. Thin nodes next to thick ones (the 1.27 cm
# unpoisoned layer above the Pyrex of wbn1-pyrex, under 21.5 cm wide nodes)
# make the undamped iteration oscillate with period two and diverge. Half
# steps converge it in 16 iterations; IAEA 3D and wbn1-zoned take about as
# many as undamped (18 / 16 against 18 / 17).
RELAX = 0.5

BASIS = (
    P([0, 1]),
    P([-0.25, 0, 3]),
    P([0, -0.25, 0, 1]),
    P([0.0125, 0, -0.3, 0, 1]),  # (xi^2 - 1/20) (xi^2 - 1/4)
)


def _integral(p: P) -> float:
    q = p.integ()
    return float(q(0.5) - q(-0.5))


# M[m, i] = int P_m P_i, K[m, i] = int P_m P_i'' (m = moments 1, 2; i = 1..4)
M = np.array([[_integral(BASIS[m] * BASIS[i]) for i in range(4)] for m in range(2)])
K = np.array(
    [[_integral(BASIS[m] * BASIS[i].deriv(2)) for i in range(4)] for m in range(2)]
)
# derivatives P_i'(xi) at the faces xi = -1/2 (index 0) and +1/2 (index 1)
DP = np.array([[float(BASIS[i].deriv()(x)) for i in range(4)] for x in (-0.5, 0.5)])


@dataclass
class Nodes:
    """Node data on the full grid ([G, nz, ny, nx]; c is [G, G', ...])."""

    d: np.ndarray
    removal: np.ndarray
    c: np.ndarray  # source into g from phi_g': chi_g nu Sigma_f,g' / k + Sigma_s,g'->g
    phi: np.ndarray
    adf: np.ndarray


def _nodes(op: diffusion.Operator, mesh: CoreMesh, phi: np.ndarray, k: float) -> Nodes:
    groups = mesh.groups
    c = np.einsum("g...,h...->gh...", op.chi, op.nu_fission) / k
    for g in range(groups):
        for h in range(groups):
            if g != h:
                c[g, h] += op.scatter[h, g]
    adf = mesh.field(lambda x: x.adf if x.adf is not None else np.ones(x.groups))
    adf = np.where(mesh.inside[None], adf, 1.0)
    return Nodes(mesh.field(lambda x: x.diffusion_coefficient), op.removal, c, phi, adf)


def _all(*index: Any) -> tuple[Any, ...]:
    """Every group, at the given node or face indices: (slice(None), *index)."""
    return (slice(None), *index)


def _pair(x: np.ndarray, axis: int) -> tuple[np.ndarray, np.ndarray]:
    """Values of the nodes below and above every face along `axis`."""
    pad = [(0, 0)] * x.ndim
    pad[axis] = (1, 1)
    p = np.pad(x, pad)
    n = x.shape[axis]
    return np.take(p, range(n + 1), axis), np.take(p, range(1, n + 2), axis)


def transverse_leakage(mesh: CoreMesh, j: list[np.ndarray]) -> list[np.ndarray]:
    """Average leakage per unit volume through the faces normal to the two
    other axes, for each axis [G, nz, ny, nx]."""
    net = []
    for a, axis in enumerate(AXES):
        n = mesh.shape[axis - 1]
        h = mesh.widths[a]
        shape = [1, 1, 1, 1]
        shape[axis] = n
        diff = np.take(j[a], range(1, n + 1), axis) - np.take(j[a], range(n), axis)
        net.append(diff / h.reshape(shape))
    return [sum(net[b] for b in range(3) if b != a) for a in range(3)]


def leakage_shape(
    mesh: CoreMesh, lbar: np.ndarray, a: int
) -> tuple[np.ndarray, np.ndarray]:
    """(l1, l2): the quadratic L = Lbar + l1 P1 + l2 P2 of every node along
    axis a that keeps the averages of the node and its two neighbours (a
    missing neighbour - boundary or void - repeats the node's own value,
    exact at a symmetry plane)."""
    axis = AXES[a]
    h = mesh.widths[a]
    shape = [1, 1, 1, 1]
    shape[axis] = len(h)
    hn = np.broadcast_to(h.reshape(shape), lbar.shape)
    inside = np.broadcast_to(mesh.inside[None], lbar.shape)

    def neighbour(step: int) -> tuple[np.ndarray, np.ndarray]:
        value = np.roll(lbar, -step, axis)
        width = np.roll(hn, -step, axis)
        ok = np.roll(inside, -step, axis)
        edge = [slice(None)] * 4
        edge[axis] = slice(-1, None) if step > 0 else slice(0, 1)
        ok = ok.copy()
        ok[tuple(edge)] = False
        return np.where(ok, value, lbar), np.where(ok, width, hn)

    (lm, hm), (lp, hp) = neighbour(-1), neighbour(+1)

    def averages(lo: np.ndarray, hi: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """Averages of P1 = xi and P2 = 3 xi^2 - 1/4 over [lo, hi] (node units)."""
        p1 = (lo + hi) / 2
        p2 = (hi**3 - lo**3) / (hi - lo) - 0.25
        return p1, p2

    r = hm / hn
    s = hp / hn
    m1, m2 = averages(-0.5 - r, np.full_like(r, -0.5))
    p1, p2 = averages(np.full_like(s, 0.5), 0.5 + s)
    # [m1 m2; p1 p2] [l1; l2] = [lm - L; lp - L]
    det = m1 * p2 - m2 * p1
    l1 = ((lm - lbar) * p2 - m2 * (lp - lbar)) / det
    l2 = (m1 * (lp - lbar) - p1 * (lm - lbar)) / det
    return l1, l2


# values of P_i at the faces xi = -1/2 (index 0) and +1/2 (index 1)
PV = np.array([[float(BASIS[i](x)) for i in range(4)] for x in (-0.5, 0.5)])


def _moments(a: np.ndarray, r: np.ndarray, rows: int, cols: int, node: dict,
             g_count: int) -> None:  # fmt: skip
    """Write the two moment equations of each group of one node into the
    batched system a [B, n, n], r [B, n]: rows from `rows`, the node's
    coefficients a_g,i in columns cols + 4 g + i."""
    d, removal, c, h = node["d"], node["removal"], node["c"], node["h"]
    for g in range(g_count):
        for m in range(2):
            row = rows + 2 * g + m
            r[:, row] = -(node["l1"][:, g] * M[m, 0] + node["l2"][:, g] * M[m, 1])
            for gp in range(g_count):
                coef = -c[:, g, gp][:, None] * M[m][None, :]  # [B, 4]
                if gp == g:
                    coef = coef + (-(d[:, g] / h**2)[:, None] * K[m][None, :]
                                   + removal[:, g][:, None] * M[m][None, :])  # fmt: skip
                a[:, row, cols + 4 * gp : cols + 4 * gp + 4] = coef


def _current_row(a, r, row, col, node, g, face, value, sign=1.0):
    """sign * J(face) of group g = value, with J = -(D / h) phi'(face)."""
    w = -sign * node["d"][:, g] / node["h"]
    a[:, row, col : col + 4] += w[:, None] * DP[face][None, :]
    r[:, row] += value


def two_node(lo: dict, hi: dict, g_count: int) -> tuple[np.ndarray, np.ndarray]:
    """Solve the two-node problems of a batch of interfaces.

    lo, hi: node data [B, ...] (keys d, removal, c, h, adf, phi, l1, l2, and
    j_far, the current in the + direction on the node's other face, from the
    coarse problem). Unknowns per group: a1..a4 of both nodes. Equations:
    two moments per node, the far currents, continuity of the flux (with
    discontinuity factors) and of the current at the interface. Returns
    (heterogeneous face flux [B, G], current in the + direction [B, G])."""
    b = lo["d"].shape[0]
    n = 8 * g_count
    a = np.zeros((b, n, n))
    r = np.zeros((b, n))
    col_l, col_r = 0, 4 * g_count
    _moments(a, r, 0, col_l, lo, g_count)
    _moments(a, r, 2 * g_count, col_r, hi, g_count)
    for g in range(g_count):
        base = 4 * g_count + 4 * g
        cl, cr = col_l + 4 * g, col_r + 4 * g
        _current_row(a, r, base, cl, lo, g, 0, lo["j_far"][:, g])  # left face of lo
        _current_row(
            a, r, base + 1, cr, hi, g, 1, hi["j_far"][:, g]
        )  # right face of hi
        # interface current: J_lo(+1/2) - J_hi(-1/2) = 0
        _current_row(a, r, base + 2, cl, lo, g, 1, 0.0)
        _current_row(a, r, base + 2, cr, hi, g, 0, 0.0, sign=-1.0)
        # flux: f_lo (phi_lo + sum P_i(1/2) a_lo) = f_hi (phi_hi + sum P_i(-1/2) a_hi)
        f_l, f_h = lo["adf"][:, g], hi["adf"][:, g]
        a[:, base + 3, cl : cl + 4] = f_l[:, None] * PV[1][None, :]
        a[:, base + 3, cr : cr + 4] = -f_h[:, None] * PV[0][None, :]
        r[:, base + 3] = f_h * hi["phi"][:, g] - f_l * lo["phi"][:, g]
    x = np.linalg.solve(a, r[..., None])[..., 0]
    coef_l = x[:, : 4 * g_count].reshape(b, g_count, 4)
    phi_s = lo["adf"] * (lo["phi"] + coef_l @ PV[1])
    current = -(lo["d"] / lo["h"][:, None]) * (coef_l @ DP[1])
    return phi_s, current


def one_node(node: dict, sigma: float, beta: np.ndarray, g_count: int,
             ) -> tuple[np.ndarray, np.ndarray]:  # fmt: skip
    """A node with an external face at xi = sigma / 2, where the outgoing
    current is beta phi_b (beta = c f on vacuum, 0 on reflective faces); its
    other face carries j_far (+ direction) from the coarse problem.
    Returns (homogeneous face flux [B, G], outgoing current [B, G])."""
    b = node["d"].shape[0]
    n = 4 * g_count
    a = np.zeros((b, n, n))
    r = np.zeros((b, n))
    _moments(a, r, 0, 0, node, g_count)
    face, far = (1, 0) if sigma > 0 else (0, 1)
    for g in range(g_count):
        row = 2 * g_count + 2 * g
        col = 4 * g
        _current_row(a, r, row, col, node, g, far, node["j_far"][:, g])
        # J_out - beta phi_b = 0, J_out = sigma J(face), phi_b = phi + sum P_i a_i
        _current_row(a, r, row + 1, col, node, g, face, 0.0, sign=sigma)
        a[:, row + 1, col : col + 4] -= beta[:, g][:, None] * PV[face][None, :]
        r[:, row + 1] = beta[:, g] * node["phi"][:, g]
    x = np.linalg.solve(a, r[..., None])[..., 0].reshape(b, g_count, 4)
    phi_b = node["phi"] + x @ PV[face]
    j_out = sigma * (-(node["d"] / node["h"][:, None]) * (x @ DP[face]))
    return phi_b, j_out


def effective(base: Couplings, phi: np.ndarray) -> Couplings:
    """Couplings with non-negative coefficients that carry the same currents.

    Where |D^| > D~ an interior face would give the coarse operator a
    negative coefficient, and its iterations can diverge. There the face
    current J (at the flux phi) is written one-sidedly instead: for J > 0,
    D~ = -D^ = J / (2 phi_lo), so that J = 2 D~ phi_lo; for J < 0,
    D~ = D^ = -J / (2 phi_hi). The current is unchanged, so is the solution
    once the iteration has converged."""
    tilde, hat = [], []
    for a, axis in enumerate(AXES):
        t0, h0 = base.tilde[a], base.hat[a]
        lo, hi = _pair(phi, axis)
        kind = np.broadcast_to(base.kind[a], t0.shape)
        j = -t0 * (hi - lo) - h0 * (hi + lo)
        bad = (kind == 1) & (np.abs(h0) > t0)
        with np.errstate(divide="ignore", invalid="ignore"):
            up = j / (2 * lo)
            down = -j / (2 * hi)
        t = np.where(bad, np.where(j > 0, up, down), t0)
        h = np.where(bad, np.where(j > 0, -up, down), h0)
        tilde.append(np.nan_to_num(t))
        hat.append(np.nan_to_num(h))
    return Couplings(tilde, hat, base.kind)


def _gather(mesh: CoreMesh, nodes: Nodes, hn: np.ndarray, l1: np.ndarray, l2: np.ndarray,
            where: np.ndarray, j_far: np.ndarray) -> dict:  # fmt: skip
    """Node data at flat grid indices `where` (into [nz, ny, nx]), per batch entry."""
    flat = np.unravel_index(where, mesh.shape)

    def per_group(x):
        return x[_all(*flat)].T

    return {
        "d": per_group(nodes.d),
        "removal": per_group(nodes.removal),
        "c": np.moveaxis(nodes.c[_all(slice(None), *flat)], -1, 0),
        "h": hn[flat],
        "adf": per_group(nodes.adf),
        "phi": per_group(nodes.phi),
        "j_far": per_group(j_far),
        "l1": per_group(l1),
        "l2": per_group(l2),
    }


def _update(
    base: Couplings, a: int, faces: tuple, value: np.ndarray, relax: float
) -> None:
    """D^ of the faces <- (1 - relax) old + relax new (value is [B, G])."""
    old = base.hat[a][_all(*faces)]
    base.hat[a][_all(*faces)] = (1 - relax) * old + relax * value.T


@dataclass
class Result:
    k: float
    flux: np.ndarray  # node averages [G, nz, ny, nx]
    iterations: int
    history: list[float] = field(default_factory=list)
    couplings: Couplings | None = None
    seconds: float = 0.0
    hat: list[np.ndarray] | None = None  # D^ per axis (relative to the FD D~)


def solve(mesh: CoreMesh, tol_k: float = 1e-7, tol_flux: float = 1e-6,
          max_iterations: int = 500, relax: float = RELAX, omega: float = 1.0,
          verbose: bool = False, start: Result | None = None) -> Result:  # fmt: skip
    """NEM eigenvalue: alternate two-node problems and the corrected coarse
    finite-difference problem until the eigenvalue and the fluxes settle.
    relax < 1 under-relaxes the update of D^ (see RELAX). start: a solution of a nearby
    problem (same mesh, other constants) to begin from (flux, k and D^)."""
    t0 = time.perf_counter()
    base = diffusion.couplings(mesh)  # D~ of finite difference, D^ of NEM
    if start is not None and start.hat is not None:
        base.hat = [h.copy() for h in start.hat]
        c = effective(base, start.flux)
        sol = diffusion.solve(mesh, c, flux=start.flux, k=start.k, tol_k=tol_k / 10,
                              tol_source=tol_flux / 10, omega=omega)  # fmt: skip
    else:
        sol = diffusion.solve(mesh, base, tol_k=tol_k / 10, tol_source=tol_flux / 10)
        c = base
    groups = mesh.groups
    inside = mesh.inside
    history = [sol.k]
    names = (("west", "east"), ("north", "south"), ("bottom", "top"))
    index = np.arange(int(np.prod(mesh.shape))).reshape(mesh.shape)
    it = 0
    for it in range(1, max_iterations + 1):
        op = diffusion.operator(mesh, c)
        nodes = _nodes(op, mesh, sol.flux, sol.k)
        j = diffusion.currents(mesh, c, sol.flux)
        lbar = transverse_leakage(mesh, j)
        for a, axis in enumerate(AXES):
            l1, l2 = leakage_shape(mesh, lbar[a], a)
            h = mesh.widths[a]
            shape = [1, 1, 1]
            shape[axis - 1] = len(h)
            hn = np.broadcast_to(h.reshape(shape), mesh.shape)
            n = mesh.shape[axis - 1]
            j_minus = np.take(j[a], range(n), axis)  # current on each node's - face
            j_plus = np.take(j[a], range(1, n + 1), axis)

            kind = c.kind[a]
            idx_lo, idx_hi = _pair(index, axis - 1)

            f = np.nonzero(kind == 1)  # interior faces
            if len(f[0]):
                lo_i, hi_i = idx_lo[f], idx_hi[f]
                data = (mesh, nodes, hn, l1, l2)
                _, cur = two_node(_gather(*data, lo_i, j_minus), _gather(*data, hi_i, j_plus),
                                  groups)  # fmt: skip
                p_lo = sol.flux[_all(*np.unravel_index(lo_i, mesh.shape))].T
                p_hi = sol.flux[_all(*np.unravel_index(hi_i, mesh.shape))].T
                t = base.tilde[a][_all(*f)].T
                _update(base, a, f, -(cur + t * (p_hi - p_lo)) / (p_hi + p_lo), relax)
            # external faces: kind 2 (node below the face), 3 (node above)
            for k_value, sigma in ((2, 1.0), (3, -1.0)):
                f = np.nonzero(kind == k_value)
                if not len(f[0]):
                    continue
                node_i = (idx_lo if k_value == 2 else idx_hi)[f]
                pos = f[axis - 1]
                edge = names[a][1] if k_value == 2 else names[a][0]
                at_edge = pos == (n if k_value == 2 else 0)
                vacuum = ~at_edge | (mesh.boundary[edge] == "vacuum")
                node_n = np.unravel_index(node_i, mesh.shape)
                adf = nodes.adf[_all(*node_n)].T
                beta = np.where(vacuum[:, None], mesh.vacuum * adf, 0.0)
                far = j_minus if k_value == 2 else j_plus
                _, j_out = one_node(_gather(mesh, nodes, hn, l1, l2, node_i, far), sigma,
                                    beta, groups)  # fmt: skip
                p = sol.flux[_all(*node_n)].T
                t = base.tilde[a][_all(*f)].T
                _update(base, a, f, j_out / p - t, relax)
        old = sol
        c = effective(base, old.flux)
        sol = diffusion.solve(mesh, c, flux=old.flux, k=old.k, tol_k=tol_k / 10,
                              tol_source=tol_flux / 10, omega=omega)  # fmt: skip
        history.append(sol.k)
        change = np.max(np.abs(sol.flux - old.flux)[:, inside]) / np.max(
            sol.flux[:, inside]
        )
        if verbose:
            print(
                f"  NEM {it:3d}: k = {sol.k:.7f}, flux change {change:.1e}", flush=True
            )
        if abs(sol.k - old.k) < tol_k and change < tol_flux:
            break
    return Result(sol.k, sol.flux, it, history, c, time.perf_counter() - t0,
                  [h.copy() for h in base.hat])  # fmt: skip
