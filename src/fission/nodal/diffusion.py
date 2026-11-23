"""Coarse-mesh finite-difference diffusion: the operator and its eigenvalue.

Balance of node i in group g (volume V, faces f of area A):

    sum_f A J_f,out + Sigma_r,g V phi_g = V (chi_g / k sum_g' nu Sigma_f,g' phi_g'
                                             + sum_g' != g Sigma_s,g'->g phi_g')

with Sigma_r = Sigma_a + out-scattering. Through a face between nodes lo and
hi (hi further along the axis), the current in the + direction is

    J = -D~ (phi_hi - phi_lo) - D^ (phi_hi + phi_lo).

D~ = 2 D_lo D_hi / (D_lo h_hi + D_hi h_lo) is the finite-difference coupling
(the current of a flux that is linear in each node, with the current
continuous at the face). D^ is a correction that a better method (NEM) can
compute so that the coarse problem reproduces its currents; it is 0 for
plain finite difference. On an external face the outgoing current is
(D~_b + D^_b) phi: 0 on a reflective face, and with the Robin condition
d phi / dn = -(c / D) phi, eliminating the surface flux gives
D~_b = 2 D c / (2 D + c h).

Solved by power iteration; each group's diffusion equation by red-black
successive over-relaxation (SOR), the nodes coloured like a chess board in
3D so that every node of one colour depends only on the other colour.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .mesh import CoreMesh

AXES = (3, 2, 1)  # array axis of x, y, z in [G, nz, ny, nx]


@dataclass
class Couplings:
    """D~ and D^ on every face: per axis a, arrays [G, ...] with n_a + 1
    faces along that axis (face f between nodes f - 1 and f)."""

    tilde: list[np.ndarray]
    hat: list[np.ndarray]
    kind: list[np.ndarray]  # 0 none, 1 interior, 2 boundary with lo inside, 3 hi inside


def _faces(inside: np.ndarray, axis: int, low: str, high: str, boundary: dict[str, str],
           ) -> tuple[np.ndarray, np.ndarray]:  # fmt: skip
    """Face kinds along one axis of the [nz, ny, nx] grid, and which boundary
    faces are vacuum (void neighbours always are)."""
    pad = [(0, 0)] * 3
    pad[axis] = (1, 1)
    p = np.pad(inside, pad, constant_values=False)
    n = inside.shape[axis]
    lo = np.take(p, range(n + 1), axis=axis)
    hi = np.take(p, range(1, n + 2), axis=axis)
    kind = np.where(lo & hi, 1, np.where(lo & ~hi, 2, np.where(~lo & hi, 3, 0)))
    vacuum = np.ones(kind.shape, dtype=bool)
    first = [slice(None)] * 3
    first[axis] = slice(0, 1)
    last = [slice(None)] * 3
    last[axis] = slice(n, n + 1)
    vacuum[tuple(first)] = boundary[low] == "vacuum"
    vacuum[tuple(last)] = boundary[high] == "vacuum"
    return kind, vacuum


def couplings(mesh: CoreMesh) -> Couplings:
    """Finite-difference couplings D~ of the mesh (D^ = 0)."""
    d = mesh.field(lambda x: x.diffusion_coefficient)  # [G, nz, ny, nx]
    names = (("west", "east"), ("north", "south"), ("bottom", "top"))
    tilde, hat, kinds = [], [], []
    for a, axis in enumerate(AXES):
        kind, vacuum = _faces(mesh.inside, axis - 1, *names[a], mesh.boundary)
        h = mesh.widths[a]
        shape = [1, 1, 1, 1]
        shape[axis] = len(h)
        hh = np.broadcast_to(h.reshape(shape), d.shape)
        pad = [(0, 0)] * 4
        pad[axis] = (1, 1)
        dp = np.pad(d, pad)
        hp = np.pad(hh, pad, constant_values=1.0)
        n = d.shape[axis]
        d_lo, d_hi = np.take(dp, range(n + 1), axis), np.take(dp, range(1, n + 2), axis)
        h_lo, h_hi = np.take(hp, range(n + 1), axis), np.take(hp, range(1, n + 2), axis)
        k = np.broadcast_to(kind, d_lo.shape)
        v = np.broadcast_to(vacuum, d_lo.shape)
        c = mesh.vacuum
        with np.errstate(divide="ignore", invalid="ignore"):
            interior = 2 * d_lo * d_hi / (d_lo * h_hi + d_hi * h_lo)
            b_lo = 2 * d_lo * c / (2 * d_lo + c * h_lo)  # boundary, lo inside
            b_hi = 2 * d_hi * c / (2 * d_hi + c * h_hi)
        t = np.where(k == 1, interior, 0.0)
        t = np.where((k == 2) & v, b_lo, t)
        t = np.where((k == 3) & v, b_hi, t)
        tilde.append(np.nan_to_num(t))
        hat.append(np.zeros_like(t))
        kinds.append(kind)
    return Couplings(tilde, hat, kinds)


def _areas(mesh: CoreMesh) -> list[np.ndarray]:
    """Face areas per axis, broadcastable to the face arrays [G, nz, ny, n+1]."""
    dx, dy, dz = mesh.widths
    return [
        (dz[:, None] * dy[None, :])[None, :, :, None],
        (dz[:, None] * dx[None, :])[None, :, None, :],
        (dy[:, None] * dx[None, :])[None, None, :, :],
    ]


@dataclass
class Operator:
    """diag phi = b + sum over the six neighbours of coef * phi_neighbour."""

    diag: np.ndarray  # [G, nz, ny, nx]
    coef: list[
        tuple[np.ndarray, np.ndarray]
    ]  # per axis: (on the lower, on the upper neighbour)
    removal: np.ndarray
    scatter: np.ndarray  # [G', G, nz, ny, nx]: from g' to g, without self-scatter
    nu_fission: np.ndarray
    chi: np.ndarray
    volume: np.ndarray  # [nz, ny, nx]
    inside: np.ndarray


def operator(mesh: CoreMesh, c: Couplings) -> Operator:
    groups = mesh.groups
    scatter_full = np.zeros((groups, groups, *mesh.shape))
    inside = mesh.inside
    table = np.array([x.scatter for x in mesh.xs])  # [m, G, G']
    for g in range(groups):
        for h in range(groups):
            if g != h:
                scatter_full[g, h][inside] = table[mesh.material[inside], g, h]
    absorption = mesh.field(lambda x: x.absorption)
    removal = absorption + scatter_full.sum(axis=1)  # out-scattering
    volume = mesh.volumes()
    diag = removal * volume[None]
    coef = []
    for a, axis in enumerate(AXES):
        area = _areas(mesh)[a]
        t, hat, kind = c.tilde[a], c.hat[a], c.kind[a]
        n = mesh.shape[axis - 1]
        lower = np.arange(n)  # the face below node i is face i
        upper = lower + 1
        a_lo = np.take(area * np.ones_like(t), lower, axis)  # area of the face below
        a_hi = np.take(area * np.ones_like(t), upper, axis)
        t_lo, h_lo, k_lo = (
            np.take(x, lower, axis) for x in (t, hat, np.broadcast_to(kind, t.shape))
        )
        t_hi, h_hi, k_hi = (
            np.take(x, upper, axis) for x in (t, hat, np.broadcast_to(kind, t.shape))
        )
        # face below the node: the node is "hi" of that face
        diag += np.where((k_lo == 1) | (k_lo == 3), a_lo * (t_lo + h_lo), 0.0)
        c_lower = np.where(k_lo == 1, a_lo * (t_lo - h_lo), 0.0)
        # face above the node: the node is "lo" of that face
        diag += np.where(k_hi == 1, a_hi * (t_hi - h_hi), 0.0)
        diag += np.where(k_hi == 2, a_hi * (t_hi + h_hi), 0.0)
        c_upper = np.where(k_hi == 1, a_hi * (t_hi + h_hi), 0.0)
        coef.append((c_lower, c_upper))
    diag = np.where(inside[None], diag, 1.0)
    return Operator(
        diag, coef, removal, scatter_full, mesh.field(lambda x: x.nu_fission),
        mesh.field(lambda x: x.chi), volume, inside,
    )  # fmt: skip


def _sweep(diag, coef, b, phi, colour, omega) -> None:
    """One SOR half-sweep over the nodes of one colour (in place)."""
    s = b.copy()
    for a in range(3):
        axis = 2 - a  # x is the last axis of [nz, ny, nx]
        c_lower, c_upper = coef[a]
        lo = np.zeros_like(phi)
        hi = np.zeros_like(phi)
        sl = [slice(None)] * 3
        src = [slice(None)] * 3
        sl[axis], src[axis] = slice(1, None), slice(0, -1)
        lo[tuple(sl)] = phi[tuple(src)]
        sl[axis], src[axis] = slice(0, -1), slice(1, None)
        hi[tuple(sl)] = phi[tuple(src)]
        s += c_lower * lo + c_upper * hi
    new = s / diag
    phi[colour] += omega * (new[colour] - phi[colour])


def _rust_sweeps(diag, coef, inside, b, phi, count, omega) -> np.ndarray:
    """`count` red-black SOR sweeps of one group in Rust (port of `_sweep`)."""
    from fission import _rust

    shape = phi.shape
    nz, ny, nx = shape
    flat = [np.ascontiguousarray(x, dtype=float).ravel() for pair in coef for x in pair]
    out = _rust.diffusion_sweeps(
        nz, ny, nx, np.ascontiguousarray(diag).ravel(), flat,
        np.ascontiguousarray(inside, dtype=np.uint8).ravel(),
        np.ascontiguousarray(b).ravel(), np.ascontiguousarray(phi).ravel(), count, omega,
    )  # fmt: skip
    return out.reshape(shape)


@dataclass
class Solution:
    k: float
    flux: np.ndarray  # [G, nz, ny, nx], normalized to sum nu Sigma_f phi V = 1
    outer: int
    history: list[float] = field(default_factory=list)


def solve(mesh: CoreMesh, c: Couplings | None = None, flux: np.ndarray | None = None,
          k: float = 1.0, tol_k: float = 1e-7, tol_source: float = 1e-6,
          inner: int = 4, omega: float = 1.6, max_outer: int = 20000,
          use_rust: bool = True) -> Solution:  # fmt: skip
    """Eigenvalue and flux by power iteration (red-black SOR inners).
    use_rust: the inner sweeps in Rust (``_rust.diffusion_sweeps``, the same
    arithmetic, bit for bit); False runs the Python reference ``_sweep``."""
    c = couplings(mesh) if c is None else c
    op = operator(mesh, c)
    groups = mesh.groups
    nz, ny, nx = mesh.shape
    zz, yy, xx = np.meshgrid(range(nz), range(ny), range(nx), indexing="ij")
    colours = [((zz + yy + xx) % 2 == p) & op.inside for p in (0, 1)]
    if flux is None:
        phi = np.repeat(np.where(op.inside, 1.0, 0.0)[None], groups, axis=0)
    else:
        phi = flux.copy()

    def production(f: np.ndarray) -> np.ndarray:
        return (op.nu_fission * f).sum(axis=0) * op.volume

    source = production(phi)
    total = source.sum()
    phi /= total
    source /= total
    history = []
    for outer in range(1, max_outer + 1):
        for g in range(groups):
            b = op.chi[g] * source / k
            for h in range(groups):
                if h != g:
                    b = b + op.scatter[h, g] * phi[h] * op.volume
            b = np.where(op.inside, b, 0.0)
            coef = [(lower[g], upper[g]) for lower, upper in op.coef]
            if use_rust:
                phi[g] = _rust_sweeps(
                    op.diag[g], coef, op.inside, b, phi[g], inner, omega
                )
            else:
                for _ in range(inner):
                    for colour in colours:
                        _sweep(op.diag[g], coef, b, phi[g], colour, omega)
        new = production(phi)
        k_new = k * new.sum() / source.sum()
        change = np.max(np.abs(new / new.sum() - source / source.sum())) / np.max(
            source / source.sum()
        )
        history.append(k_new)
        dk = abs(k_new - k)
        k = k_new
        norm = new.sum()
        phi /= norm
        source = new / norm
        if dk < tol_k and change < tol_source:
            break
    return Solution(k, phi, outer, history)


def currents(mesh: CoreMesh, c: Couplings, phi: np.ndarray) -> list[np.ndarray]:
    """Net current in the + direction on every face, per axis [G, ...]."""
    out = []
    for a, axis in enumerate(AXES):
        pad = [(0, 0)] * 4
        pad[axis] = (1, 1)
        p = np.pad(phi, pad)
        n = phi.shape[axis]
        lo = np.take(p, range(n + 1), axis)
        hi = np.take(p, range(1, n + 2), axis)
        t, h, kind = c.tilde[a], c.hat[a], np.broadcast_to(c.kind[a], c.tilde[a].shape)
        j = np.where(kind == 1, -t * (hi - lo) - h * (hi + lo), 0.0)
        j = np.where(kind == 2, (t + h) * lo, j)
        j = np.where(kind == 3, -(t + h) * hi, j)
        out.append(j)
    return out
