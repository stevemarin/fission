"""Collision probabilities in a cylindrical cell (step C4: the CP cross-check).

An independent way to solve the pin cell: instead of following tracks
through the square cell (MOC), replace the square by a circle of the same
area (the **Wigner-Seitz cell**, outer radius R = p / sqrt(pi)) with a
**white** boundary (every neutron that leaves comes back, isotropically,
through the surface), and compute once per group the first-flight
probabilities P_ij that a neutron born uniformly and isotropically in ring i
makes its next collision in ring j. Then, with flat sources q_i,

    Sigma_j phi_j V_j = sum_i q_i V_i P_ij.

**Chords.** In an infinite cylinder only the height y of a chord matters.
Along a chord at height y the rings are crossed in the order N, N-1, ...,
k, ..., N. For a neutron moving in 3D at polar angle theta, a 2D distance s
is an optical path tau / sin(theta); integrating the polar angle turns the
exponential into the Bickley-Naylor functions

    Ki_n(x) = Int_0^{pi/2} cos^{n-1}(u) exp(-x / cos u) du,   dKi_n/dx = -Ki_{n-1}

[Bell & Glasstone sec. 2.8 introduce them for the slab and cylinder escape
problems]. Integrating Ki_1 over two segments a (region i) and b (region j,
downstream, optical thickness between them tau_ab) gives, per unit chord
width,

    N_ab = [Ki3(tau_ab) - Ki3(tau_ab + tau_a) - Ki3(tau_ab + tau_b)
            + Ki3(tau_ab + tau_a + tau_b)] / Sigma_a,

and a segment with itself N_aa = [tau_a - Ki3(0) + Ki3(tau_a)] / Sigma_a;
what passes the end of the chord escapes, [Ki3(tau_end) - Ki3(tau_end +
tau_a)] / Sigma_a. Summed over the chords these are V_i P_ij and V_i P_i,esc.
The weight: averaging the source over all directions of travel, every
direction gives the same chords in a cylinder, so one direction suffices,
and the chords at -y repeat those at y: 2 Int_0^R dy. Along one direction
only downstream pairs count, and since the ring sequence N...k...N is a
palindrome they include both i -> j and j -> i. Each source point's
collisions and escape add up to Ki2(0) = 1, so the totals are the chord
lengths: the volumes. Gauss-Legendre points in every annulus band [r_(m-1), r_m]
integrate over y; using the volumes of the same quadrature makes sum_j P_ij
+ P_i,esc = 1 exactly.

**White boundary.** A neutron entering through the surface collides first
in ring j with probability P_sj = 4 V_j Sigma_j P_j,esc / S (reciprocity,
S = 2 pi R) or crosses the cell, P_ss = 1 - sum_j P_sj. Returning every
escaping neutron as often as needed,

    P^w_ij = P_ij + P_i,esc P_sj / (1 - P_ss).
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

# Gauss-Legendre in the angle u of Ki_n (the integrand is smooth on [0, pi/2])
_U, _W = np.polynomial.legendre.leggauss(96)
_U = (_U + 1.0) * math.pi / 4.0
_W = _W * math.pi / 4.0


def bickley(n: int, x: np.ndarray) -> np.ndarray:
    """Ki_n(x) for x >= 0 (any shape)."""
    x = np.asarray(x, dtype=float)
    c = np.cos(_U)
    terms = c ** (n - 1) * np.exp(-x[..., None] / c)
    return terms @ _W


@dataclass
class CylinderCP:
    """First-flight probabilities of an annular cell, every group."""

    radii: np.ndarray  # (N,) outer radii, the last one the cell boundary
    volumes: np.ndarray  # (N,) quadrature areas (cm^2)
    sigma: np.ndarray  # (G, N) total cross sections
    collide: np.ndarray  # (G, N, N): P_ij
    escape: np.ndarray  # (G, N): P_i,esc

    def white(self) -> np.ndarray:
        """P^w_ij with a white outer boundary, shape (G, N, N)."""
        surface = 2.0 * math.pi * self.radii[-1]
        p_sj = 4.0 * self.volumes * self.sigma * self.escape / surface  # (G, N)
        p_ss = 1.0 - p_sj.sum(axis=1)
        return (
            self.collide
            + self.escape[:, :, None] * p_sj[:, None, :] / (1.0 - p_ss)[:, None, None]
        )


def chords(radii: np.ndarray, points: int = 24) -> tuple[np.ndarray, np.ndarray]:
    """Gauss-Legendre heights y in (0, R) and weights, band by band."""
    x, w = np.polynomial.legendre.leggauss(points)
    ys, ws = [], []
    inner = 0.0
    for r in radii:
        # y = inner + (r - inner) t with t = 1 - s^2 clusters points near the
        # band's outer edge, where the chord lengths vary like sqrt(r - y)
        s = 0.5 * (x + 1.0)
        t = 1.0 - s * s
        ys.append(inner + (r - inner) * t)
        ws.append((r - inner) * 2.0 * s * 0.5 * w)
        inner = r
    return np.concatenate(ys), np.concatenate(ws)


def collision_probabilities(
    radii: np.ndarray, sigma: np.ndarray, points: int = 24
) -> CylinderCP:
    """P_ij of an annular cell; sigma has shape (G, N) (or (N,))."""
    radii = np.asarray(radii, dtype=float)
    sigma = np.atleast_2d(np.asarray(sigma, dtype=float))
    groups, n = sigma.shape
    ys, ws = chords(radii, points)
    vol_collide = np.zeros((groups, n, n))
    vol_escape = np.zeros((groups, n))
    volumes = np.zeros(n)
    for y, w in zip(ys, ws, strict=True):
        # half-chord lengths in every ring the chord crosses
        half = np.sqrt(np.maximum(radii**2 - y * y, 0.0))
        inner = np.concatenate([[0.0], half[:-1]])
        k = int(np.searchsorted(radii, y, side="right"))  # innermost ring crossed
        lengths = half[k:] - inner[k:]  # rings k..N-1, one half each
        regions = list(range(n - 1, k, -1)) + [k] + list(range(k + 1, n))
        seg = np.concatenate([lengths[::-1][:-1], [2 * lengths[0]], lengths[1:]])
        tau = sigma[:, regions] * seg[None, :]  # (G, segments)
        for a, i in enumerate(regions):
            volumes[i] += 2.0 * w * seg[a]
        cum = np.concatenate([np.zeros((groups, 1)), np.cumsum(tau, axis=1)], axis=1)
        m = len(regions)
        ki0 = math.pi / 4.0
        for a in range(m):
            i = regions[a]
            ta = tau[:, a]
            sa = sigma[:, i]
            vol_collide[:, i, i] += 2.0 * w * (ta - ki0 + bickley(3, ta)) / sa
            for b in range(a + 1, m):
                j = regions[b]
                tb = tau[:, b]
                gap = cum[:, b] - cum[:, a + 1]
                combo = (
                    bickley(3, gap) - bickley(3, gap + ta)
                    - bickley(3, gap + tb) + bickley(3, gap + ta + tb)
                )  # fmt: skip
                # downstream pairs only: the chord's ring sequence is a
                # palindrome, so the pairs from j to i are among them too
                vol_collide[:, i, j] += 2.0 * w * combo / sa
            # past the end of the chord: escape
            ahead = cum[:, m] - cum[:, a + 1]
            vol_escape[:, i] += (
                2.0 * w * (bickley(3, ahead) - bickley(3, ahead + ta)) / sa
            )
    return CylinderCP(
        radii, volumes, sigma,
        vol_collide / volumes[None, :, None], vol_escape / volumes[None, :],
    )  # fmt: skip


def eigenvalue(
    cp: CylinderCP, scatter: np.ndarray, nu_fission: np.ndarray, chi: np.ndarray
) -> tuple[float, np.ndarray]:
    """k and the flux (G, N) of the cell with a white boundary.

    scatter (N, G, G') from g to g', nu_fission and chi (N, G). Dense:
    Sigma_g phi_g,j V_j = sum_i V_i P^w_g,ij q_g,i is assembled for all
    groups and rings at once and solved as A phi = F phi / k."""
    p = cp.white()
    groups, n, _ = p.shape
    v = cp.volumes
    # T[g]: phi_g = T[g] q_g, T_ji = V_i P_ij / (Sigma_j V_j)
    t = np.einsum("i,gij->gji", v, p) / (cp.sigma[:, :, None] * v[None, :, None])
    size = groups * n
    a = np.eye(size)
    f = np.zeros((size, size))
    for g in range(groups):
        rows = slice(g * n, (g + 1) * n)
        for h in range(groups):
            # in-scatter h -> g of every ring i, then transport i -> j
            s = scatter[:, h, g]  # (N,)
            if np.any(s):
                a[rows, h * n : (h + 1) * n] -= t[g] * s[None, :]
            nf = nu_fission[:, h]
            if np.any(nf):
                f[rows, h * n : (h + 1) * n] += t[g] * (chi[:, g] * nf)[None, :]
    m = np.linalg.solve(a, f)
    values, vectors = np.linalg.eig(m)
    top = int(np.argmax(values.real))
    phi = np.abs(vectors[:, top].real).reshape(groups, n)
    return float(values[top].real), phi


def pin_cell(problem) -> tuple[float, np.ndarray, list[str]]:
    """k of a single-cell MOC problem solved by CP on its Wigner-Seitz cell.

    Uses the problem's own cross sections (``moc.Problem.xs``, transport
    corrected) and rings; sectors do not matter in a 1D cell. Returns k, the
    ring fluxes (G, N) and the ring materials."""
    g = problem.geometry
    if len(g.cells) != 1:
        raise ValueError("the CP cross-check is for a single cell")
    cell = g.cells[0]
    radius = math.sqrt(g.width * g.height / math.pi)  # same area as the square
    radii = np.array([*cell.radii, radius])
    names = sorted(set(g.materials()))
    xs = [problem.xs[names.index(m)] for m in cell.materials]
    cp = collision_probabilities(radii, np.array([x.total for x in xs]).T)
    k, phi = eigenvalue(
        cp,
        np.array([x.scatter for x in xs]),
        np.array([x.nu_fission for x in xs]),
        np.array([x.chi for x in xs]),
    )
    return k, phi, list(cell.materials)
