"""Cartesian core geometry: node widths, material map, faces and boundary conditions."""

from dataclasses import dataclass

import numpy as np

from .materials import Material

SIDES = ("x-", "x+", "y-", "y+", "z-", "z+")
_BC_KEYWORDS = {"reflective": 0.0, "vacuum": 0.5, "zero_flux": np.inf}


def parse_albedo(bc) -> float:
    """Boundary condition as an albedo alpha with J_out = alpha * phi_surface.

    'reflective' -> 0, 'vacuum' (Marshak, zero incoming partial current) -> 0.5,
    'zero_flux' -> inf, or any non-negative number (e.g. 0.4692 for IAEA).
    """
    if isinstance(bc, str):
        return _BC_KEYWORDS[bc]
    return float(bc)


def refine(widths, target):
    """Split each coarse width into equal nodes of about `target` size.

    Returns (fine_widths, owner) where owner[i] is the coarse index of fine node i.
    """
    fine, owner = [], []
    for i, w in enumerate(widths):
        n = max(1, round(w / target))
        fine += [w / n] * n
        owner += [i] * n
    return np.array(fine), np.array(owner)


@dataclass
class Interfaces:
    """Faces shared by two active nodes along one direction (L at -, R at +)."""
    L: np.ndarray
    R: np.ndarray


@dataclass
class Boundaries:
    """Faces with an active node on one side only.

    side = +1 if the face is the node's + face, -1 if it is the - face.
    The condition is cJ * J_out = cphi * phi_surface, equivalently the incoming
    partial current is refl * (outgoing partial current).
    """
    node: np.ndarray
    side: np.ndarray
    cJ: np.ndarray
    cphi: np.ndarray

    @property
    def refl(self) -> np.ndarray:
        return (self.cJ - 2 * self.cphi) / (self.cJ + 2 * self.cphi)


class Core:
    """Structured 3D Cartesian node mesh.

    widths:    [dx, dy, dz] node widths (cm). Use dz=[1.0] for a 2D problem.
    matmap:    int array (nx, ny, nz) of material indices; -1 marks nodes
               outside the core (lets you model a staircase core boundary).
    materials: list of Material.
    bc:        dict keyed by 'x-','x+','y-','y+','z-','z+'. Faces next to a
               -1 node use the condition of the domain side they face.
    """

    def __init__(self, widths, matmap, materials: list[Material], bc: dict):
        self.widths = [np.asarray(w, float) for w in widths]
        matmap = np.asarray(matmap, int)
        self.shape = tuple(len(w) for w in self.widths)
        if matmap.shape != self.shape:
            raise ValueError(f"matmap shape {matmap.shape} != mesh {self.shape}")
        self.matmap = matmap
        self.materials = materials
        self.G = materials[0].G
        unknown = set(bc) - set(SIDES)
        if unknown:
            raise ValueError(f"unknown boundary sides {unknown}")
        self.albedo = {s: parse_albedo(bc.get(s, "reflective")) for s in SIDES}

        active = matmap >= 0
        self.N = int(active.sum())
        self.idx = np.full(self.shape, -1)
        self.idx[active] = np.arange(self.N)
        self.ijk = np.argwhere(active)  # C order, matches idx numbering
        self.h = np.stack([self.widths[d][self.ijk[:, d]] for d in range(3)], axis=1)
        self.vol = self.h.prod(axis=1)
        self.mat = matmap[active]

        def table(attr):
            return np.array([getattr(m, attr) for m in materials])[self.mat]

        self.D = table("D")
        self.sigma_r = table("sigma_r")
        self.nu_sigma_f = table("nu_sigma_f")
        self.kappa_sigma_f = table("kappa_sigma_f")
        self.chi = table("chi")
        self.sigma_s = table("sigma_s")
        self.adf = table("adf")
        self._build_faces()

    def _build_faces(self) -> None:
        self.interfaces: list[Interfaces] = []
        self.boundaries: list[Boundaries] = []
        self.nbr = np.full((self.N, 3, 2), -1)  # [node, dir, 0=minus/1=plus]
        self.reflective = np.zeros((self.N, 3, 2), bool)  # face is a reflective boundary
        self.active_dirs = []
        for d in range(3):
            pad = [(0, 0)] * 3
            pad[d] = (1, 1)
            p = np.pad(self.idx, pad, constant_values=-1)
            lo_sl, hi_sl = [slice(None)] * 3, [slice(None)] * 3
            lo_sl[d], hi_sl[d] = slice(0, -1), slice(1, None)
            lo, hi = p[tuple(lo_sl)].ravel(), p[tuple(hi_sl)].ravel()

            both = (lo >= 0) & (hi >= 0)
            L, R = lo[both], hi[both]
            self.interfaces.append(Interfaces(L, R))
            self.nbr[L, d, 1] = R
            self.nbr[R, d, 0] = L

            bp, bm = (lo >= 0) & (hi < 0), (lo < 0) & (hi >= 0)
            node = np.concatenate([lo[bp], hi[bm]])
            side = np.concatenate([np.ones(bp.sum(), int), -np.ones(bm.sum(), int)])
            a_plus, a_minus = self.albedo[SIDES[2 * d + 1]], self.albedo[SIDES[2 * d]]
            alpha = np.where(side > 0, a_plus, a_minus)
            inf = np.isinf(alpha)
            cJ = np.where(inf, 0.0, 1.0)
            cphi = np.where(inf, 1.0, alpha)
            self.boundaries.append(Boundaries(node, side, cJ, cphi))
            self.reflective[node, d, (side > 0).astype(int)] = cphi == 0

            # A direction with no interfaces and only reflective faces carries no
            # leakage (e.g. z in a 2D problem); skip it in the nodal solve.
            if L.size or np.any(cphi > 0):
                self.active_dirs.append(d)

    def to_grid(self, values):
        """Scatter per-node values (N, ...) onto the (nx, ny, nz, ...) mesh; NaN outside."""
        values = np.asarray(values)
        out = np.full(self.shape + values.shape[1:], np.nan)
        out[tuple(self.ijk.T)] = values
        return out
