"""The 2D geometry seen by the transport solver: flat-source regions.

The method of characteristics assumes the source is constant in each of a
set of small regions, the **flat-source regions** (FSRs). A lattice is a
rectangular grid of square cells of one pitch. Every cell is a pin
(concentric rings, the last material filling the square) cut into:

* **rings**: each material ring can be split into n rings of equal area
  (radii r_k^2 = r_in^2 + k/n (r_out^2 - r_in^2)), because the flux varies
  most with radius (it dips inside the fuel at resonance energies);
* **sectors**: the whole cell is cut by n_s rays from its centre at angles
  2 pi k / n_s, so the water in the corners of the cell, far from the pin,
  gets its own regions.

FSR numbering: cell (ix, iy), counted from the bottom-left with ix fastest,
owns regions offset[cell] + ring * sectors + sector, ring 0 at the centre
and sector 0 starting at angle 0 (the +x direction).
"""

from __future__ import annotations

import bisect
import math
from dataclasses import dataclass

import numpy as np

from fission.model import Model

TWO_PI = 2.0 * math.pi


@dataclass(frozen=True)
class CellMesh:
    """Flat-source regions of one square cell: rings x sectors."""

    radii: tuple[float, ...]  # ring outer radii (cm), increasing
    materials: tuple[str, ...]  # one per ring, len(radii) + 1 (last fills the square)
    sectors: int = 1

    @property
    def regions(self) -> int:
        return (len(self.radii) + 1) * self.sectors

    def region(self, dx: float, dy: float) -> int:
        """Local FSR of the point at (dx, dy) from the cell centre."""
        r2 = dx * dx + dy * dy
        ring = bisect.bisect_right([r * r for r in self.radii], r2)
        if self.sectors == 1:
            return ring
        theta = math.atan2(dy, dx)
        if theta < 0.0:
            theta += TWO_PI
        sector = min(int(theta / (TWO_PI / self.sectors)), self.sectors - 1)
        return ring * self.sectors + sector

    def areas(self, pitch: float) -> list[float]:
        """Exact area (cm^2) of every local FSR."""
        n = self.sectors
        out = []
        inner = 0.0
        for r in self.radii:  # annuli: every sector of a circle has 1/n of its area
            out.extend([math.pi * (r * r - inner * inner) / n] * n)
            inner = r
        for k in range(
            n
        ):  # outside the last circle: wedge of the square minus wedge of the disc
            wedge = _wedge_of_square(TWO_PI * k / n, TWO_PI * (k + 1) / n, pitch / 2)
            out.append(wedge - math.pi * inner * inner / n)
        return out


def _wedge_of_square(theta1: float, theta2: float, half: float) -> float:
    """Area of the part of the square |x|, |y| <= half between two angles."""

    def edge(theta: float) -> tuple[float, float]:
        c, s = math.cos(theta), math.sin(theta)
        scale = half / max(abs(c), abs(s))
        return c * scale, s * scale

    points = [(0.0, 0.0), edge(theta1)]
    for k in range(8):  # the corners, at 45 + 90 k degrees
        corner = math.pi / 4 + k * math.pi / 2
        if theta1 < corner < theta2:
            points.append(edge(corner))
    points.append(edge(theta2))
    area = 0.0
    for (x1, y1), (x2, y2) in zip(points, points[1:] + points[:1], strict=True):
        area += x1 * y2 - x2 * y1  # shoelace formula
    return area / 2


def split_ring(r_in: float, r_out: float, n: int) -> list[float]:
    """Outer radii of n rings of equal area between r_in and r_out."""
    return [
        math.sqrt(r_in * r_in + (k / n) * (r_out * r_out - r_in * r_in))
        for k in range(1, n + 1)
    ]


@dataclass(frozen=True)
class Geometry:
    """A rectangular lattice of square cells, bottom-left corner at (0, 0)."""

    pitch: float
    nx: int
    ny: int
    cells: tuple[CellMesh, ...]  # index iy * nx + ix

    @property
    def width(self) -> float:
        return self.nx * self.pitch

    @property
    def height(self) -> float:
        return self.ny * self.pitch

    @property
    def offsets(self) -> list[int]:
        out, n = [], 0
        for c in self.cells:
            out.append(n)
            n += c.regions
        return out

    @property
    def regions(self) -> int:
        return sum(c.regions for c in self.cells)

    def cell_at(self, x: float, y: float) -> tuple[int, int]:
        ix = min(max(int(x / self.pitch), 0), self.nx - 1)
        iy = min(max(int(y / self.pitch), 0), self.ny - 1)
        return ix, iy

    def region(self, x: float, y: float) -> int:
        """Global FSR of a point inside the lattice."""
        ix, iy = self.cell_at(x, y)
        cell = iy * self.nx + ix
        dx = x - (ix + 0.5) * self.pitch
        dy = y - (iy + 0.5) * self.pitch
        return self.offsets[cell] + self.cells[cell].region(dx, dy)

    def materials(self) -> list[str]:
        """Material of every FSR."""
        out = []
        for c in self.cells:
            for ring in range(len(c.radii) + 1):
                out.extend([c.materials[ring]] * c.sectors)
        return out

    def areas(self) -> np.ndarray:
        """Exact area of every FSR (cm^2)."""
        return np.concatenate([c.areas(self.pitch) for c in self.cells])

    @classmethod
    def from_lattice(
        cls,
        model: Model,
        lattice: str,
        rings: dict[str, int] | None = None,
        sectors: int = 8,
    ) -> Geometry:
        """FSR mesh of a model lattice.

        rings: material -> number of equal-area rings per ring of that material
        (default 1); sectors: rays per cell (1 = no sectors).
        """
        lat = model.lattices[lattice]
        if lat.pins is None or lat.pitch is None:
            raise ValueError(
                f"lattice {lattice!r} is homogeneous: it has no pitch to trace"
            )
        if sectors < 1:
            raise ValueError("sectors must be >= 1")
        meshes: dict[str, CellMesh] = {}
        for name in {p for row in lat.pins for p in row}:
            pin = model.pins[name]
            radii, materials, inner = [], [], 0.0
            for r, m in zip(pin.radii, pin.materials, strict=False):
                n = (rings or {}).get(m, 1)
                radii.extend(split_ring(inner, r, n))
                materials.extend([m] * n)
                inner = r
            materials.append(pin.materials[-1])
            meshes[name] = CellMesh(tuple(radii), tuple(materials), sectors)
        n = len(lat.pins)
        # model maps list rows from the top; cells count from the bottom
        cells = tuple(
            meshes[lat.pins[n - 1 - iy][ix]] for iy in range(n) for ix in range(n)
        )
        return cls(lat.pitch, n, n, cells)
