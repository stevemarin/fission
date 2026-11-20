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

An assembly lattice may be surrounded by a water **gap** (half the space
between neighbouring assemblies). It adds no regions: a point in the gap
belongs to the outermost region of the nearest edge cell (same sector, by
angle from that cell's centre), whose area grows by its share of the gap.
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

    def areas(
        self, pitch: float, extents: tuple[float, float, float, float] | None = None
    ) -> list[float]:
        """Exact area (cm^2) of every local FSR.

        extents: distances (left, right, bottom, top) from the centre to the
        edges of the region the cell covers, if not the square of side
        `pitch` (edge cells of a lattice with a gap)."""
        n = self.sectors
        h = pitch / 2
        left, right, bottom, top = extents or (h, h, h, h)
        out = []
        inner = 0.0
        for r in self.radii:  # annuli: every sector of a circle has 1/n of its area
            out.extend([math.pi * (r * r - inner * inner) / n] * n)
            inner = r
        for k in range(
            n
        ):  # outside the last circle: wedge of the square minus wedge of the disc
            wedge = _wedge_of_rectangle(
                TWO_PI * k / n, TWO_PI * (k + 1) / n, left, right, bottom, top
            )
            out.append(wedge - math.pi * inner * inner / n)
        return out


def _wedge_of_rectangle(
    theta1: float, theta2: float, left: float, right: float, bottom: float, top: float
) -> float:
    """Area of the part of the rectangle -left <= x <= right, -bottom <= y <= top
    (which contains the origin) between the rays at two angles (0 <= theta1 <
    theta2 <= 2 pi)."""

    def edge(theta: float) -> tuple[float, float]:
        c, s = math.cos(theta), math.sin(theta)
        scale = math.inf
        if c > 0:
            scale = min(scale, right / c)
        elif c < 0:
            scale = min(scale, left / -c)
        if s > 0:
            scale = min(scale, top / s)
        elif s < 0:
            scale = min(scale, bottom / -s)
        return c * scale, s * scale

    corners = [(right, top), (-left, top), (-left, -bottom), (right, -bottom)]
    inside = sorted(
        (math.atan2(y, x) % TWO_PI, x, y)
        for x, y in corners
        if theta1 < math.atan2(y, x) % TWO_PI < theta2
    )
    points = [(0.0, 0.0), edge(theta1), *((x, y) for _, x, y in inside), edge(theta2)]
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
    gap: float = 0.0  # water around the lattice, on every side (cm)

    @property
    def width(self) -> float:
        return self.nx * self.pitch + 2 * self.gap

    @property
    def height(self) -> float:
        return self.ny * self.pitch + 2 * self.gap

    def centre(self, ix: int, iy: int) -> tuple[float, float]:
        return self.gap + (ix + 0.5) * self.pitch, self.gap + (iy + 0.5) * self.pitch

    def cell_widths(self) -> tuple[np.ndarray, np.ndarray]:
        """Widths of the cell columns and heights of the rows (edge cells
        include the gap)."""
        wx = np.full(self.nx, self.pitch)
        wy = np.full(self.ny, self.pitch)
        wx[[0, -1]] += self.gap
        wy[[0, -1]] += self.gap
        if self.nx == 1:
            wx[0] = self.pitch + 2 * self.gap
        if self.ny == 1:
            wy[0] = self.pitch + 2 * self.gap
        return wx, wy

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
        ix = min(max(int((x - self.gap) / self.pitch), 0), self.nx - 1)
        iy = min(max(int((y - self.gap) / self.pitch), 0), self.ny - 1)
        return ix, iy

    def region(self, x: float, y: float) -> int:
        """Global FSR of a point inside the lattice."""
        ix, iy = self.cell_at(x, y)
        cell = iy * self.nx + ix
        cx, cy = self.centre(ix, iy)
        dx, dy = x - cx, y - cy
        return self.offsets[cell] + self.cells[cell].region(dx, dy)

    def materials(self) -> list[str]:
        """Material of every FSR."""
        out = []
        for c in self.cells:
            for ring in range(len(c.radii) + 1):
                out.extend([c.materials[ring]] * c.sectors)
        return out

    def split_rings(self, material: str) -> Geometry:
        """The same mesh with the rings of `material` relabelled ``material/0``
        (innermost), ``material/1``... in every cell, so that each ring can
        get its own cross sections (spatial self-shielding)."""
        cells = []
        for c in self.cells:
            k, labels = 0, []
            for m in c.materials:
                if m == material:
                    labels.append(f"{m}/{k}")
                    k += 1
                else:
                    labels.append(m)
            cells.append(CellMesh(c.radii, tuple(labels), c.sectors))
        return Geometry(self.pitch, self.nx, self.ny, tuple(cells), self.gap)

    def areas(self) -> np.ndarray:
        """Exact area of every FSR (cm^2)."""
        h, g = self.pitch / 2, self.gap
        out = []
        for i, c in enumerate(self.cells):
            ix, iy = i % self.nx, i // self.nx
            extents = (
                h + (g if ix == 0 else 0.0), h + (g if ix == self.nx - 1 else 0.0),
                h + (g if iy == 0 else 0.0), h + (g if iy == self.ny - 1 else 0.0),
            )  # fmt: skip
            out.append(c.areas(self.pitch, extents))
        return np.concatenate(out)

    def cell_of_regions(self) -> np.ndarray:
        """Index of the cell that contains each FSR."""
        return np.repeat(np.arange(len(self.cells)), [c.regions for c in self.cells])

    @classmethod
    def from_core(
        cls,
        model: Model,
        rings: dict[str, int] | None = None,
        sectors: int = 8,
        z: float | None = None,
        outer_ring: float | None = None,
    ) -> Geometry:
        """FSR mesh of the whole core map, cut at height z (default: mid-height).

        rings, sectors, outer_ring: as in ``from_lattice``.

        Every assembly's lattice must be exactly one assembly pitch wide and
        all lattices must share one pin pitch, so that the core is one grid
        of cells. The core map is unfolded to the full map first.
        """
        core = model.core
        if core is None:
            raise ValueError(f"model {model.name!r} has no core")
        full = core.full_map()
        names = {a for row in full for a in row if a is not None}
        if any(a is None for row in full for a in row):
            raise ValueError(
                "core positions outside the core ('.') are not supported yet"
            )
        height = model.assemblies[next(iter(names))].height
        level = height / 2 if z is None else z
        slices: dict[str, Geometry] = {}
        for a in names:
            bottom = 0.0
            for seg in model.assemblies[a].segments:
                if bottom <= level < bottom + seg.height:
                    slices[a] = cls.from_lattice(
                        model, seg.lattice, rings, sectors, outer_ring
                    )
                bottom += seg.height
        if any(g.gap for g in slices.values()):
            raise ValueError("lattices with a gap are not supported in a core map yet")
        pitches = {g.pitch for g in slices.values()}
        if len(pitches) != 1:
            raise ValueError(f"lattices at z = {level:g} cm have different pin pitches")
        for a, g in slices.items():
            if abs(g.width - core.assembly_pitch) > 1e-9 * core.assembly_pitch:
                raise ValueError(f"assembly {a!r} is {g.width:g} cm wide, not the "
                                 f"assembly pitch {core.assembly_pitch:g} cm")  # fmt: skip
        n = slices[next(iter(names))].nx  # pins per assembly side
        rows, cols = len(full), len(full[0])
        cells = []
        for iy in range(rows * n):  # from the bottom; the map's rows run from the top
            a_row = full[rows - 1 - iy // n]
            for ix in range(cols * n):
                g = slices[a_row[ix // n]]  # type: ignore[index]
                cells.append(g.cells[(iy % n) * n + ix % n])
        return cls(pitches.pop(), cols * n, rows * n, tuple(cells))

    @classmethod
    def from_lattice(
        cls,
        model: Model,
        lattice: str,
        rings: dict[str, int] | None = None,
        sectors: int = 8,
        outer_ring: float | None = None,
    ) -> Geometry:
        """FSR mesh of a model lattice.

        rings: material -> number of equal-area rings per ring of that material
        (default 1); sectors: rays per cell (1 = no sectors); outer_ring: an
        extra circle (cm) in the material that fills the square, outside the
        pin, for cells whose rings end inside it. It changes only the mesh: the
        coolant near the pin and in the corners get separate flat sources.
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
            if outer_ring is not None and inner < outer_ring < lat.pitch / 2:
                radii.append(outer_ring)
                materials.append(pin.materials[-1])
            materials.append(pin.materials[-1])
            meshes[name] = CellMesh(tuple(radii), tuple(materials), sectors)
        n = len(lat.pins)
        # model maps list rows from the top; cells count from the bottom
        cells = tuple(
            meshes[lat.pins[n - 1 - iy][ix]] for iy in range(n) for ix in range(n)
        )
        return cls(lat.pitch, n, n, cells, lat.gap)
