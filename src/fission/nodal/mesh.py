"""The core of a model as a 3D grid of homogeneous nodes.

Axes: x along the map's columns, y along its rows (top to bottom), z up the
assembly stacks; arrays are indexed [z, y, x]. Every assembly position of
the map is split into `refine` x `refine` nodes; axially, every interval
between the segment boundaries of the stacks is split into nodes of at most
`axial_size` cm. Positions outside the map (".") are void: the faces
towards them are external boundaries.

With quarter symmetry the map's first row and column lie on the symmetry
planes, through the middle of the centre assembly: those positions are half
an assembly wide, and the two symmetry faces are reflective.

Control rods: `banks` gives the steps withdrawn of some banks of the
model's [core.banks]. At each position of such a bank the rod's absorber
zones cut the stack, and inside a zone a lattice is replaced by its
``rodded`` lattice for that absorber (a lattice without one, like the
reflector above the fuel, stays as it is). The zone edges become axial
node boundaries everywhere.
"""

from __future__ import annotations

import itertools
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from fission.model import AxialSegment, MacroXS, Model

FACES = ("west", "east", "north", "south", "bottom", "top")  # -x +x -y +y -z +z
TOUCH = 1e-6  # cm: axial edges closer than this are one edge
TABLES = Path("endf-data/tables")


@dataclass
class CoreMesh:
    widths: tuple[np.ndarray, np.ndarray, np.ndarray]  # dx [nx], dy [ny], dz [nz]
    material: np.ndarray  # [nz, ny, nx] index into xs, -1 outside the core
    xs: list[MacroXS]
    names: list[str]
    boundary: dict[str, str]  # face -> "vacuum" | "reflective"
    vacuum: float  # c in d phi / dn = -(c / D) phi on vacuum faces
    position: np.ndarray  # [ny, nx, 2]: (map row, map column) of each node column
    symmetry: str  # "full" | "quarter"
    base: np.ndarray | None = (
        None  # [nz, ny, nx] model material of each node (index into names)
    )
    tables: dict[int, object] | None = None  # names index -> XSTable of table materials

    @property
    def shape(self) -> tuple[int, int, int]:
        return self.material.shape  # type: ignore[return-value]

    @property
    def groups(self) -> int:
        return self.xs[0].groups

    @property
    def inside(self) -> np.ndarray:
        return self.material >= 0

    def volumes(self) -> np.ndarray:
        dx, dy, dz = self.widths
        return dz[:, None, None] * dy[None, :, None] * dx[None, None, :]

    def field(self, f) -> np.ndarray:
        """[G, nz, ny, nx] array of f(xs) (a per-group vector) per node, 0 outside."""
        table = np.array([f(x) for x in self.xs])  # [materials, G]
        out = np.zeros((self.groups, *self.shape))
        inside = self.inside
        out[:, inside] = table[self.material[inside]].T
        return out


def _split(width: float, count: int) -> list[float]:
    return [width / count] * count


def rodded(
    segments: tuple[AxialSegment, ...], model: Model, bottom: float,
    zones: tuple[tuple[str, float], ...],
) -> tuple[AxialSegment, ...]:  # fmt: skip
    """The stack with a rod's absorber column inserted from `bottom` (cm):
    segments cut at the zone edges, each piece inside a zone given its
    lattice's rodded lattice for that absorber (if it has one)."""
    edges = [bottom]
    for _, height in zones:
        edges.append(edges[-1] + height)
    out: list[AxialSegment] = []
    lo = 0.0
    for seg in segments:
        hi = lo + seg.height
        cuts = sorted({lo, hi} | {e for e in edges if lo + TOUCH < e < hi - TOUCH})
        for a, b in itertools.pairwise(cuts):
            mid = 0.5 * (a + b)
            lattice = seg.lattice
            for k, (absorber, _) in enumerate(zones):
                if edges[k] <= mid < edges[k + 1]:
                    lattice = model.lattices[seg.lattice].rodded.get(absorber, lattice)
            if out and out[-1].lattice == lattice:
                out[-1] = AxialSegment(lattice, out[-1].height + (b - a))
            else:
                out.append(AxialSegment(lattice, b - a))
        lo = hi
    return tuple(out)


def stacks(
    model: Model, banks: dict[str, float] | None = None
) -> dict[tuple[int, int], tuple[AxialSegment, ...]]:
    """Axial stack of every map position, with the given banks inserted
    (bank -> steps withdrawn)."""
    core = model.core
    assert core is not None
    banks = banks or {}
    if banks and core.banks is None:
        raise ValueError(f"{model.name} has no control rod banks ([core.banks])")
    rod = model.control_rods[core.banks.rod] if core.banks is not None else None
    for bank, steps in banks.items():
        assert core.banks is not None and rod is not None
        if bank not in core.banks.positions:
            raise ValueError(
                f"no bank {bank!r} (have {', '.join(core.banks.positions)})"
            )
        if not 0 <= steps <= rod.steps:
            raise ValueError(f"bank {bank}: {steps} steps, outside 0..{rod.steps}")
    out = {}
    for r, row in enumerate(core.map):
        for c, a in enumerate(row):
            if a is None:
                continue
            segments = model.assemblies[a].segments
            here = core.banks.at(r, c) if core.banks is not None else None
            if here is not None and here in banks and rod is not None:
                bottom = rod.bottom + banks[here] * rod.step
                segments = rodded(segments, model, bottom, rod.zones)
            out[r, c] = segments
    return out


def core_mesh(
    model: Model, refine: int = 1, axial_size: float | None = None,
    banks: dict[str, float] | None = None,
) -> CoreMesh:  # fmt: skip
    """Nodes of the model's core (macroscopic, homogeneous lattices only).

    refine: nodes per assembly width (a half assembly on a symmetry plane
    gets max(1, refine // 2)); axial_size: largest axial node (default: the
    radial node size); banks: control rod banks inserted (bank -> steps
    withdrawn), the others are out."""
    core = model.core
    if core is None:
        raise ValueError(f"{model.name} has no core")
    grid = core.map
    pitch = core.assembly_pitch
    quarter = core.symmetry == "quarter"
    rows, cols = len(grid), len(grid[0])

    def widths(count: int) -> tuple[list[float], list[int]]:
        out, owner = [], []
        for i in range(count):
            half = quarter and i == 0
            n = max(1, refine // 2) if half else refine
            w = _split(pitch / 2 if half else pitch, n)
            out += w
            owner += [i] * n
        return out, owner

    dx, col_of = widths(cols)
    dy, row_of = widths(rows)
    stack = stacks(model, banks)
    if not stack:
        raise ValueError("the core map is empty")
    breaks: list[float] = []
    for z in sorted({0.0} | {
        float(z) for segs in stack.values() for z in np.cumsum([s.height for s in segs])
    }):  # fmt: skip
        if (
            not breaks or z - breaks[-1] > TOUCH
        ):  # edges that coincide (rod top = fuel top)
            breaks.append(z)
    size = axial_size or pitch / refine
    dz: list[float] = []
    for lo, hi in itertools.pairwise(breaks):
        dz += _split(hi - lo, max(1, math.ceil((hi - lo) / size - 1e-9)))
    z_mid = np.cumsum(dz) - np.array(dz) / 2

    names = sorted({
        model.lattices[s.lattice].material or "" for segs in stack.values() for s in segs
    })  # fmt: skip
    # a table material may be listed by several lattices; nodes keep the name
    if "" in names:
        raise ValueError("pin lattices need homogenized constants (tables, M4b)")
    xs: list[MacroXS] = []
    tables: dict[int, object] = {}
    for n in names:
        mat = model.materials[n]
        if mat.table is not None:
            from fission.tables import XSTable

            if core.state is None:
                raise ValueError(
                    "table materials need a core state ([core] state = ...)"
                )
            path = (
                Path(mat.table)
                if mat.table.endswith(".json")
                else TABLES / f"{mat.table}.json"
            )
            table = XSTable.read(path)
            tables[len(xs)] = table
            xs.append(table.macro(model.states[core.state], n))
            continue
        if mat.macroscopic is None:
            raise ValueError(f"material {n!r} has no macroscopic cross sections")
        xs.append(model.multigroup[mat.macroscopic])
    index = {n: i for i, n in enumerate(names)}

    material = np.full((len(dz), len(dy), len(dx)), -1, dtype=np.int64)
    position = np.zeros((len(dy), len(dx), 2), dtype=np.int64)
    for j, r in enumerate(row_of):
        for i, c in enumerate(col_of):
            position[j, i] = (r, c)
            if grid[r][c] is None:
                continue
            segs = stack[r, c]
            tops = np.cumsum([s.height for s in segs])
            for k, z in enumerate(z_mid):
                seg = segs[min(int(np.searchsorted(tops, z)), len(segs) - 1)]
                material[k, j, i] = index[model.lattices[seg.lattice].material or ""]

    b = core.boundary
    boundary = {
        "west": "reflective" if quarter else b["left"],
        "east": b["right"],
        "north": "reflective" if quarter else b["top"],
        "south": b["bottom"],
        "bottom": "vacuum",
        "top": "vacuum",
    }
    return CoreMesh(
        (np.array(dx), np.array(dy), np.array(dz)), material, xs, names, boundary,
        core.vacuum_coefficient, position, core.symmetry, material.copy(), tables,
    )  # fmt: skip


def with_states(mesh: CoreMesh, states: dict[str, np.ndarray]) -> CoreMesh:
    """The mesh with the table materials evaluated at every node's own
    state (states: BRANCH_AXES name -> [nz, ny, nx]); other nodes keep their
    material's constants. Every table node gets its own MacroXS."""
    from dataclasses import replace

    assert mesh.base is not None and mesh.tables is not None
    xs = list(mesh.xs[: len(mesh.names)])
    material = mesh.base.copy()
    for k, j, i in zip(*np.nonzero(mesh.base >= 0), strict=True):
        m = int(mesh.base[k, j, i])
        table = mesh.tables.get(m)
        if table is None:
            continue
        state = {a: float(v[k, j, i]) for a, v in states.items()}
        material[k, j, i] = len(xs)
        xs.append(table.macro(state, mesh.names[m]))  # type: ignore[attr-defined]
    return replace(mesh, xs=xs, material=material)


def parse_banks(items: list[str] | None) -> dict[str, float]:
    """Command-line bank positions, "D=0" (steps withdrawn) or "D" (fully in)."""
    out = {}
    for item in items or []:
        bank, _, steps = item.partition("=")
        out[bank] = float(steps) if steps else 0.0
    return out
