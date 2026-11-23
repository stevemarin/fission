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
"""

from __future__ import annotations

import itertools
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from fission.model import MacroXS, Model

FACES = ("west", "east", "north", "south", "bottom", "top")  # -x +x -y +y -z +z
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


def core_mesh(
    model: Model, refine: int = 1, axial_size: float | None = None
) -> CoreMesh:
    """Nodes of the model's core (macroscopic, homogeneous lattices only).

    refine: nodes per assembly width (a half assembly on a symmetry plane
    gets max(1, refine // 2)); axial_size: largest axial node (default: the
    radial node size)."""
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
    used = {a for row in grid for a in row if a is not None}
    if not used:
        raise ValueError("the core map is empty")
    breaks = sorted({0.0} | {
        z for a in used
        for z in np.cumsum([s.height for s in model.assemblies[a].segments])
    })  # fmt: skip
    size = axial_size or pitch / refine
    dz: list[float] = []
    for lo, hi in itertools.pairwise(breaks):
        dz += _split(hi - lo, max(1, math.ceil((hi - lo) / size - 1e-9)))
    z_mid = np.cumsum(dz) - np.array(dz) / 2

    names = sorted({
        model.lattices[s.lattice].material or "" for a in used
        for s in model.assemblies[a].segments
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
            a = grid[r][c]
            if a is None:
                continue
            tops = np.cumsum([s.height for s in model.assemblies[a].segments])
            for k, z in enumerate(z_mid):
                seg = model.assemblies[a].segments[int(np.searchsorted(tops, z))]
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
