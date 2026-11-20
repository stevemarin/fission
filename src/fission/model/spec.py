"""The reactor model as plain dataclasses: what a model file describes.

The layers each refer only to the layer before them:

    materials -> pins -> lattices (2D, xy) -> assemblies (stacks in z) -> core (xy map)

The 3D core is built from 2D slices: a lattice is one xy cross-section,
assumed infinitely tall (what a 2D lattice calculation sees), and an
assembly stacks lattices in z. Beside the layers sit the **states**, the
operating conditions (temperatures, water density, soluble boron) at which
the materials are evaluated. A material stores composition only; its
temperature and, for the coolant, its density come from the state.

`load.py` reads a TOML file into these classes and checks every reference.
"""

from __future__ import annotations

import itertools
import math
from dataclasses import dataclass, field
from typing import Literal

from .multigroup import MacroXS

SIDES = ("left", "right", "bottom", "top")

type Basis = Literal["atom", "weight"]


@dataclass(frozen=True)
class Isotopics:
    """Isotopic make-up of one element: percent by atom or by weight."""

    basis: Basis
    percent: dict[str, float]  # nuclide -> percent (sums to 100)


@dataclass(frozen=True)
class Material:
    """A composition, given in one of three ways.

    * macroscopic: ``macroscopic`` names an entry of the model's multigroup
      library (group constants given directly, as in C5G7);
    * explicit: ``atom_densities`` in atoms per barn-cm;
    * a compound: mass ``density`` (g/cm3, or "moderator" for the state's
      water density), a chemical ``formula`` (element -> atoms per molecule)
      and the ``isotopics`` of every element; optionally ``soluble_boron``,
      dissolved at the state's ppm by weight.
    """

    name: str
    temperature: Literal["fuel", "moderator"]  # which state temperature applies
    atom_densities: dict[str, float] | None = None
    density: float | Literal["moderator"] | None = None
    formula: dict[str, float] | None = None
    isotopics: dict[str, Isotopics] | None = None
    soluble_boron: Isotopics | None = None
    thermal: dict[str, str] | None = None  # nuclide -> S(alpha, beta) law
    macroscopic: str | None = None  # name in the model's multigroup library

    def thermal_law(self, nuclide: str) -> str:
        return (self.thermal or {}).get(nuclide, "free-gas")


@dataclass(frozen=True)
class Pin:
    """Concentric rings in a square cell; the last material fills the rest."""

    name: str
    radii: tuple[float, ...]  # cm, increasing
    materials: tuple[str, ...]  # len(radii) + 1

    def areas(self, pitch: float) -> list[tuple[str, float]]:
        """(material, area in cm^2) of every ring, the last one up to the cell edge."""
        out, inner = [], 0.0
        for r, m in zip(self.radii, self.materials, strict=False):
            out.append((m, math.pi * (r * r - inner * inner)))
            inner = r
        out.append((self.materials[-1], pitch * pitch - math.pi * inner * inner))
        return out


@dataclass(frozen=True)
class Lattice:
    """One xy slice: a square map of pins, or a single homogeneous material."""

    name: str
    geometry: Literal["square"]
    pitch: float | None = None  # pin pitch, cm
    pins: tuple[tuple[str, ...], ...] | None = None  # rows of pin names
    material: str | None = None  # homogeneous slice (reflectors)
    gap: float = (
        0.0  # water around the pins on every side (half the inter-assembly gap)
    )

    @property
    def width(self) -> float | None:
        """Side length (cm) of a pin lattice; None for a homogeneous slice."""
        if self.pins is None or self.pitch is None:
            return None
        return len(self.pins) * self.pitch + 2 * self.gap


@dataclass(frozen=True)
class AxialSegment:
    lattice: str
    height: float  # cm


@dataclass(frozen=True)
class Assembly:
    """A stack of lattices, listed from bottom to top."""

    name: str
    segments: tuple[AxialSegment, ...]

    @property
    def height(self) -> float:
        return sum(s.height for s in self.segments)


@dataclass(frozen=True)
class Core:
    """Assemblies on a square map; None marks a position outside the core.

    ``boundary`` gives the condition on each side of the mapped region
    ("reflective" or "vacuum"; default vacuum).

    With ``symmetry = "quarter"`` the map is the south-east quadrant: its
    first row and column lie on the symmetry lines (an odd-sized core, with
    the centre assembly in the top-left corner).
    """

    geometry: Literal["square"]
    symmetry: Literal["full", "quarter"]
    assembly_pitch: float  # cm
    map: tuple[tuple[str | None, ...], ...]
    boundary: dict[str, str] = field(
        default_factory=lambda: dict.fromkeys(SIDES, "vacuum")
    )

    def full_map(self) -> tuple[tuple[str | None, ...], ...]:
        if self.symmetry == "full":
            return self.map
        rows = [tuple(reversed(r[1:])) + r for r in self.map]
        return tuple(reversed(rows[1:])) + tuple(rows)


@dataclass(frozen=True)
class Reference:
    """A reference result for a state (a benchmark value, not an input)."""

    keff: float
    sigma: float
    source: str


@dataclass(frozen=True)
class State:
    name: str
    fuel_temperature: float  # K
    moderator_temperature: float  # K (clad and gap follow it)
    moderator_density: float  # g/cm3
    boron_ppm: float  # parts per million by weight
    reference: Reference | None = None

    def temperature(self, which: Literal["fuel", "moderator"]) -> float:
        return self.fuel_temperature if which == "fuel" else self.moderator_temperature


# The state variables a branch can change, in the order of the table axes.
BRANCH_AXES = (
    "fuel_temperature",
    "moderator_temperature",
    "moderator_density",
    "boron_ppm",
)


@dataclass(frozen=True)
class Branches:
    """A state matrix: the lattice calculation is repeated at every point of
    the tensor product of the axes, so that a nodal code can interpolate its
    constants in the local conditions (M2, docs chapter 20).

    axes holds a tuple for every name of BRANCH_AXES; an axis the file does
    not list has the base state's value alone."""

    name: str
    base: str  # the state the matrix is built around
    lattice: str | None  # None: the model's only lattice
    axes: dict[str, tuple[float, ...]]

    @property
    def shape(self) -> tuple[int, ...]:
        return tuple(len(self.axes[a]) for a in BRANCH_AXES)

    def states(self) -> list[State]:
        """Every point as a State, in C order of the axes (boron fastest)."""
        out = []
        for index in itertools.product(*(range(n) for n in self.shape)):
            fuel, mod_t, mod_rho, ppm = (
                self.axes[a][i] for a, i in zip(BRANCH_AXES, index, strict=True)
            )
            label = "-".join(str(i) for i in index)
            out.append(State(f"{self.name}-{label}", fuel, mod_t, mod_rho, ppm))
        return out


@dataclass(frozen=True)
class Model:
    name: str
    description: str
    materials: dict[str, Material]
    pins: dict[str, Pin]
    lattices: dict[str, Lattice]
    assemblies: dict[str, Assembly]
    core: Core | None
    states: dict[str, State]
    multigroup: dict[str, MacroXS] = field(default_factory=dict)  # macroscopic data
    branches: dict[str, Branches] = field(default_factory=dict)  # state matrices
