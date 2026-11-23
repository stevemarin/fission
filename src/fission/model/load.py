"""Read a model from a TOML file and check it.

Every error names where in the file it is, e.g.

    vera-p1.toml: lattices.cell.map row 2: 3 entries, row 1 has 2

Checks: unknown keys (typos), wrong types, missing or unknown references
(a pin that names an undefined material, ...), ring radii that do not
increase or do not fit in the cell, ragged maps, isotopic percentages that
do not add up to 100, isotopes that do not belong to their element, and
options that are reserved but not implemented yet (hexagonal geometry,
eighth-core symmetry).
"""

from __future__ import annotations

import itertools
import tomllib
from pathlib import Path
from typing import Any, NoReturn, cast

from .materials import element_of
from .multigroup import LibraryError, MacroXS, load_multigroup
from .spec import (
    BRANCH_AXES,
    SIDES,
    Assembly,
    AxialSegment,
    Banks,
    Branches,
    ControlRod,
    Core,
    CoreThermal,
    Isotopics,
    Lattice,
    Material,
    Model,
    Pin,
    Reference,
    State,
)

FORMAT = 1
MODELS_DIR = Path(__file__).resolve().parents[3] / "models"


class ModelError(ValueError):
    pass


_REQUIRED: Any = object()


class _Table:
    """A TOML table being read. Keys that are never asked for are errors."""

    def __init__(self, data: object, where: str):
        if not isinstance(data, dict):
            _fail(where, f"expected a table, got {_kind(data)}")
        self.data: dict[str, Any] = data
        self.where = where
        self.known: set[str] = set()

    def at(self, key: str) -> str:
        return f"{self.where}.{key}" if self.where else key

    def has(self, key: str) -> bool:
        self.known.add(key)
        return key in self.data

    def raw(self, key: str, default: Any = _REQUIRED) -> Any:
        self.known.add(key)
        if key not in self.data:
            if default is _REQUIRED:
                _fail(self.where or "(top level)", f"missing required key '{key}'")
            return default
        return self.data[key]

    def number(
        self, key: str, default: Any = _REQUIRED, minimum: float | None = 0.0
    ) -> float:
        v = self.raw(key, default)
        if v is default and default is not _REQUIRED:
            return cast(float, v)
        if isinstance(v, bool) or not isinstance(v, (int, float)):
            _fail(self.at(key), f"expected a number, got {_kind(v)}")
        if minimum is not None and not v > minimum:
            _fail(self.at(key), f"must be > {minimum:g}, got {v:g}")
        return float(v)

    def string(
        self, key: str, default: Any = _REQUIRED, choices: tuple[str, ...] = ()
    ) -> str:
        v = self.raw(key, default)
        if v is default and default is not _REQUIRED:
            return cast(str, v)
        if not isinstance(v, str):
            _fail(self.at(key), f"expected a string, got {_kind(v)}")
        if choices and v not in choices:
            _fail(self.at(key), f"{v!r} is not one of {', '.join(map(repr, choices))}")
        return v

    def boolean(self, key: str, default: Any = _REQUIRED) -> bool:
        v = self.raw(key, default)
        if not isinstance(v, bool):
            _fail(self.at(key), f"expected true or false, got {_kind(v)}")
        return v

    def table(self, key: str, default: Any = _REQUIRED) -> _Table:
        v = self.raw(key, default)
        return _Table(v, self.at(key))

    def numbers(self, key: str) -> dict[str, float]:
        """A table of non-negative numbers (nuclide -> value, element -> count)."""
        t = self.table(key)
        out = {k: t.number(k, minimum=None) for k in t.data}
        for k, v in out.items():
            if v < 0:
                _fail(t.at(k), f"must be >= 0, got {v:g}")
        if not out:
            _fail(t.where, "is empty")
        return out

    def finish(self) -> None:
        extra = sorted(set(self.data) - self.known)
        if extra:
            expected = ", ".join(sorted(self.known)) or "nothing"
            _fail(
                self.where or "(top level)",
                f"unknown key(s) {extra}; expected {expected}",
            )


def _kind(v: object) -> str:
    return {dict: "a table", list: "an array", str: "a string", bool: "a boolean"}.get(
        type(v), type(v).__name__
    )


def _fail(where: str, message: str) -> NoReturn:
    raise ModelError(f"{where}: {message}")


def _reserved(where: str, value: str, implemented: tuple[str, ...]) -> None:
    if value not in implemented:
        _fail(
            where,
            f"{value!r} is reserved but not implemented yet (have: {', '.join(implemented)})",
        )


# ---------------------------------------------------------------- materials
def _isotopics(t: _Table, element: str, basis: str) -> Isotopics:
    percent = t.numbers(element)
    for nuclide in percent:
        try:
            owner = element_of(nuclide)
        except ValueError:
            _fail(f"{t.at(element)}.{nuclide}", "not a nuclide name (like 'U235')")
        if owner != element:
            _fail(f"{t.at(element)}.{nuclide}", f"is not an isotope of {element}")
    total = sum(percent.values())
    if abs(total - 100.0) > 1e-3:
        _fail(t.at(element), f"percentages add up to {total:g}, not 100")
    return Isotopics("atom" if basis == "atom_percent" else "weight", percent)


def _element_isotopics(m: _Table) -> dict[str, Isotopics]:
    out: dict[str, Isotopics] = {}
    for basis in ("atom_percent", "weight_percent"):
        if m.has(basis):
            t = m.table(basis)
            for element in t.data:
                if element in out:
                    _fail(t.at(element), "isotopics given both by atom and by weight")
                out[element] = _isotopics(t, element, basis)
            t.finish()
    return out


def _material(name: str, m: _Table) -> Material:
    if m.has("table"):
        temperature_t: Any = m.string("temperature", "moderator", ("fuel", "moderator"))
        table = m.string("table")
        m.finish()
        return Material(name, temperature_t, table=table)
    if m.has("macroscopic"):
        temperature: Any = m.string("temperature", "moderator", ("fuel", "moderator"))
        macroscopic = m.string("macroscopic")
        m.finish()  # anything else (densities, formula, ...) is an unknown key here
        return Material(name, temperature, macroscopic=macroscopic)
    temperature = m.string("temperature", choices=("fuel", "moderator"))
    shielded = m.boolean("shielded", False)
    thermal = None
    if m.has("thermal"):
        t = m.table("thermal")
        thermal = {k: t.string(k) for k in t.data}
    if m.has("atom_densities"):
        for key in (
            "density",
            "formula",
            "atom_percent",
            "weight_percent",
            "soluble_boron",
        ):
            if m.has(key):
                _fail(m.at(key), "not allowed together with atom_densities")
        n = m.numbers("atom_densities")
        for nuclide in n:
            try:
                element_of(nuclide)
            except ValueError:
                _fail(
                    m.at(f"atom_densities.{nuclide}"),
                    "not a nuclide name (like 'U235')",
                )
        material = Material(name, temperature, atom_densities=n, thermal=thermal,
                            shielded=shielded)  # fmt: skip
    else:
        rho: float | str
        if isinstance(m.raw("density"), str):
            rho = m.string("density", choices=("moderator",))
        else:
            rho = m.number("density")
        formula = m.numbers("formula")
        isotopics = _element_isotopics(m)
        for element in formula:
            if element not in isotopics:
                _fail(m.at("formula"), f"no isotopics for element {element!r} "
                      f"(add atom_percent.{element} or weight_percent.{element})")  # fmt: skip
        for element in isotopics:
            if element not in formula:
                _fail(
                    m.where, f"isotopics for {element!r}, which is not in the formula"
                )
        boron = None
        if m.has("soluble_boron"):
            b = m.table("soluble_boron")
            found = _element_isotopics(b)
            boron = found.get("B")
            if boron is None or len(found) != 1:
                _fail(b.where, "expected atom_percent.B or weight_percent.B only")
            b.finish()
            if rho != "moderator":
                _fail(
                    m.at("soluble_boron"),
                    'needs density = "moderator" (boron is a state variable)',
                )
        material = Material(
            name, temperature, density=cast(Any, rho), formula=formula,
            isotopics=isotopics, soluble_boron=boron, thermal=thermal, shielded=shielded,
        )  # fmt: skip
    m.finish()
    return material


def _nuclides(material: Material) -> set[str]:
    if material.atom_densities is not None:
        return set(material.atom_densities)
    out = {n for iso in (material.isotopics or {}).values() for n in iso.percent}
    if material.soluble_boron is not None:
        out |= set(material.soluble_boron.percent)
    return out


# ---------------------------------------------------------------- geometry
def _grid(
    t: _Table, legend: dict[str, str | None]
) -> tuple[tuple[str | None, ...], ...]:
    text = t.string("map")
    rows = [line.split() for line in text.splitlines() if line.strip()]
    if not rows:
        _fail(t.at("map"), "is empty")
    for i, row in enumerate(rows, 1):
        if len(row) != len(rows[0]):
            _fail(
                t.at("map"), f"row {i} has {len(row)} entries, row 1 has {len(rows[0])}"
            )
        for token in row:
            if token not in legend:
                _fail(t.at("map"), f"row {i}: {token!r} is not in the legend")
    return tuple(tuple(legend[token] for token in row) for row in rows)


def _legend(t: _Table) -> dict[str, str | None]:
    lt = t.table("legend")
    legend: dict[str, str | None] = {k: lt.string(k) for k in lt.data}
    lt.finish()
    if "." in legend:
        _fail(lt.at("."), "'.' is reserved for an empty position")
    return legend


def _pin(name: str, t: _Table, materials: dict[str, Material]) -> Pin:
    radii_raw = t.raw("radii")
    names = t.raw("materials")
    if not isinstance(radii_raw, list) or not all(
        isinstance(r, (int, float)) and not isinstance(r, bool) for r in radii_raw
    ):
        _fail(t.at("radii"), "expected an array of numbers")
    if not isinstance(names, list) or not all(isinstance(n, str) for n in names):
        _fail(t.at("materials"), "expected an array of material names")
    radii = tuple(float(r) for r in radii_raw)
    if any(r <= 0 for r in radii) or any(b <= a for a, b in itertools.pairwise(radii)):
        _fail(t.at("radii"), f"must be positive and increasing, got {list(radii)}")
    if len(names) != len(radii) + 1:
        _fail(t.at("materials"), f"{len(radii)} radii need {len(radii) + 1} materials "
              f"(the last fills the cell), got {len(names)}")  # fmt: skip
    for n in names:
        if n not in materials:
            _fail(t.at("materials"), f"unknown material {n!r}")
    t.finish()
    return Pin(name, radii, tuple(names))


def _lattice(
    name: str, t: _Table, materials: dict[str, Material], pins: dict[str, Pin]
) -> Lattice:
    geometry = t.string("geometry", "square")
    _reserved(t.at("geometry"), geometry, ("square",))
    if t.has("material"):
        for key in ("pitch", "legend", "map", "gap"):
            if t.has(key):
                _fail(
                    t.at(key),
                    "not allowed in a homogeneous lattice (one with 'material')",
                )
        material = t.string("material")
        if material not in materials:
            _fail(t.at("material"), f"unknown material {material!r}")
        rodded = {}
        if t.has("rodded"):
            r = t.table("rodded")
            rodded = {k: r.string(k) for k in r.data}
            r.finish()
        t.finish()
        return Lattice(name, "square", material=material, rodded=rodded)
    pitch = t.number("pitch")
    gap = t.number("gap", 0.0, minimum=None)
    if gap < 0:
        _fail(t.at("gap"), f"must be >= 0, got {gap:g}")
    if t.has("fill"):
        for key in ("legend", "map"):
            if t.has(key):
                _fail(t.at(key), "not allowed together with 'fill'")
        fill = t.string("fill")
        size = t.raw("size")
        if isinstance(size, bool) or not isinstance(size, int) or size < 1:
            _fail(t.at("size"), "expected a positive integer (pins per side)")
        legend: dict[str, str | None] = {"fill": fill}
    else:
        legend = _legend(t)
    for key, pin in legend.items():
        where = t.at("fill" if key == "fill" else f"legend.{key}")
        if pin is None or pin not in pins:
            _fail(where, f"unknown pin {pin!r}")
        if pins[pin].radii and pins[pin].radii[-1] > pitch / 2:
            _fail(where, f"pin {pin!r} (radius {pins[pin].radii[-1]:g} cm) "
                  f"does not fit in the pitch {pitch:g} cm")  # fmt: skip
    if "fill" in legend:
        grid: tuple[tuple[str | None, ...], ...] = ((legend["fill"],) * size,) * size
    else:
        grid = _grid(t, legend)
    if gap and len(grid) != len(grid[0]):
        _fail(t.at("map"), f"a lattice with a gap must be square, got {len(grid)} rows "
              f"of {len(grid[0])}")  # fmt: skip
    t.finish()
    return Lattice(
        name, "square", pitch, cast(tuple[tuple[str, ...], ...], grid), gap=gap
    )


def _assembly(name: str, t: _Table, lattices: dict[str, Lattice]) -> Assembly:
    axial = t.raw("axial")
    if not isinstance(axial, list) or not axial:
        _fail(t.at("axial"), "expected a non-empty array of {lattice, height} tables")
    segments = []
    for i, entry in enumerate(axial, 1):
        s = _Table(entry, f"{t.at('axial')}[{i}]")
        lattice = s.string("lattice")
        if lattice not in lattices:
            _fail(s.at("lattice"), f"unknown lattice {lattice!r}")
        segments.append(AxialSegment(lattice, s.number("height")))
        s.finish()
    t.finish()
    return Assembly(name, tuple(segments))


def _core_thermal(t: _Table) -> CoreThermal:
    th = CoreThermal(
        t.number("assembly_power"),
        t.number("assembly_flow"),
        t.number("inlet_temperature"),
        t.number("pressure"),
        t.number("fuel_heat_fraction"),
        int(t.number("rods")),
        t.number("pellet_radius"),
        t.number("clad_inner"),
        t.number("clad_outer"),
        t.number("gap_conductance"),
        t.number("rod_pitch", None),
        t.number("flow_area", None),
        t.number("hydraulic_diameter", None),
    )
    film = (th.rod_pitch, th.flow_area, th.hydraulic_diameter)
    if any(x is None for x in film) and any(x is not None for x in film):
        _fail(t.where, "rod_pitch, flow_area and hydraulic_diameter go together")
    if th.fuel_heat_fraction > 1:
        _fail(t.at("fuel_heat_fraction"), "must be at most 1")
    if not th.pellet_radius < th.clad_inner < th.clad_outer:
        _fail(t.where, "need pellet_radius < clad_inner < clad_outer")
    t.finish()
    return th


def _control_rod(name: str, t: _Table) -> ControlRod:
    zones_raw = t.raw("zones")
    if not isinstance(zones_raw, list) or not zones_raw:
        _fail(t.at("zones"), "expected a non-empty array of {absorber, height} tables")
    zones = []
    for i, entry in enumerate(zones_raw, 1):
        z = _Table(entry, f"{t.at('zones')}[{i}]")
        zones.append((z.string("absorber"), z.number("height")))
        z.finish()
    steps = t.raw("steps")
    if isinstance(steps, bool) or not isinstance(steps, int) or steps < 1:
        _fail(t.at("steps"), "expected a positive integer")
    rod = ControlRod(name, t.number("bottom", minimum=None), t.number("step"), steps,
                     tuple(zones))  # fmt: skip
    t.finish()
    return rod


def _banks(
    t: _Table, grid: tuple[tuple[str | None, ...], ...], rods: dict[str, ControlRod]
) -> Banks:
    rod = t.string("rod")
    if rod not in rods:
        _fail(t.at("rod"), f"unknown control rod {rod!r} (see [control_rods])")
    text = t.string("map")
    rows = [line.split() for line in text.splitlines() if line.strip()]
    if len(rows) != len(grid) or any(len(r) != len(grid[0]) for r in rows):
        _fail(t.at("map"), f"must have the core map's shape, {len(grid)} rows of "
              f"{len(grid[0])}")  # fmt: skip
    positions: dict[str, list[tuple[int, int]]] = {}
    for i, row in enumerate(rows):
        for j, token in enumerate(row):
            if token == ".":
                continue
            if grid[i][j] is None:
                _fail(t.at("map"), f"row {i + 1}: bank {token!r} outside the core")
            positions.setdefault(token, []).append((i, j))

    def worths(key: str) -> dict[str, float]:
        if not t.has(key):
            return {}
        w = t.table(key)
        out = {k: w.number(k) for k in w.data}
        for k in out:
            if k not in positions:
                _fail(w.at(k), f"no bank {k!r} in the map")
        return out

    banks = Banks(rod, {k: tuple(v) for k, v in positions.items()}, worths("measured"),
                  worths("reference"), t.string("source", ""))  # fmt: skip
    t.finish()
    return banks


def _core(
    t: _Table, assemblies: dict[str, Assembly], lattices: dict[str, Lattice],
    rods: dict[str, ControlRod],
) -> Core:  # fmt: skip
    geometry = t.string("geometry", "square")
    _reserved(t.at("geometry"), geometry, ("square",))
    symmetry = t.string("symmetry", "full")
    _reserved(t.at("symmetry"), symmetry, ("full", "quarter"))
    pitch = t.number("assembly_pitch")
    boundary = dict.fromkeys(SIDES, "vacuum")
    if t.has("boundary"):
        b = t.table("boundary")
        for side in SIDES:
            boundary[side] = b.string(side, "vacuum", ("vacuum", "reflective"))
        b.finish()
    vacuum = t.number("vacuum_coefficient", 0.5)
    state = t.string("state", None)
    thermal = _core_thermal(t.table("thermal")) if t.has("thermal") else None
    legend = _legend(t)
    for key, a in legend.items():
        if a not in assemblies:
            _fail(t.at(f"legend.{key}"), f"unknown assembly {a!r}")
    legend["."] = None
    grid = _grid(t, legend)
    used = {a for row in grid for a in row if a is not None}
    heights = {a: assemblies[a].height for a in used}
    if len(set(heights.values())) > 1:
        _fail(t.at("map"), "assemblies in the core must have the same total height, got "
              + ", ".join(f"{a} {h:g} cm" for a, h in sorted(heights.items())))  # fmt: skip
    for a in sorted(used):
        for s in assemblies[a].segments:
            width = lattices[s.lattice].width
            if width is not None and width > pitch:
                _fail(t.at("assembly_pitch"), f"lattice {s.lattice!r} of assembly {a!r} "
                      f"is {width:g} cm wide, more than the assembly pitch {pitch:g} cm")  # fmt: skip
    banks = _banks(t.table("banks"), grid, rods) if t.has("banks") else None
    t.finish()
    return Core(cast(Any, geometry), cast(Any, symmetry), pitch, grid, boundary, vacuum,
                state, thermal, banks)  # fmt: skip


def _state(name: str, t: _Table) -> State:
    reference = None
    if t.has("reference"):
        r = t.table("reference")
        reference = Reference(
            r.number("keff"), r.number("sigma", minimum=None), r.string("source")
        )
        r.finish()
    state = State(
        name,
        t.number("fuel_temperature"),
        t.number("moderator_temperature"),
        t.number("moderator_density"),
        t.number("boron_ppm", minimum=None),
        reference,
    )
    if state.boron_ppm < 0:
        _fail(t.at("boron_ppm"), "must be >= 0")
    t.finish()
    return state


def _branches(
    name: str, t: _Table, states: dict[str, State], lattices: dict[str, Lattice]
) -> Branches:
    base_name = t.string("base")
    if base_name not in states:
        _fail(t.at("base"), f"unknown state {base_name!r} (have {', '.join(states)})")
    base = states[base_name]
    lattice = t.string("lattice", None)
    if lattice is not None and lattice not in lattices:
        _fail(t.at("lattice"), f"unknown lattice {lattice!r}")
    axes: dict[str, tuple[float, ...]] = {}
    for axis in BRANCH_AXES:
        if not t.has(axis):
            axes[axis] = (getattr(base, axis),)
            continue
        values = t.raw(axis)
        if not isinstance(values, list) or not values:
            _fail(
                t.at(axis), f"expected a non-empty list of numbers, got {_kind(values)}"
            )
        for v in values:
            if isinstance(v, bool) or not isinstance(v, (int, float)):
                _fail(t.at(axis), f"expected numbers, got {_kind(v)}")
        lowest = 0.0 if axis == "boron_ppm" else 1e-300
        if min(values) < lowest:
            _fail(
                t.at(axis), "values must be positive" + (" or 0" if lowest == 0 else "")
            )
        if any(b <= a for a, b in itertools.pairwise(values)):
            _fail(t.at(axis), "values must increase (they are interpolation nodes)")
        axes[axis] = tuple(float(v) for v in values)
        if not any(
            abs(v - getattr(base, axis)) <= 1e-9 * max(1.0, abs(v)) for v in values
        ):
            _fail(t.at(axis), f"must contain the base state's value {getattr(base, axis):g}"
                  " (the coefficients are taken around it)")  # fmt: skip
    t.finish()
    return Branches(name, base_name, lattice, axes)


# ---------------------------------------------------------------- model
def _section(top: _Table, key: str) -> _Table:
    return top.table(key, {}) if top.has(key) else _Table({}, key)


def parse(data: dict[str, Any], base: Path | None = None) -> Model:
    """Build a Model from parsed TOML (a dict), checking everything.

    base: directory that relative paths in the file (``multigroup``) refer to.
    """
    top = _Table(data, "")
    head = top.table("model")
    fmt = head.raw("format")
    if fmt != FORMAT:
        _fail(head.at("format"), f"this code reads format {FORMAT}, got {fmt!r}")
    name = head.string("name")
    description = head.string("description", "")
    library: dict[str, MacroXS] = {}
    if head.has("multigroup"):
        path = Path(head.string("multigroup"))
        path = (base or MODELS_DIR) / path if not path.is_absolute() else path
        try:
            library = load_multigroup(path)
        except (OSError, LibraryError) as e:
            _fail(head.at("multigroup"), str(e))
    head.finish()

    sec = _section(top, "materials")
    materials = {k: _material(k, sec.table(k)) for k in sec.data}
    for mat in materials.values():
        for nuclide in mat.thermal or {}:
            if nuclide not in _nuclides(mat):
                _fail(
                    f"materials.{mat.name}.thermal.{nuclide}",
                    "nuclide not in this material",
                )
    for mat in materials.values():
        if mat.macroscopic is not None and mat.macroscopic not in library:
            have = ", ".join(library) or "none: set [model] multigroup"
            _fail(f"materials.{mat.name}.macroscopic",
                  f"{mat.macroscopic!r} is not in the multigroup library (have {have})")  # fmt: skip
    sec = _section(top, "pins")
    pins = {k: _pin(k, sec.table(k), materials) for k in sec.data}
    sec = _section(top, "lattices")
    lattices = {k: _lattice(k, sec.table(k), materials, pins) for k in sec.data}
    for lat in lattices.values():
        for absorber, other in lat.rodded.items():
            if other not in lattices:
                _fail(
                    f"lattices.{lat.name}.rodded.{absorber}",
                    f"unknown lattice {other!r}",
                )
    sec = _section(top, "assemblies")
    assemblies = {k: _assembly(k, sec.table(k), lattices) for k in sec.data}
    sec = _section(top, "control_rods")
    rods = {k: _control_rod(k, sec.table(k)) for k in sec.data}
    core = (
        _core(top.table("core"), assemblies, lattices, rods)
        if top.has("core")
        else None
    )
    sec = _section(top, "states")
    states = {k: _state(k, sec.table(k)) for k in sec.data}
    if core is not None and core.state is not None and core.state not in states:
        _fail("core.state", f"unknown state {core.state!r}")
    for mat in materials.values():
        if mat.density == "moderator" and not states:
            _fail(
                f"materials.{mat.name}.density", '"moderator" needs at least one state'
            )
    sec = _section(top, "branches")
    branches = {k: _branches(k, sec.table(k), states, lattices) for k in sec.data}
    top.finish()
    return Model(
        name, description, materials, pins, lattices, assemblies, core, states, library,
        branches, rods,
    )  # fmt: skip


def _resolve(path: str | Path, near: Path | None = None) -> Path:
    """A model file: as given, next to `near`, or a bare name in models/."""
    p = Path(path)
    for candidate in (p, (near / p) if near else None, MODELS_DIR / p,
                      MODELS_DIR / f"{p}.toml", (near / f"{p}.toml") if near else None):  # fmt: skip
        if candidate is not None and candidate.exists():
            return candidate
    return p


def merge(base: dict, override: dict) -> dict:
    """Tables merged key by key (recursively); any other value replaces."""
    out = dict(base)
    for key, value in override.items():
        if isinstance(value, dict) and isinstance(out.get(key), dict):
            out[key] = merge(out[key], value)
        else:
            out[key] = value
    return out


def _read(p: Path, seen: tuple[Path, ...] = ()) -> dict:
    """The TOML of a model file, with `[model] extends = "<base model>"`
    resolved: the base file's tables, with this file's merged over them."""
    if p.resolve() in seen:
        raise ModelError(f"{p.name}: model.extends: circular ({p.name} extends itself)")
    with open(p, "rb") as f:
        data = tomllib.load(f)
    head = data.get("model", {})
    if not isinstance(head, dict) or "extends" not in head:
        return data
    base_path = _resolve(str(head["extends"]), p.resolve().parent)
    if not base_path.exists():
        raise ModelError(f"model.extends: no model {head['extends']!r}")
    base = _read(base_path, (*seen, p.resolve()))
    data = {**data, "model": {k: v for k, v in head.items() if k != "extends"}}
    return merge(base, data)


def load_model(path: str | Path) -> Model:
    """Read and check a model file. A bare name is looked up in models/.
    A file with ``[model] extends = "<other model>"`` adds to (and overrides,
    table by table) the other model."""
    p = _resolve(path)
    try:
        data = _read(p)
        return parse(data, p.resolve().parent)
    except tomllib.TOMLDecodeError as e:
        raise ModelError(f"{p.name}: {e}") from None
    except ModelError as e:
        raise ModelError(f"{p.name}: {e}") from None
