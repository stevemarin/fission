"""The state matrix (M2): the lattice calculation at every branch point.

A nodal code solves the whole core with a few homogenized nodes per
assembly. Each node sits at its own conditions: fuel temperature (power),
coolant temperature and density (how far up the channel it is), soluble
boron (the same everywhere, but changed through the cycle). Its few-group
constants must follow them, and these four variables are where the
reactivity feedbacks of a PWR come from:

* fuel temperature: Doppler broadening of the U-238 resonances raises
  resonance capture (chapter 8), the prompt negative feedback;
* moderator density: less water moderates less and absorbs less; in an
  undermoderated lattice the net effect lowers k;
* moderator temperature: at fixed density, a hotter Maxwellian shifts the
  thermal spectrum up, where U-235 fission and the 1/v absorbers drop;
* boron: a thermal absorber in the water.

The lattice code computes constants at a set of points, and the nodal code
interpolates between them. Here the points are the tensor product of the
values of the model file's ``[branches.<name>]`` (``Branches``), so that
M3 can interpolate multilinearly. Every point is an ordinary lattice
calculation of chapter 18: solve, homogenize, B1 critical spectrum,
condense to few groups.

Two lattice solvers: ``infinite`` (the homogenized cell with narrow-
resonance shielding of chapter 12, seconds per point, the first one to use)
and ``moc`` (the heterogeneous pin of chapters 16-19, minutes per point).

    uv run python -m fission.branches vera
    uv run python -m fission.branches vera --solver moc --rings --sph --clad
"""

from __future__ import annotations

import json
import time
from collections.abc import Callable, Iterator
from dataclasses import dataclass, replace
from pathlib import Path

import numpy as np

from fission.data import EndfLibrary
from fission.fewgroup.b1 import critical, mode
from fission.fewgroup.condense import THERMAL_CUT, FewGroup, condense
from fission.fewgroup.homogenize import FineConstants, from_mixture, homogenize
from fission.model import BRANCH_AXES, Branches, Model, State

FORMAT = 1
SOLVERS = ("infinite", "moc")
OUT = Path("endf-data/branches")


def lattice_constants(
    model: Model, lattice: str, state: State, solver: str, options: dict
) -> tuple[float, FineConstants]:
    """Solve the lattice at one state: (k-infinity, fine-group constants)."""
    m = replace(model, states={**model.states, state.name: state})
    if solver == "infinite":
        from fission.pincell.run import run

        sol = run(state.name, options.get("shielding", "nr"), m, lattice)
        return sol.k, from_mixture(sol, state.name)
    if solver == "moc":
        from fission.pincell.heterogeneous import run as moc

        r = moc(state.name, model=m, lattice=lattice, **options)
        return r.solution.k, homogenize(r, state.name)
    raise ValueError(f"unknown solver {solver!r} (have {', '.join(SOLVERS)})")


def few_group(fine: FineConstants, spectrum: str, cuts: tuple[float, ...]) -> FewGroup:
    weight = critical(fine) if spectrum == "b1" else mode(fine, 0.0)
    return condense(fine, weight, cuts)


def _point_record(index: tuple[int, ...], state: State, k: float, fg: FewGroup,
                  seconds: float) -> dict:  # fmt: skip
    return {
        "index": list(index),
        "state": {a: getattr(state, a) for a in BRANCH_AXES},
        "k_infinity": k,
        "few_group_k_infinity": fg.k(0.0),
        "buckling_cm-2": fg.buckling,
        "flux": fg.flux.tolist(),
        "D": fg.D.tolist(),
        "total": fg.total.tolist(),
        "sigma_a": fg.sigma_a.tolist(),
        "sigma_s": fg.sigma_s.tolist(),
        "nu_sigma_f": fg.nu_sigma_f.tolist(),
        "sigma_f": fg.sigma_f.tolist(),
        "kappa_sigma_f": fg.kappa_sigma_f.tolist(),
        "chi": fg.chi.tolist(),
        "adf": fg.adf.tolist(),
        "seconds": seconds,
    }


@dataclass
class Table:
    """Few-group constants on the grid of a state matrix."""

    header: dict  # model, branches, solver, options, axes, ...
    points: dict[tuple[int, ...], dict]  # index -> record

    @property
    def axes(self) -> dict[str, list[float]]:
        return self.header["axes"]

    @property
    def shape(self) -> tuple[int, ...]:
        return tuple(len(self.axes[a]) for a in BRANCH_AXES)

    @property
    def complete(self) -> bool:
        return len(self.points) == int(np.prod(self.shape))

    def array(self, key: str) -> np.ndarray:
        """A quantity on the grid: shape = axes shape + the quantity's own
        shape (e.g. groups, or groups x groups for sigma_s)."""
        if not self.complete:
            raise ValueError("the table is not complete")
        first = np.asarray(self.points[next(iter(self.points))][key], dtype=float)
        out = np.zeros(self.shape + first.shape)
        for index, record in self.points.items():
            out[index] = record[key]
        return out

    def write(self, path: Path) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        data = {**self.header, "points": [self.points[i] for i in sorted(self.points)]}
        tmp = path.with_suffix(".tmp")
        tmp.write_text(json.dumps(data, indent=1))
        tmp.replace(path)  # never leave a half-written table

    @classmethod
    def read(cls, path: Path) -> Table:
        data = json.loads(path.read_text())
        points = {tuple(p["index"]): p for p in data.pop("points")}
        return cls(data, points)


def header(model: Model, branches: Branches, lattice: str, solver: str,
           options: dict, spectrum: str, cuts: tuple[float, ...]) -> dict:  # fmt: skip
    base = model.states[branches.base]
    return {
        "format": FORMAT,
        "model": model.name,
        "branches": branches.name,
        "lattice": lattice,
        "solver": solver,
        "options": options,
        "spectrum": spectrum,
        "boundaries_eV": None,  # filled from the first point
        "cuts_eV": list(cuts),
        "axes": {a: list(branches.axes[a]) for a in BRANCH_AXES},
        "base": {"state": base.name, **{a: getattr(base, a) for a in BRANCH_AXES}},
        "units": "1/cm; D in cm; kappa_sigma_f in eV/cm; temperatures K; "
        "moderator density g/cm3; boron ppm",
    }


def _same_run(a: dict, b: dict) -> bool:
    keys = ("format", "model", "branches", "lattice", "solver", "options",
            "spectrum", "cuts_eV", "axes")  # fmt: skip
    return all(a.get(k) == b.get(k) for k in keys)


def run_matrix(
    model: Model,
    name: str,
    solver: str = "infinite",
    options: dict | None = None,
    spectrum: str = "b1",
    cuts: tuple[float, ...] = (THERMAL_CUT,),
    path: Path | None = None,
    fresh: bool = False,
    keep_libraries: bool = False,
    progress: Callable[[str], None] | None = None,
) -> Table:
    """Run every point of ``model.branches[name]`` that is not in `path` yet.

    The table is rewritten after every point, so an interrupted run resumes
    where it stopped (unless the settings changed, or `fresh`). The shaped
    libraries of a point are deleted after it, unless `keep_libraries`:
    every point has its own self-shielding, so they would only fill the
    cache (about 2 MB per nuclide and point)."""
    if solver not in SOLVERS:
        raise ValueError(f"unknown solver {solver!r} (have {', '.join(SOLVERS)})")
    branches = model.branches[name]
    lattice = branches.lattice or next(iter(model.lattices))
    options = dict(options or {})
    head = header(model, branches, lattice, solver, options, spectrum, cuts)
    table = Table(head, {})
    if path is not None and path.exists() and not fresh:
        old = Table.read(path)
        if _same_run(old.header, head):
            table = old
    for index, state in _indexed(branches):
        if index in table.points:
            continue
        t0 = time.perf_counter()
        k, fine = lattice_constants(model, lattice, state, solver, options)
        fg = few_group(fine, spectrum, cuts)
        table.header["boundaries_eV"] = list(fg.boundaries)
        table.points[index] = _point_record(index, state, k, fg,
                                            time.perf_counter() - t0)  # fmt: skip
        if not keep_libraries:
            _drop_shaped_libraries(state.name)
        if path is not None:
            table.write(path)
        if progress is not None:
            progress(
                f"{state.name}: fuel {state.fuel_temperature:g} K, moderator "
                f"{state.moderator_temperature:g} K {state.moderator_density:g} g/cm3, "
                f"boron {state.boron_ppm:g} ppm: k-inf {k:.5f} "
                f"({time.perf_counter() - t0:.0f} s, {len(table.points)}/"
                f"{int(np.prod(branches.shape))})"
            )
    return table


# ---------------------------------------------------------------- check points
def sample_states(branches: Branches, count: int, seed: int = 1) -> list[State]:
    """Random states inside the box of the matrix (uniform in each variable,
    off the grid), to test interpolation. Moderator temperatures are rounded
    to 5 K, the spacing of the H in H2O scattering law, fuel temperatures to
    1 K."""
    rng = np.random.default_rng(seed)
    out = []
    for n in range(count):
        values = []
        for axis in BRANCH_AXES:
            lo, hi = min(branches.axes[axis]), max(branches.axes[axis])
            v = float(rng.uniform(lo, hi)) if hi > lo else lo
            if axis == "moderator_temperature":
                v = float(min(max(5.0 * round(v / 5.0), lo), hi))
            elif axis == "fuel_temperature":
                v = float(round(v))  # whole kelvin: one Doppler table each
            values.append(round(v, 4))
        fuel, mod_t, mod_rho, ppm = values
        out.append(State(f"{branches.name}-check{n}", fuel, mod_t, mod_rho, ppm))
    return out


def run_states(
    model: Model,
    name: str,
    states: list[State],
    solver: str = "infinite",
    options: dict | None = None,
    spectrum: str = "b1",
    cuts: tuple[float, ...] = (THERMAL_CUT,),
    path: Path | None = None,
    progress: Callable[[str], None] | None = None,
) -> list[dict]:
    """The lattice calculation at arbitrary states (same records as the
    matrix points, without an index); resumable like ``run_matrix``."""
    branches = model.branches[name]
    lattice = branches.lattice or next(iter(model.lattices))
    options = dict(options or {})
    head = header(model, branches, lattice, solver, options, spectrum, cuts)
    done: dict[str, dict] = {}
    if path is not None and path.exists():
        old = json.loads(path.read_text())
        if _same_run(old, head):
            done = {r["name"]: r for r in old["points"]}
    for state in states:
        if state.name in done:
            continue
        t0 = time.perf_counter()
        k, fine = lattice_constants(model, lattice, state, solver, options)
        fg = few_group(fine, spectrum, cuts)
        record = _point_record((), state, k, fg, time.perf_counter() - t0)
        del record["index"]
        done[state.name] = {"name": state.name, **record}
        _drop_shaped_libraries(state.name)
        if path is not None:
            path.parent.mkdir(parents=True, exist_ok=True)
            data = {**head, "boundaries_eV": list(fg.boundaries),
                    "points": list(done.values())}  # fmt: skip
            tmp = path.with_suffix(".tmp")
            tmp.write_text(json.dumps(data, indent=1))
            tmp.replace(path)
        if progress is not None:
            progress(f"{state.name}: k-inf {k:.5f} ({time.perf_counter() - t0:.0f} s)")
    return [done[s.name] for s in states]


def _indexed(branches: Branches) -> Iterator[tuple[tuple[int, ...], State]]:
    for state in branches.states():
        label = state.name.removeprefix(f"{branches.name}-")
        yield tuple(int(i) for i in label.split("-")), state


def _drop_shaped_libraries(state_name: str) -> None:
    """Delete the cached libraries whose flux shape belongs to this point
    (their file names end in ``-<state name>-<key>.npz``)."""
    directory = EndfLibrary().root / "groups"
    for f in directory.glob(f"*-{state_name}-*.npz"):
        f.unlink()


# ---------------------------------------------------------------- coefficients
@dataclass
class Coefficient:
    axis: str
    low: float
    high: float
    value: float  # d rho / d x, pcm per unit of the axis (central if possible)


def coefficients(table: Table) -> list[Coefficient]:
    """Reactivity coefficients at the base point: d rho / dx along each axis,
    rho = (k - 1) / k with the lattice k-infinity, from the neighbours of
    the base (central difference where it has two)."""
    k = table.array("k_infinity")
    rho = 1e5 * (k - 1.0) / k  # pcm
    base = [
        int(np.argmin(np.abs(np.asarray(table.axes[a]) - table.header["base"][a])))
        for a in BRANCH_AXES
    ]
    out = []
    for n, axis in enumerate(BRANCH_AXES):
        values = table.axes[axis]
        if len(values) < 2:
            continue
        lo, hi = max(base[n] - 1, 0), min(base[n] + 1, len(values) - 1)
        i_lo, i_hi = list(base), list(base)
        i_lo[n], i_hi[n] = lo, hi
        slope = (rho[tuple(i_hi)] - rho[tuple(i_lo)]) / (values[hi] - values[lo])
        out.append(Coefficient(axis, values[lo], values[hi], float(slope)))
    return out
