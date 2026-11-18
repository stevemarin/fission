"""Full-energy-range pointwise cross sections, linear-linear on one shared grid.

This is what a lattice code consumes: for one nuclide, every important
reaction tabulated from 1e-5 eV to the top of the evaluation, such that linear
interpolation reproduces the evaluated data to a tolerance.

    sigma(E) = MF3(E)                          everywhere
             + resonance part (MF2)             in the resolved range, for
                                                elastic, capture and fission
                                                (and therefore the total)

In the unresolved range all our evaluations set LSSF = 1: MF3 already holds
the infinitely dilute average there (self-shielding comes later).

Three practical issues:

* Redundancy. ENDF stores sums as well as their parts: MT1 total, MT3
  nonelastic, MT4 = sum of inelastic levels 51-91, MT18 = sum of fission
  chances 19-21, 38, MT101 = sum of absorption channels, MT103-107 = sums of
  their level series. We keep MT1 and one non-redundant set of partials, so
  MT1 = sum(partials) can be checked.
* Discontinuities. ENDF marks a jump by repeating an energy. A lin-lin table
  cannot hold two values at one energy, so (like NJOY) the pair is "shaded"
  into E(1 - 1e-7) and E(1 + 1e-7).
* Linearization. MF3 tables may use log or histogram interpolation laws (H-1
  capture is log-log, i.e. 1/v); adaptive bisection adds the points that make
  linear interpolation accurate.
"""

from __future__ import annotations

import bisect
import math
from collections.abc import Callable
from dataclasses import dataclass

from .evaluation import Evaluation
from .mlbw import mlbw
from .reconstruct import adaptive_grid, reconstruction_seeds
from .reich_moore import reich_moore
from .resonances import ResolvedRange
from .slbw import slbw

SHADE = 1e-7  # relative offset used to split a discontinuity
TOTAL, ELASTIC, FISSION, CAPTURE = 1, 2, 18, 102

# summary MTs kept, and the component MTs they make redundant
_SUMMARIES = {
    4: range(51, 92),
    18: (19, 20, 21, 38),
    103: range(600, 650),
    104: range(650, 700),
    105: range(700, 750),
    106: range(750, 800),
    107: range(800, 850),
}
_ALWAYS_REDUNDANT = {1, 3, 101}  # total, nonelastic, disappearance


def partial_mts(ev: Evaluation) -> list[int]:
    """A non-redundant set of MF3 reactions whose sum should equal MT1.

    MT5 ("anything") is a genuine partial: everything not given explicitly,
    used by high-energy evaluations above ~20 MeV.
    """
    present = set(ev.reactions)
    keep = set()
    for mt in present:
        if mt in _ALWAYS_REDUNDANT or mt > 891:  # >891: production/derived data
            continue
        covered = any(
            mt in parts and summary in present for summary, parts in _SUMMARIES.items()
        )
        if not covered:
            keep.add(mt)
    return sorted(keep)


def evaluated_total_inconsistency(ev: Evaluation, energies: list[float]) -> float:
    """Largest |MT1 - sum of partials| / MT1 of the raw MF3 data at `energies`.

    Evaluations are not perfectly self-consistent (H-1: 0.6% at 1e-5 eV), which
    is why `pointwise` rebuilds the total from the partials.
    """
    total = ev.cross_section(TOTAL)
    parts = [ev.cross_section(mt) for mt in partial_mts(ev)]
    worst = 0.0
    for e in energies:
        t = total(e)
        if t > 0:
            worst = max(worst, abs(t - sum(p(e) for p in parts)) / t)
    return worst


@dataclass
class PointwiseTable:
    """Cross sections (barns) on one energy grid (eV), linear-linear between points."""

    energies: list[float]
    xs: dict[int, list[float]]  # MT -> values on the grid

    def __len__(self) -> int:
        return len(self.energies)

    @property
    def mts(self) -> list[int]:
        return sorted(self.xs)

    def __call__(self, mt: int, energy: float) -> float:
        e, y = self.energies, self.xs[mt]
        if not e[0] <= energy <= e[-1]:
            return 0.0
        i = min(bisect.bisect_right(e, energy) - 1, len(e) - 2)
        t = (energy - e[i]) / (e[i + 1] - e[i])
        return y[i] + t * (y[i + 1] - y[i])

    def integral_over_e(self, mt: int, lo: float, hi: float) -> float:
        """Integral of sigma(E) dE / E over [lo, hi], exact for lin-lin data."""
        e, y = self.energies, self.xs[mt]
        total = 0.0
        for i in range(len(e) - 1):
            a, b = max(e[i], lo), min(e[i + 1], hi)
            if b <= a:
                continue
            slope = (y[i + 1] - y[i]) / (e[i + 1] - e[i])
            ya = y[i] + slope * (a - e[i])
            total += (ya - slope * a) * math.log(b / a) + slope * (b - a)
        return total


def shaded_points(x: list[float]) -> tuple[set[float], set[float]]:
    """Grid points of an ENDF table with every repeated energy split in two.

    Returns (points, jumps): the shaded points, and the jump energies that
    must not appear on the final grid themselves.
    """
    points: set[float] = set()
    jumps: set[float] = set()
    for i, e in enumerate(x):
        repeated = (i > 0 and x[i - 1] == e) or (i + 1 < len(x) and x[i + 1] == e)
        if repeated:
            jumps.add(e)
            points.update((e * (1 - SHADE), e * (1 + SHADE)))
        else:
            points.add(e)
    return points, jumps


def _resonance_evaluator(
    rr: ResolvedRange | None, use_rust: bool, python_f: Callable[[float], tuple]
) -> Callable[[float], tuple[float, float, float]] | None:
    """Resonance-only (elastic, capture, fission), in Rust when available."""
    if rr is None:
        return None
    if use_rust and rr.lrf in (2, 3):
        from .fast import resonance_model

        model = resonance_model(rr)
        return lambda e: model(e)
    return python_f


def pointwise(
    ev: Evaluation,
    mts: list[int] | None = None,
    lo: float | None = None,
    hi: float | None = None,
    tol: float = 1e-3,
    atol: float = 1e-5,
    use_rust: bool = True,
) -> PointwiseTable:
    """Linear-linear table of MT1 and the non-redundant partials.

    MT1 is the sum of the partials. Extra MTs (e.g. redundant sums such as 4 or
    18 components) can be requested with `mts`. lo/hi restrict the energy range
    (default: 1e-5 eV to the evaluation's EMAX).
    """
    partials = partial_mts(ev)
    mts = sorted({TOTAL, *partials, *(mts or [])})
    tables = {mt: ev.cross_section(mt) for mt in mts if mt in ev.reactions}
    lo = 1e-5 if lo is None else lo
    hi = ev.header.emax if hi is None else hi

    rr = None
    if ev.header.lrp == 1:
        (iso,) = ev.resonances
        rr = iso.resolved
    resonance = None
    if rr is not None:
        formalism = {1: slbw, 2: mlbw, 3: reich_moore}[rr.lrf]

        def python_resonance(e: float) -> tuple[float, float, float]:
            x = formalism(rr, e)
            return x.elastic, x.capture, x.fission

        resonance = _resonance_evaluator(rr, use_rust, python_resonance)
    slots = {mt: i for i, mt in enumerate(mts)}

    def f(e: float) -> tuple[float, ...]:
        values = [tables[mt](e) if mt in tables and mt != TOTAL else 0.0 for mt in mts]
        if resonance is not None and rr is not None and rr.el <= e <= rr.eh:
            el, ca, fi = resonance(e)
            for mt, v in ((ELASTIC, el), (CAPTURE, ca), (FISSION, fi)):
                if mt in slots:
                    values[slots[mt]] += v
        # the total is rebuilt from the partials (as NJOY does), so it is
        # consistent by construction even where the evaluation's MT1 is not
        values[slots[TOTAL]] = math.fsum(values[slots[mt]] for mt in partials)
        return tuple(values)

    # seeds: every MF3 breakpoint (jumps shaded), resonance seeds, decades
    seeds: set[float] = {lo, hi}
    jumps: set[float] = set()
    for t in tables.values():
        points, j = shaded_points(t.table.x)
        seeds |= points
        jumps |= j
    if rr is not None:
        seeds.update(reconstruction_seeds(ev, max(lo, rr.el), min(hi, rr.eh)))
        for edge in (rr.el, rr.eh):  # the resonance part switches on/off here
            if lo < edge < hi:
                seeds |= {edge * (1 - SHADE), edge * (1 + SHADE)}
                jumps.add(edge)
    seeds = {e for e in seeds if lo <= e <= hi and e not in jumps}
    # a shaded pair straddles a jump: refining it would only put the jump back
    shaded_pairs = {(e * (1 - SHADE), e * (1 + SHADE)) for e in jumps}
    if use_rust and (rr is None or rr.lrf in (2, 3)):
        # the same refinement, entirely in Rust and on all cores
        from .. import _rust
        from .fast import resonance_model

        grid, columns = _rust.pointwise_table(
            resonance_model(rr) if rr is not None else None,
            rr.el if rr is not None else 0.0,
            rr.eh if rr is not None else 0.0,
            [slots.get(mt) for mt in (ELASTIC, CAPTURE, FISSION)],
            [tables[mt].table if mt in tables and mt != TOTAL else None for mt in mts],
            slots[TOTAL],
            [slots[mt] for mt in partials],
            sorted(seeds),
            sorted(shaded_pairs),
            tol,
            atol,
        )
        return PointwiseTable(
            energies=grid.tolist(),
            xs={mt: col.tolist() for mt, col in zip(mts, columns, strict=True)},
        )
    energies, values = adaptive_grid(
        f, sorted(seeds), tol, atol, unsplittable=shaded_pairs
    )
    return PointwiseTable(
        energies=energies, xs={mt: [v[slots[mt]] for v in values] for mt in mts}
    )
