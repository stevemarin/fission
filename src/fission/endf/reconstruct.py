"""Pointwise (linearly interpolable) cross sections in the resolved resonance range.

A transport code wants sigma(E) as a table it can interpolate linearly. Resonances
are narrow (U-238 at 6.67 eV is 0.024 eV wide) while the range spans 1e-5 eV to
20 keV, so a uniform grid would be absurd. Instead, as NJOY's RECONR does:

1. seed a grid: range ends, every resonance energy and points at +-G/2 and +-3G
   around it, the MF3 background grid, and a few points per decade;
2. for every interval, evaluate sigma at the midpoint and compare with linear
   interpolation from the ends; if any reaction misses by more than
   tol * |sigma| + atol, split the interval and repeat (depth first).

The result is a table on which linear interpolation reproduces the Reich-Moore
cross sections to about tol, with points only where sigma curves.

Within the resolved range the tabulated value is resonance + MF3 background,
since for LRP = 1 the file splits the cross section into those two parts.
"""

from __future__ import annotations

import bisect
import functools
import math
from collections.abc import Callable, Collection, Sequence
from dataclasses import dataclass

from .evaluation import Evaluation
from .mlbw import mlbw
from .reich_moore import reich_moore
from .resonances import ResolvedRange, RMLRange
from .rml import RMatrixLimited, levels_of
from .slbw import ResonanceXS, single_resonance, slbw

REACTIONS = ("elastic", "capture", "fission")
BACKGROUND_MT = {"elastic": 2, "capture": 102, "fission": 18}

Values = tuple[float, ...]  # one value per reaction, e.g. (elastic, capture, fission)


@dataclass
class PointwiseXS:
    """Cross sections on a shared energy grid, linear-linear between points (barns, eV)."""

    energies: list[float]
    elastic: list[float]
    capture: list[float]
    fission: list[float]

    def __len__(self) -> int:
        return len(self.energies)

    def reaction(self, name: str) -> list[float]:
        return getattr(self, name)

    def __call__(self, name: str, energy: float) -> float:
        """Linear interpolation of one reaction."""
        e, y = self.energies, self.reaction(name)
        if not e[0] <= energy <= e[-1]:
            raise ValueError(f"{energy} eV outside [{e[0]}, {e[-1]}]")
        i = min(bisect.bisect_right(e, energy) - 1, len(e) - 2)
        t = (energy - e[i]) / (e[i + 1] - e[i])
        return y[i] + t * (y[i + 1] - y[i])

    def integral_over_e(self, name: str, lo: float, hi: float) -> float:
        """Integral of sigma(E) dE / E over [lo, hi], exact for lin-lin data.

        For capture above 0.5 eV this is the resonance integral.
        """
        e, y = self.energies, self.reaction(name)
        total = 0.0
        for i in range(len(e) - 1):
            a, b = max(e[i], lo), min(e[i + 1], hi)
            if b <= a:
                continue
            slope = (y[i + 1] - y[i]) / (e[i + 1] - e[i])
            ya = y[i] + slope * (a - e[i])
            total += (ya - slope * a) * math.log(b / a) + slope * (b - a)
        return total


def seed_grid(
    rr: ResolvedRange | RMLRange, lo: float, hi: float, extra: Sequence[float] = ()
) -> list[float]:
    """Starting energies: resonance peaks and half-widths, decades, extras."""
    points = {lo, hi, *(e for e in extra if lo < e < hi)}
    decade = 10 ** math.floor(math.log10(lo))
    while decade < hi:  # 10 points per decade so smooth regions start sampled
        for i in range(10):
            e = decade * 10 ** (i / 10)
            if lo < e < hi:
                points.add(e)
        decade *= 10
    levels: list[tuple[float, float]] = []  # (energy, total width)
    if isinstance(rr, RMLRange):
        for sg in rr.groups:
            levels += levels_of(sg)
    else:
        for group in rr.groups:
            for res in group.resonances:
                if lo < res.energy < hi:
                    width = single_resonance(rr, group, res, res.energy).width
                    levels.append((res.energy, width))
    for energy, width in levels:
        if not lo < energy < hi:
            continue
        for offset in (-3.0, -0.5, 0.0, 0.5, 3.0):
            e = energy + offset * width
            if lo < e < hi:
                points.add(e)
    return sorted(points)


def _converged(
    left: Values, mid: Values, right: Values, tol: float, atol: float
) -> bool:
    return all(
        abs(m - 0.5 * (a + b)) <= tol * abs(m) + atol
        for a, m, b in zip(left, mid, right, strict=True)
    )


def adaptive_grid(
    f: Callable[[float], Values],
    seeds: Sequence[float],
    tol: float = 1e-3,
    atol: float = 1e-5,
    min_rel_width: float = 1e-9,
    unsplittable: Collection[tuple[float, float]] = (),
) -> tuple[list[float], list[Values]]:
    """Refine `seeds` until linear interpolation of f is within tolerance.

    Depth-first bisection: `left` is the last accepted point, `stack` holds
    the pending right ends, nearest on top. Intervals listed in `unsplittable`
    (pairs of adjacent seeds straddling a discontinuity) are kept as they are.
    """
    energies = [seeds[0]]
    values = [f(seeds[0])]
    for right_seed in seeds[1:]:
        stack = [(right_seed, f(right_seed))]
        while stack:
            e1, v1 = energies[-1], values[-1]
            e2, v2 = stack[-1]
            mid = 0.5 * (e1 + e2)
            if e2 - e1 > min_rel_width * e2 and (e1, e2) not in unsplittable:
                vm = f(mid)
                if not _converged(v1, vm, v2, tol, atol):
                    stack.append((mid, vm))
                    continue
            energies.append(e2)
            values.append(v2)
            stack.pop()
    return energies, values


def resolved_range_xs(
    ev: Evaluation,
    formalism: Callable[[ResolvedRange, float], ResonanceXS] | None = None,
) -> tuple[ResolvedRange | RMLRange, Callable[[float], Values]]:
    """sigma(E) = resonance part + MF3 background, for each reaction."""
    (iso,) = ev.resonances
    rr = iso.resolved
    if rr is None:
        raise ValueError("no resolved resonance range")
    resonance: Callable[[float], ResonanceXS]
    if isinstance(rr, RMLRange):
        resonance = RMatrixLimited(rr)
    else:
        model = formalism or {1: slbw, 2: mlbw, 3: reich_moore}[rr.lrf]
        resonance = functools.partial(model, rr)
    backgrounds = {
        name: ev.cross_section(mt) if mt in ev.reactions else None
        for name, mt in BACKGROUND_MT.items()
    }

    def f(energy: float) -> Values:
        xs = resonance(energy)
        out = []
        for name in REACTIONS:
            bg = backgrounds[name]
            out.append(getattr(xs, name) + (bg(energy) if bg is not None else 0.0))
        return out[0], out[1], out[2]

    return rr, f


def reconstruction_seeds(
    ev: Evaluation, lo: float | None = None, hi: float | None = None
) -> list[float]:
    """Seed grid for an evaluation's resolved range, clipped to [lo, hi]."""
    (iso,) = ev.resonances
    rr = iso.resolved
    if rr is None:
        raise ValueError("no resolved resonance range")
    lo = rr.el if lo is None else max(lo, rr.el)
    hi = rr.eh if hi is None else min(hi, rr.eh)
    background_grid: list[float] = []
    for mt in BACKGROUND_MT.values():
        if mt in ev.reactions:
            background_grid += ev.cross_section(mt).table.x
    return seed_grid(rr, lo, hi, extra=[0.0253, *background_grid])


def reconstruct(
    ev: Evaluation,
    lo: float | None = None,
    hi: float | None = None,
    tol: float = 1e-3,
    atol: float = 1e-5,
) -> PointwiseXS:
    """Pointwise resolved-range cross sections of an evaluation (default: whole range).

    Pure-Python reference; `fission.endf.fast.reconstruct` is the Rust port.
    """
    _, f = resolved_range_xs(ev)
    energies, values = adaptive_grid(f, reconstruction_seeds(ev, lo, hi), tol, atol)
    return PointwiseXS(
        energies=energies,
        elastic=[v[0] for v in values],
        capture=[v[1] for v in values],
        fission=[v[2] for v in values],
    )
