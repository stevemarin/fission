"""Interpolating distributions between incident energies.

A distribution such as a fission spectrum is given at a list of incident
energies E_1 < E_2 < ... (a TAB2 record). At an energy in between, ENDF
specifies how to combine the two neighbouring distributions with the TAB2
interpolation code [ENDF-102 §0.5.2]:

    INT = 1, 2       ordinary ("Cartesian") interpolation at fixed E'
    INT = 21, 22     unit-base interpolation: map each distribution onto
                     x = (E' - E'_low) / (E'_high - E'_low) in [0, 1],
                     interpolate at fixed x, map back to the interpolated
                     range, and multiply by the Jacobian dx/dE' so that the
                     result stays normalized.
    INT = 11, 12     "method of corresponding points". ENDF-102 (§0.5.2)
                     splits each distribution into bins of equal integral
                     and interpolates corresponding bins by unit base; it
                     does not fix the number of bins. With one bin, the
                     whole range, this is unit-base interpolation, and that
                     is what is done here. (Silicon's MF6 continua use
                     INT = 12; their neighbouring panels have different
                     numbers of points, so point-by-point correspondence is
                     not possible.)

Unit base matters when the range of E' moves with E (e.g. inelastic
scattering, E' <= E - threshold): Cartesian interpolation would then mix a
distribution ending at E'_1 with one ending at E'_2 and produce a spurious
two-step shape.
"""

from __future__ import annotations

import math
from bisect import bisect_left, bisect_right
from collections.abc import Callable, Sequence
from dataclasses import dataclass


@dataclass(frozen=True)
class Bracket:
    """Where an incident energy falls: between grid points i and i+1."""

    i: int
    t: float  # 0 at E_i, 1 at E_{i+1}
    unit_base: bool


@dataclass(frozen=True)
class IncidentGrid:
    """The incident energies of a TAB2 record and their interpolation codes."""

    energies: tuple[float, ...]
    breaks: tuple[int, ...]  # NBT, 1-based
    laws: tuple[int, ...]  # INT

    @classmethod
    def from_tab2(
        cls, energies: Sequence[float], breaks: Sequence[int], laws: Sequence[int]
    ):
        return cls(tuple(energies), tuple(breaks), tuple(laws))

    def bracket(self, energy: float) -> Bracket:
        """Neighbours of `energy`; outside the grid, the end distribution is used."""
        es = self.energies
        if len(es) == 1 or energy <= es[0]:
            return Bracket(0, 0.0, False)
        if energy >= es[-1]:
            return Bracket(len(es) - 2, 1.0, False)
        i = bisect_right(es, energy) - 1  # a repeated energy: take the upper side
        law = self.laws[bisect_left(self.breaks, i + 2)]
        base = law % 10 if law > 10 else law  # 11-15 and 21-25 follow laws 1-5
        if base == 1:
            t = 0.0
        elif base == 2:
            t = (energy - es[i]) / (es[i + 1] - es[i])
        elif base == 3:
            t = math.log(energy / es[i]) / math.log(es[i + 1] / es[i])
        else:
            raise NotImplementedError(f"incident-energy interpolation law {law}")
        return Bracket(i, t, law > 10)  # corresponding points: unit base (above)


def interpolate_panels(
    f1: Callable[[float], float],
    range1: tuple[float, float],
    f2: Callable[[float], float],
    range2: tuple[float, float],
    t: float,
    unit_base: bool,
    e_out: float,
) -> float:
    """Value at E' of the distribution interpolated a fraction t from panel 1 to 2.

    f1, f2 are the neighbouring distributions and range1, range2 their
    secondary-energy ranges (E'_low, E'_high).
    """
    if not unit_base:
        return (1.0 - t) * f1(e_out) + t * f2(e_out)
    lo = range1[0] + t * (range2[0] - range1[0])
    hi = range1[1] + t * (range2[1] - range1[1])
    if not lo <= e_out <= hi or hi <= lo:
        return 0.0
    x = (e_out - lo) / (hi - lo)
    w1, w2 = range1[1] - range1[0], range2[1] - range2[0]
    g1 = f1(range1[0] + x * w1) * w1  # unit-base density g(x) = f(E') dE'/dx
    g2 = f2(range2[0] + x * w2) * w2
    return ((1.0 - t) * g1 + t * g2) / (hi - lo)


def interpolate_parameter(
    f1: Callable[[float], float],
    range1: tuple[float, float],
    f2: Callable[[float], float],
    range2: tuple[float, float],
    t: float,
    unit_base: bool,
    e_out: float,
) -> float:
    """Like `interpolate_panels`, for a quantity that is not a density in E'
    (e.g. Kalbach's pre-compound fraction r): no Jacobian in unit base."""
    if not unit_base:
        return (1.0 - t) * f1(e_out) + t * f2(e_out)
    lo = range1[0] + t * (range2[0] - range1[0])
    hi = range1[1] + t * (range2[1] - range1[1])
    if not lo <= e_out <= hi or hi <= lo:
        return 0.0
    x = (e_out - lo) / (hi - lo)
    v1 = f1(range1[0] + x * (range1[1] - range1[0]))
    v2 = f2(range2[0] + x * (range2[1] - range2[0]))
    return (1.0 - t) * v1 + t * v2
