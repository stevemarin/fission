"""MF5: secondary energy distributions f(E -> E').

A section describes the outgoing-energy spectrum of one reaction as a sum of
partial distributions [ENDF-102 §5]:

    f(E -> E') = sum_k p_k(E) f_k(E -> E'),       sum_k p_k(E) = 1,

each normalized to 1 over E'. The pin-cell nuclides use two laws:

    LF = 1  tabulated: f_k given as a table of E' at a list of incident
            energies (prompt fission spectra, MT18)
    LF = 5  general evaporation: f_k(E -> E') = g(E'/theta(E)), one shape g
            scaled by a temperature-like theta(E) (delayed neutrons, MT455,
            where each partial is one precursor family and p_k its abundance)
"""

from __future__ import annotations

import math
from bisect import bisect_right
from dataclasses import dataclass

from ..evaluation import Evaluation
from ..interpolation import HISTOGRAM, LIN_LIN, Tabulated1D
from .incident import IncidentGrid, interpolate_panels


def _table(data) -> Tabulated1D:
    return Tabulated1D(data.x, data.y, data.breaks, data.laws)


@dataclass(frozen=True)
class Panel:
    """A distribution of E' at one incident energy, histogram or lin-lin.

    Integrals and means are exact for these two laws.
    """

    x: tuple[float, ...]
    y: tuple[float, ...]
    law: int = LIN_LIN

    def __post_init__(self):
        if self.law not in (HISTOGRAM, LIN_LIN):
            raise NotImplementedError(f"secondary-energy law {self.law}")

    @property
    def range(self) -> tuple[float, float]:
        return self.x[0], self.x[-1]

    def __call__(self, e_out: float) -> float:
        x = self.x
        if not x[0] <= e_out <= x[-1]:
            return 0.0
        i = min(bisect_right(x, e_out) - 1, len(x) - 2)
        if self.law == HISTOGRAM or x[i + 1] == x[i]:
            return self.y[i]
        t = (e_out - x[i]) / (x[i + 1] - x[i])
        return self.y[i] + t * (self.y[i + 1] - self.y[i])

    def _segments(self):
        for i in range(len(self.x) - 1):
            x1, x2, y1 = self.x[i], self.x[i + 1], self.y[i]
            y2 = y1 if self.law == HISTOGRAM else self.y[i + 1]
            yield x1, x2, y1, y2

    def integral(self) -> float:
        return math.fsum(
            0.5 * (y1 + y2) * (x2 - x1) for x1, x2, y1, y2 in self._segments()
        )

    def first_moment(self) -> float:
        """Int E' f(E') dE' (exact for a linear segment)."""
        return math.fsum(
            (x2 - x1) * (y1 * (2 * x1 + x2) + y2 * (x1 + 2 * x2)) / 6
            for x1, x2, y1, y2 in self._segments()
        )

    def mean(self) -> float:
        return self.first_moment() / self.integral()


@dataclass(frozen=True)
class TabulatedLaw:
    """LF = 1: f(E -> E') tabulated at a grid of incident energies."""

    grid: IncidentGrid
    panels: tuple[Panel, ...]

    def __call__(self, energy: float, e_out: float) -> float:
        b = self.grid.bracket(energy)
        p1, p2 = self.panels[b.i], self.panels[min(b.i + 1, len(self.panels) - 1)]
        return interpolate_panels(p1, p1.range, p2, p2.range, b.t, b.unit_base, e_out)

    def mean(self, energy: float) -> float:
        """Mean E'. Both interpolations are linear in the panels, so the mean
        is the interpolated mean (in unit base, of the mean position x)."""
        b = self.grid.bracket(energy)
        p1, p2 = self.panels[b.i], self.panels[min(b.i + 1, len(self.panels) - 1)]
        if not b.unit_base:
            return (1 - b.t) * p1.mean() + b.t * p2.mean()
        (a1, z1), (a2, z2) = p1.range, p2.range
        x = (1 - b.t) * (p1.mean() - a1) / (z1 - a1) + b.t * (p2.mean() - a2) / (
            z2 - a2
        )
        lo, hi = a1 + b.t * (a2 - a1), z1 + b.t * (z2 - z1)
        return lo + x * (hi - lo)


@dataclass(frozen=True)
class EvaporationLaw:
    """LF = 5: f(E -> E') = g(E'/theta(E)) / theta(E), for 0 <= E' <= E - U.

    g(x) is renormalized to 1 over its own range, so that each partial
    distribution integrates to 1 exactly (the file's g is normalized only to
    the precision of its digits).
    """

    u: float
    theta: Tabulated1D
    g: Panel

    def __call__(self, energy: float, e_out: float) -> float:
        if e_out > energy - self.u:
            return 0.0
        th = self.theta(energy)
        return self.g(e_out / th) / (th * self.g.integral())

    def mean(self, energy: float) -> float:
        return self.theta(energy) * self.g.mean()


type Law = TabulatedLaw | EvaporationLaw


@dataclass(frozen=True)
class Partial:
    probability: Tabulated1D  # p_k(E)
    law: Law
    lf: int


@dataclass(frozen=True)
class EnergyDistribution:
    """All partial distributions of one MF5 section."""

    mt: int
    partials: tuple[Partial, ...]

    def __call__(self, energy: float, e_out: float) -> float:
        """f(E -> E') in 1/eV."""
        return math.fsum(
            p.probability(energy) * p.law(energy, e_out) for p in self.partials
        )

    def mean(self, energy: float) -> float:
        """Average secondary energy (eV)."""
        return math.fsum(
            p.probability(energy) * p.law.mean(energy) for p in self.partials
        )


def energy_distribution(ev: Evaluation, mt: int) -> EnergyDistribution:
    """Read MF5 section `mt` (LF = 1 and 5)."""
    cur = ev.cursor(5, mt)
    head = cur.cont()
    partials = []
    for _ in range(head.n1):
        c, p = cur.tab1()
        lf = c.l2
        law: Law
        if lf == 1:
            c2, breaks, laws = cur.tab2()
            energies, panels = [], []
            for _ in range(c2.n2):
                ci, t = cur.tab1()
                if len(set(t.laws)) != 1:
                    raise NotImplementedError("mixed secondary-energy interpolation")
                energies.append(ci.c2)
                panels.append(Panel(tuple(t.x), tuple(t.y), t.laws[0]))
            law = TabulatedLaw(
                IncidentGrid.from_tab2(energies, breaks, laws), tuple(panels)
            )
        elif lf == 5:
            _, theta = cur.tab1()
            _, g = cur.tab1()
            if len(set(g.laws)) != 1:
                raise NotImplementedError("mixed interpolation in g(x)")
            law = EvaporationLaw(
                c.c1, _table(theta), Panel(tuple(g.x), tuple(g.y), g.laws[0])
            )
        else:
            raise NotImplementedError(f"MF5 LF={lf}")
        partials.append(Partial(_table(p), law, lf))
    return EnergyDistribution(mt, tuple(partials))
