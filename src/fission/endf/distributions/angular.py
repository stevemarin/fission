"""MF4: angular distributions of the emitted neutron.

For each incident energy, f(mu, E) is the probability density of the cosine
mu of the scattering angle, normalized over [-1, 1], in the LAB or the
centre-of-mass (CM) frame [ENDF-102 §4]:

    LTT = 0  isotropic
    LTT = 1  Legendre coefficients a_l(E), l = 1..NL (a_0 = 1)
    LTT = 2  tabulated f(mu, E)
    LTT = 3  Legendre at low energies, tabulated above

Between incident energies the distributions are interpolated linearly
(INT = 2 in every pin-cell file), which is the same as interpolating the
Legendre moments linearly.
"""

from __future__ import annotations

from dataclasses import dataclass

from ..evaluation import Evaluation
from ..interpolation import Tabulated1D
from .incident import IncidentGrid
from .legendre import legendre_series, moments_of_linear

LAB, CM = 1, 2


@dataclass(frozen=True)
class LegendrePanel:
    """f(mu) = sum_l (2l+1)/2 a_l P_l(mu) with a_0 = 1."""

    coefficients: tuple[float, ...]  # a_1..a_NL

    def __call__(self, mu: float) -> float:
        return legendre_series((1.0, *self.coefficients), mu)

    def moments(self, lmax: int) -> list[float]:
        a = [1.0, *self.coefficients]
        return (a + [0.0] * (lmax + 1))[: lmax + 1]


@dataclass(frozen=True)
class TabulatedPanel:
    """f(mu) as a table, linear between points."""

    table: Tabulated1D

    def __call__(self, mu: float) -> float:
        return self.table(mu)

    def moments(self, lmax: int) -> list[float]:
        if set(self.table.laws) != {2}:
            raise NotImplementedError("non-linear interpolation in mu")
        return moments_of_linear(self.table.x, self.table.y, lmax)


type AngularPanel = LegendrePanel | TabulatedPanel


@dataclass(frozen=True)
class AngularDistribution:
    """One MF4 section (or an isotropic stand-in)."""

    mt: int
    frame: int  # LAB = 1, CM = 2
    grid: IncidentGrid | None  # None: isotropic at all energies
    panels: tuple[AngularPanel, ...]

    @property
    def isotropic(self) -> bool:
        return self.grid is None

    def moments(self, energy: float, lmax: int) -> list[float]:
        """a_0..a_lmax at an incident energy (a_0 = 1, a_1 = mean cosine)."""
        if self.grid is None:
            return [1.0] + [0.0] * lmax
        b = self.grid.bracket(energy)
        m1 = self.panels[b.i].moments(lmax)
        m2 = self.panels[min(b.i + 1, len(self.panels) - 1)].moments(lmax)
        return [(1 - b.t) * x + b.t * y for x, y in zip(m1, m2, strict=True)]

    def __call__(self, energy: float, mu: float) -> float:
        """f(mu, E)."""
        if self.grid is None:
            return 0.5
        b = self.grid.bracket(energy)
        f1 = self.panels[b.i](mu)
        f2 = self.panels[min(b.i + 1, len(self.panels) - 1)](mu)
        return (1 - b.t) * f1 + b.t * f2

    def mean_cosine(self, energy: float) -> float:
        """Average cosine in the frame of the data."""
        return self.moments(energy, 1)[1]


def _legendre_block(cur):
    c, breaks, laws = cur.tab2()
    energies, panels = [], []
    for _ in range(c.n2):
        ci, a = cur.list_record()
        energies.append(ci.c2)
        panels.append(LegendrePanel(tuple(a)))
    return energies, breaks, laws, panels


def _tabulated_block(cur):
    c, breaks, laws = cur.tab2()
    energies, panels = [], []
    for _ in range(c.n2):
        ci, t = cur.tab1()
        energies.append(ci.c2)
        panels.append(TabulatedPanel(Tabulated1D(t.x, t.y, t.breaks, t.laws)))
    return energies, breaks, laws, panels


def angular_distribution(ev: Evaluation, mt: int) -> AngularDistribution:
    """Read MF4 section `mt`."""
    cur = ev.cursor(4, mt)
    ltt = cur.cont().l2
    c = cur.cont()
    li, lct = c.l1, c.l2
    if ltt == 0 or li == 1:
        return AngularDistribution(mt, lct, None, ())
    if ltt in (1, 2):
        block = _legendre_block if ltt == 1 else _tabulated_block
        energies, breaks, laws, panels = block(cur)
    elif ltt == 3:
        e1, _, l1, p1 = _legendre_block(cur)
        e2, _, l2, p2 = _tabulated_block(cur)
        if set(l1) | set(l2) != {2}:
            raise NotImplementedError("LTT=3 with non-linear energy interpolation")
        # the boundary energy appears in both blocks: keep it twice so that the
        # tabulated panel is used from the boundary up (bracket takes the upper side)
        energies, panels = e1 + e2, p1 + p2
        breaks, laws = [len(energies)], [2]
    else:
        raise NotImplementedError(f"MF4 LTT={ltt}")
    grid = IncidentGrid.from_tab2(energies, breaks, laws)
    return AngularDistribution(mt, lct, grid, tuple(panels))
