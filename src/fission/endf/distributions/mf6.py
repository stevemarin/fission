"""MF6: energy-angle distributions of reaction products (neutrons only).

A section lists the products of one reaction [ENDF-102 §6]. Each has a yield
y(E) (neutrons per reaction: 1 for inelastic scattering, 2 for (n,2n)) and a
distribution law:

    LAW = 0   no distribution (only the yield)
    LAW < 0   distribution in MF4/MF5 (fission, with JP > 0; see nu.py)
    LAW = 1   continuum: f(mu, E -> E') at a grid of incident energies, with
              the angular part as Legendre coefficients (LANG = 1) or
              Kalbach-Mann systematics (LANG = 2)
    LAW = 2   discrete two-body: an angular distribution; E' follows from
              kinematics and the Q value of the level (MF3)

Products appear in the order particles, residuals, photons; only the leading
neutron products (ZAP = 1) are read.
"""

from __future__ import annotations

from dataclasses import dataclass

from ..evaluation import Evaluation
from ..interpolation import Tabulated1D
from .angular import AngularPanel, LegendrePanel, TabulatedPanel
from .energy import Panel
from .incident import IncidentGrid, interpolate_panels, interpolate_parameter
from .kinematics import kalbach_mann, kalbach_slope
from .legendre import legendre_series

LEGENDRE, KALBACH = 1, 2


@dataclass(frozen=True)
class ContinuumPanel:
    """One incident energy of LAW = 1: per outgoing energy E', the
    parameters b_0..b_NA (f_0 = emission probability density, then Legendre
    coefficients f_l or the Kalbach pre-compound fraction r)."""

    e_out: tuple[float, ...]
    b: tuple[tuple[float, ...], ...]  # b[i] = (b_0, ..., b_NA) at e_out[i]
    lep: int  # 1 histogram, 2 lin-lin in E'

    def component(self, j: int) -> Panel:
        return Panel(
            self.e_out,
            tuple(row[j] if j < len(row) else 0.0 for row in self.b),
            self.lep,
        )

    @property
    def range(self) -> tuple[float, float]:
        return self.e_out[0], self.e_out[-1]


@dataclass(frozen=True)
class Continuum:
    """LAW = 1."""

    lang: int
    grid: IncidentGrid
    panels: tuple[ContinuumPanel, ...]
    za: int  # target, for Kalbach's slope
    awr: float

    def _neighbours(self, energy: float):
        b = self.grid.bracket(energy)
        p1 = self.panels[b.i]
        p2 = self.panels[min(b.i + 1, len(self.panels) - 1)]
        return b, p1, p2

    def _density(self, energy: float, j: int, e_out: float) -> float:
        """f_0, or a Legendre component f_l: densities in E'."""
        b, p1, p2 = self._neighbours(energy)
        return interpolate_panels(
            p1.component(j),
            p1.range,
            p2.component(j),
            p2.range,
            b.t,
            b.unit_base,
            e_out,
        )

    def _parameter(self, energy: float, j: int, e_out: float) -> float:
        """Kalbach's r: a fraction, not a density."""
        b, p1, p2 = self._neighbours(energy)
        return interpolate_parameter(
            p1.component(j),
            p1.range,
            p2.component(j),
            p2.range,
            b.t,
            b.unit_base,
            e_out,
        )

    def energy_pdf(self, energy: float, e_out: float) -> float:
        """f_0(E -> E'): the outgoing-energy density, integrated over angle."""
        return self._density(energy, 0, e_out)

    def mean_energy(self, energy: float) -> float:
        """<E'>: exact at the grid energies, interpolated like the panels between."""
        b, p1, p2 = self._neighbours(energy)
        f1, f2 = p1.component(0), p2.component(0)
        if not b.unit_base:
            return (1 - b.t) * f1.mean() + b.t * f2.mean()
        (a1, z1), (a2, z2) = p1.range, p2.range
        x = (1 - b.t) * (f1.mean() - a1) / (z1 - a1) + b.t * (f2.mean() - a2) / (
            z2 - a2
        )
        lo, hi = a1 + b.t * (a2 - a1), z1 + b.t * (z2 - z1)
        return lo + x * (hi - lo)

    def __call__(self, energy: float, e_out: float, mu: float) -> float:
        """f(mu, E -> E') per eV per unit cosine, in the frame of the data."""
        if self.lang == LEGENDRE:
            na = max(len(row) for row in self._neighbours(energy)[1].b) - 1
            f = [self._density(energy, j, e_out) for j in range(na + 1)]
            return legendre_series(f, mu)
        f0 = self.energy_pdf(energy, e_out)
        if f0 == 0.0:
            return 0.0
        r = self._parameter(energy, 1, e_out)
        return f0 * kalbach_mann(mu, kalbach_slope(energy, e_out, self.za, self.awr), r)


@dataclass(frozen=True)
class TwoBody:
    """LAW = 2: CM angular distributions of a two-body reaction."""

    grid: IncidentGrid
    panels: tuple[AngularPanel, ...]

    def __call__(self, energy: float, mu: float) -> float:
        b = self.grid.bracket(energy)
        f1 = self.panels[b.i](mu)
        f2 = self.panels[min(b.i + 1, len(self.panels) - 1)](mu)
        return (1 - b.t) * f1 + b.t * f2

    def moments(self, energy: float, lmax: int) -> list[float]:
        b = self.grid.bracket(energy)
        m1 = self.panels[b.i].moments(lmax)
        m2 = self.panels[min(b.i + 1, len(self.panels) - 1)].moments(lmax)
        return [(1 - b.t) * x + b.t * y for x, y in zip(m1, m2, strict=True)]


type Distribution = Continuum | TwoBody | None


@dataclass(frozen=True)
class Product:
    zap: int
    awp: float
    lip: int
    law: int
    yield_: Tabulated1D  # neutrons per reaction (or P~(nu, E) for fission, JP > 0)
    distribution: Distribution


@dataclass(frozen=True)
class EnergyAngle:
    """The neutron products of one MF6 section."""

    mt: int
    za: int
    awr: float
    jp: int
    frame: int  # LCT: 1 LAB, 2 CM, 3 CM for light particles and LAB for heavy ones
    products: tuple[Product, ...]


def _continuum(cur, za: int, awr: float) -> Continuum:
    c, breaks, laws = cur.tab2()
    lang, lep = c.l1, c.l2
    if lang not in (LEGENDRE, KALBACH):
        raise NotImplementedError(f"MF6 LAW=1 LANG={lang}")
    energies, panels = [], []
    for _ in range(c.n2):
        ci, values = cur.list_record()
        nd, na, nep = ci.l1, ci.l2, ci.n2
        if nd:
            raise NotImplementedError("discrete lines in MF6 LAW=1 (ND > 0)")
        width = na + 2
        rows = [values[k * width : (k + 1) * width] for k in range(nep)]
        energies.append(ci.c2)
        panels.append(
            ContinuumPanel(
                tuple(r[0] for r in rows), tuple(tuple(r[1:]) for r in rows), lep
            )
        )
    return Continuum(
        lang, IncidentGrid.from_tab2(energies, breaks, laws), tuple(panels), za, awr
    )


def _two_body(cur) -> TwoBody:
    c, breaks, laws = cur.tab2()
    energies: list[float] = []
    panels: list[AngularPanel] = []
    for _ in range(c.n2):
        ci, values = cur.list_record()
        lang = ci.l1
        energies.append(ci.c2)
        if lang == 0:
            panels.append(LegendrePanel(tuple(values)))
        elif lang == 12:
            mu, p = values[0::2], values[1::2]
            panels.append(TabulatedPanel(Tabulated1D(mu, p)))
        else:
            raise NotImplementedError(f"MF6 LAW=2 LANG={lang}")
    return TwoBody(IncidentGrid.from_tab2(energies, breaks, laws), tuple(panels))


def energy_angle(ev: Evaluation, mt: int) -> EnergyAngle:
    """Read the neutron products of MF6 section `mt`."""
    cur = ev.cursor(6, mt)
    head = cur.cont()
    za, awr, jp, lct = int(head.c1), head.c2, head.l1, head.l2
    products = []
    for _ in range(head.n1):
        c, y = cur.tab1()
        zap, awp, lip, law = int(c.c1), c.c2, c.l1, c.l2
        if zap != 1:
            break
        dist: Distribution
        if law <= 0:
            dist = None
        elif law == 1:
            dist = _continuum(cur, za, awr)
        elif law == 2:
            dist = _two_body(cur)
        else:
            raise NotImplementedError(f"MF6 LAW={law} for neutrons")
        products.append(
            Product(zap, awp, lip, law, Tabulated1D(y.x, y.y, y.breaks, y.laws), dist)
        )
    return EnergyAngle(mt, za, awr, jp, lct, tuple(products))
