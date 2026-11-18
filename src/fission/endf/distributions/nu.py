"""Neutrons per fission: nu-bar, delayed-neutron groups, multiplicity distribution.

MF1 [ENDF-102 §1.2-1.4]:

    MT452  nu_t(E)  average total neutrons per fission
    MT455  nu_d(E)  average delayed neutrons, plus the decay constants
                    lambda_i of the precursor families
    MT456  nu_p(E)  average prompt neutrons;  nu_t = nu_p + nu_d

Each is a polynomial (LNU = 1) or a table (LNU = 2). The family abundances
and the delayed spectra live in MF5 MT455 (one partial per family).

MF6 MT18 may also hold, with JP > 0, the probability P(nu, E) that a fission
emits exactly nu prompt neutrons [ENDF-102 §6.3.9]. The file stores
P~(nu, E) = P(nu, E) / nu_p(E) as the "yield" of the nu-th product, so that

    sum_nu nu * P~(nu, E) = 1   and   sum_nu P~(nu, E) = 1 / nu_p(E).
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from ..evaluation import Evaluation
from ..interpolation import Tabulated1D
from .energy import EnergyDistribution, energy_distribution


@dataclass(frozen=True)
class NuBar:
    """nu(E) from one MF1 section: a polynomial or a table."""

    mt: int
    coefficients: tuple[float, ...] | None = None  # LNU = 1: nu = sum C_n E^(n-1)
    table: Tabulated1D | None = None  # LNU = 2

    def __call__(self, energy: float) -> float:
        if self.table is not None:
            return self.table(energy)
        assert self.coefficients is not None
        return math.fsum(c * energy**n for n, c in enumerate(self.coefficients))


def _read_nu(cur, lnu: int, mt: int) -> NuBar:
    if lnu == 1:
        _, coefficients = cur.list_record()
        return NuBar(mt, coefficients=tuple(coefficients))
    _, t = cur.tab1()
    return NuBar(mt, table=Tabulated1D(t.x, t.y, t.breaks, t.laws))


def nu_bar(ev: Evaluation, mt: int = 452) -> NuBar:
    """nu_t (452), nu_d (455) or nu_p (456)."""
    cur = ev.cursor(1, mt)
    head = cur.cont()
    lnu = head.l2
    if mt == 455:
        if head.l1 != 0:
            raise NotImplementedError("energy-dependent delayed constants (LDG=1)")
        cur.list_record()  # decay constants
    return _read_nu(cur, lnu, mt)


@dataclass(frozen=True)
class DelayedGroups:
    """Delayed-neutron precursor families."""

    decay_constants: tuple[float, ...]  # lambda_i (1/s), from MF1 MT455
    spectra: EnergyDistribution  # MF5 MT455: p_k(E) = abundance, f_k = spectrum
    nu_d: NuBar

    @property
    def half_lives(self) -> tuple[float, ...]:
        return tuple(math.log(2) / lam for lam in self.decay_constants)

    def abundances(self, energy: float) -> tuple[float, ...]:
        return tuple(p.probability(energy) for p in self.spectra.partials)

    def mean_energies(self, energy: float) -> tuple[float, ...]:
        return tuple(p.law.mean(energy) for p in self.spectra.partials)


def delayed_groups(ev: Evaluation) -> DelayedGroups:
    cur = ev.cursor(1, 455)
    head = cur.cont()
    if head.l1 != 0:
        raise NotImplementedError("energy-dependent delayed constants (LDG=1)")
    _, lambdas = cur.list_record()
    spectra = energy_distribution(ev, 455)
    if len(spectra.partials) != len(lambdas):
        raise ValueError("MF1 and MF5 disagree on the number of delayed families")
    return DelayedGroups(tuple(lambdas), spectra, nu_bar(ev, 455))


@dataclass(frozen=True)
class MultiplicityDistribution:
    """P~(nu, E) from MF6 MT18 (JP > 0), for nu = 0..nu_max."""

    tables: tuple[Tabulated1D, ...]  # P~(nu, E), index nu

    @property
    def nu_max(self) -> int:
        return len(self.tables) - 1

    def scaled(self, energy: float) -> list[float]:
        """P~(nu, E) as stored."""
        return [t(energy) for t in self.tables]

    def probabilities(self, energy: float) -> list[float]:
        """P(nu, E), normalized so that sum_nu P = 1."""
        s = self.scaled(energy)
        total = math.fsum(s)
        return [v / total for v in s]

    def mean(self, energy: float) -> float:
        """nu_p implied by the distribution: 1 / sum_nu P~(nu, E)."""
        return 1.0 / math.fsum(self.scaled(energy))


def multiplicity_distribution(ev: Evaluation) -> MultiplicityDistribution | None:
    """The prompt-neutron multiplicity distribution, if the file has one."""
    if (6, 18) not in ev.sections:
        return None
    cur = ev.cursor(6, 18)
    head = cur.cont()
    if head.l1 % 10 == 0:  # JPN = 0: no neutron multiplicity distribution
        return None
    tables = []
    for _ in range(head.n1):
        c, y = cur.tab1()
        if int(c.c1) != 1:
            break  # neutrons come first; photons (ZAP = 0) follow
        if c.l2 > 0:
            raise NotImplementedError("multiplicity-dependent spectra (JPN = 2)")
        tables.append(Tabulated1D(y.x, y.y, y.breaks, y.laws))
    return MultiplicityDistribution(tuple(tables))
