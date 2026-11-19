"""Recoverable energy per fission, kappa(E), from ENDF MF1/MT458.

MT458 gives nine components of the energy released in fission [ENDF-102
sec. 1.5]: fragment kinetic energy EFR, prompt and delayed neutrons ENP/END,
prompt and delayed gammas EGP/EGD, delayed betas EB, neutrinos ENU, ER = ET -
ENU and the total ET. ER is what stays in the reactor; kappa Sigma_f is the
power density the nodal code normalizes to.

Three layouts (ENDF-102 sec. 1.5.1): thermal values (LFC=0, NPLY=0), a
polynomial in E (LFC=0, NPLY>0, coefficients in powers of eV), or thermal
values plus TAB1 records for some components (LFC=1). For LFC=1 the energy
dependence of ER is the sum of the changes of the tabulated components
other than ENU.

(Kept outside ``fission.endf`` so that adding it does not invalidate the
cached pointwise tables, whose key includes that package's code.)
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from functools import cache

from fission.endf import Evaluation
from fission.endf.interpolation import Tabulated1D

COMPONENTS = ("EFR", "ENP", "END", "EGP", "EGD", "EB", "ENU", "ER", "ET")
ER, ENU = 7, 6  # 0-based positions in COMPONENTS


@dataclass(frozen=True)
class EnergyRelease:
    """ER(E) in eV: thermal values plus polynomial terms or tabulated changes."""

    thermal: tuple[float, ...]  # the nine c0 values (eV)
    polynomial: tuple[tuple[float, ...], ...] = ()  # c1, c2... per component
    tables: tuple[tuple[int, Tabulated1D], ...] = ()  # (IFC - 1, E_IFC(E))

    def component(self, index: int, energy: float) -> float:
        value = self.thermal[index]
        for power, row in enumerate(self.polynomial, start=1):
            value += row[index] * energy**power
        return value

    def __call__(self, energy: float) -> float:
        """Recoverable energy ER(E) (eV) for fission induced at `energy` (eV)."""
        if not self.tables:
            return self.component(ER, energy)
        change = math.fsum(
            table(energy) - self.thermal[i]
            for i, table in self.tables
            if i not in (ENU, ER, 8)
        )
        return self.thermal[ER] + change


def energy_release(ev: Evaluation) -> EnergyRelease:
    cur = ev.cursor(1, 458)
    head = cur.cont()
    lfc, nfc = head.l2, head.n2
    c, values = cur.list_record()
    npoly = c.l2
    pairs = [values[2 * i] for i in range(len(values) // 2)]  # drop uncertainties
    rows = [tuple(pairs[9 * p : 9 * p + 9]) for p in range(npoly + 1)]
    tables = []
    if lfc == 1:
        for _ in range(nfc):
            t_head, t = cur.tab1()
            tables.append((t_head.l2 - 1, Tabulated1D(t.x, t.y, t.breaks, t.laws)))
    return EnergyRelease(rows[0], tuple(rows[1:]), tuple(tables))


@cache
def recoverable_energy(nuclide: str) -> EnergyRelease:
    from fission.data import EndfLibrary

    return energy_release(Evaluation.from_file(EndfLibrary().neutron(nuclide)))
