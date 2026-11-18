"""One ENDF-6 evaluation: split into sections, decode the header and cross sections.

An evaluation file holds one material (MAT). Its data is organized as
files (MF) of sections (MT). The ones that matter so far:

    MF1 / MT451   descriptive header: what this is, who made it, flags, directory
    MF2 / MT151   resonance parameters                       (next step)
    MF3 / MT*     "smooth" cross sections sigma(E), one section per reaction

Important for LWR physics: when the header flag LRP = 1, the cross sections
in the resolved resonance range are NOT in MF3. MF3 then holds only a smooth
"background" there (often zero), and the resonances must be reconstructed from
the MF2 parameters and added on top. So MF3 alone gives the right answer for
hydrogen or oxygen, but not for uranium in its resonance range.
"""

from __future__ import annotations

from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path

from .interpolation import Tabulated1D
from .reactions import mt_name
from .records import Cursor, control

SUBLIBRARIES = {0: "photonuclear", 10: "incident neutron", 12: "thermal neutron scattering",
                3: "photo-atomic", 4: "radioactive decay", 11: "neutron fission yields"}
LRP_MEANING = {
    -1: "no MF2",
    0: "no resonance parameters (MF2 holds only a scattering radius)",
    1: "resonance parameters given; add their contribution to MF3",
    2: "resonance parameters given for information only; MF3 is complete",
}


@dataclass(frozen=True)
class Header:
    """The MF1/MT451 descriptive header."""
    za: int          # 1000 * Z + A
    awr: float       # mass in units of the neutron mass
    lrp: int         # resonance-parameter flag (see LRP_MEANING)
    lfi: int         # 1 if the material is fissionable
    nlib: int        # library identifier (0 = ENDF/B)
    nmod: int        # modification number
    elis: float      # excitation energy of the target state (eV)
    sta: float       # 0 stable, 1 radioactive target
    lis: int         # excited-state number of the target
    liso: int        # isomeric-state number of the target
    nfor: int        # format version (6 = ENDF-6)
    awi: float       # projectile mass in neutron masses (1 for neutrons)
    emax: float      # upper energy limit of the evaluation (eV)
    lrel: int        # release number
    nsub: int        # sublibrary (10 incident neutron, 12 thermal scattering, ...)
    nver: int        # library version (8 for ENDF/B-VIII)
    temp: float      # target temperature (K); 0 for incident-neutron data
    ldrv: int        # derived-evaluation flag
    description: tuple[str, ...]                  # the free-text lines
    directory: tuple[tuple[int, int, int, int], ...]  # (MF, MT, lines, mod) per section

    @property
    def z(self) -> int:
        return self.za // 1000

    @property
    def a(self) -> int:
        return self.za % 1000

    @property
    def zsymam(self) -> str:
        """Nuclide label, e.g. ' 92-U -238'."""
        return self.description[0][0:11].strip()

    @property
    def lab(self) -> str:
        return self.description[0][11:22].strip()

    @property
    def edate(self) -> str:
        return self.description[0][22:32].strip()

    @property
    def authors(self) -> str:
        return self.description[0][33:66].strip()

    @property
    def sublibrary(self) -> str:
        return SUBLIBRARIES.get(self.nsub, f"NSUB={self.nsub}")


@dataclass(frozen=True)
class CrossSection:
    """An MF3 section: sigma(E) in barns, E in eV."""
    mt: int
    qm: float        # mass-difference Q value (eV)
    qi: float        # reaction Q value for the lowest state reached (eV)
    lr: int          # complex breakup flag
    table: Tabulated1D

    @property
    def name(self) -> str:
        return mt_name(self.mt)

    @property
    def threshold(self) -> float:
        return self.table.x[0]

    def __call__(self, energy: float) -> float:
        return self.table(energy)


class Evaluation:
    """An ENDF-6 evaluation file, split into its (MF, MT) sections."""

    def __init__(self, sections: dict[tuple[int, int], list[str]], mat: int):
        self.sections = sections
        self.mat = mat

    @classmethod
    def from_file(cls, path: str | Path) -> Evaluation:
        with open(path) as f:
            return cls.from_lines(f)

    @classmethod
    def from_lines(cls, lines) -> Evaluation:
        """Group lines into sections by (MF, MT).

        End-of-section (SEND, MT=0), end-of-file (FEND, MF=0), end-of-material
        (MEND, MAT=0) and the tape label line carry zeros in those fields, so
        keeping only MF > 0 and MT > 0 drops them all.
        """
        sections: dict[tuple[int, int], list[str]] = {}
        mats = set()
        for line in lines:
            mat, mf, mt = control(line)
            if mat > 0 and mf > 0 and mt > 0:
                sections.setdefault((mf, mt), []).append(line)
                mats.add(mat)
        if len(mats) != 1:
            raise ValueError(f"expected exactly one material, found MAT numbers {sorted(mats)}")
        return cls(sections, mats.pop())

    def cursor(self, mf: int, mt: int) -> Cursor:
        if (mf, mt) not in self.sections:
            raise KeyError(f"no MF{mf}/MT{mt} section in MAT {self.mat}")
        return Cursor(self.sections[(mf, mt)])

    # --------------------------------------------------------------- MF1/MT451
    @property
    def header(self) -> Header:
        c = self.cursor(1, 451)
        h1, h2, h3, h4 = c.head(), c.cont(), c.cont(), c.cont()
        description = tuple(c.text() for _ in range(h4.n1))         # NWD lines
        directory = tuple((d.l1, d.l2, d.n1, d.n2) for d in (c.cont() for _ in range(h4.n2)))
        return Header(
            za=round(h1.c1), awr=h1.c2, lrp=h1.l1, lfi=h1.l2, nlib=h1.n1, nmod=h1.n2,
            elis=h2.c1, sta=h2.c2, lis=h2.l1, liso=h2.l2, nfor=h2.n2,
            awi=h3.c1, emax=h3.c2, lrel=h3.l1, nsub=h3.n1, nver=h3.n2,
            temp=h4.c1, ldrv=h4.l1, description=description, directory=directory,
        )

    # ---------------------------------------------------------------------- MF3
    @property
    def reactions(self) -> list[int]:
        """MT numbers that have an MF3 cross section."""
        return sorted(mt for mf, mt in self.sections if mf == 3)

    def cross_section(self, mt: int) -> CrossSection:
        c = self.cursor(3, mt)
        c.head()  # ZA, AWR (repeated from the header)
        cont, tab = c.tab1()
        return CrossSection(mt=mt, qm=cont.c1, qi=cont.c2, lr=cont.l2,
                            table=Tabulated1D(tab.x, tab.y, tab.breaks, tab.laws))

    def cross_sections(self) -> Iterator[CrossSection]:
        for mt in self.reactions:
            yield self.cross_section(mt)
