"""MF2 / MT151: resonance parameters (reading only, no cross sections yet).

Structure of the section::

    HEAD   ZA, AWR, -, -, NIS, -                    NIS isotopes (1 for our files)
    per isotope:
      CONT ZAI, ABN, -, LFW, NER, -                 NER energy ranges
      per energy range:
        CONT EL, EH, LRU, LRF, NRO, NAPS
             LRU  0 no resonances, 1 resolved, 2 unresolved
             LRF  formalism: 1 SLBW, 2 MLBW, 3 Reich-Moore, 7 R-Matrix Limited
                  (for unresolved ranges LRF says which widths depend on energy)
             NRO  1 if the scattering radius depends on energy (then a TAB1 follows)
             NAPS how to get the channel radius used in penetrabilities (see physics step)
        ... formalism-specific records ...

Units: energies and widths in eV, radii in 10^-12 cm (= 10 fm).

What a resonance is: the neutron and the target briefly form a compound
nucleus. It only exists at certain excitation energies; when the neutron energy
hits one, the cross section peaks. Each resonance is described by its energy
E_r, the compound-state spin J, and partial widths Gamma for each way the state
can decay (neutron, gamma, fission). Resonances at negative energy are bound
states just below the neutron separation energy; their tails shape the
cross section near thermal energies.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from .interpolation import Tabulated1D
from .records import Cursor

FORMALISMS = {
    1: "SLBW",
    2: "MLBW",
    3: "Reich-Moore",
    4: "Adler-Adler",
    7: "R-Matrix Limited",
}


@dataclass(frozen=True)
class Resonance:
    """One resolved resonance. Widths in eV.

    For Reich-Moore, fission goes through two channels. gfa and gfb are the
    partial widths of the two channels, but carry a sign: the sign of the
    channel's reduced-width amplitude, which decides whether neighbouring
    resonances interfere constructively or destructively. For SLBW/MLBW, gfa is
    the fission width and gfb is 0.
    """

    energy: float
    j: float  # compound-state spin J (some evaluations use the sign to pick
    #                   the channel spin; kept as given)
    gn: float  # neutron width
    gg: float  # radiation (capture) width
    gfa: float  # fission width (BW) / first fission channel amplitude (RM)
    gfb: float  # 0 (BW) / second fission channel amplitude (RM)

    @property
    def total_width(self) -> float:
        return self.gn + self.gg + abs(self.gfa) + abs(self.gfb)


@dataclass(frozen=True)
class LGroup:
    """All resonances with one orbital angular momentum l."""

    l: int
    awri: float  # mass ratio of this isotope (target / neutron)
    apl: float  # l-dependent scattering radius (RM), 0 means "use AP"
    resonances: tuple[Resonance, ...]
    lrx: int = 0  # BW only: 1 if a competitive (inelastic) width is present


@dataclass(frozen=True)
class ResolvedRange:
    el: float
    eh: float
    lrf: int
    naps: int
    spi: float  # target spin I
    ap: float  # scattering radius (10^-12 cm)
    lad: int  # RM only: 1 if parameters may also be used for angular distributions
    nlsc: int  # RM only: number of l values needed for convergence of potential scattering
    groups: tuple[LGroup, ...]
    ap_table: Tabulated1D | None = None  # energy-dependent radius when NRO = 1

    @property
    def formalism(self) -> str:
        return FORMALISMS.get(self.lrf, f"LRF={self.lrf}")

    @property
    def resonances(self) -> list[Resonance]:
        return [r for g in self.groups for r in g.resonances]


@dataclass(frozen=True)
class ParticlePair:
    """Two particles a channel can hold (LRF=7): e.g. gamma + compound
    nucleus (MT=102) or neutron + target in its ground state (MT=2)."""

    ma: float  # masses in neutron masses
    mb: float
    za: float  # charges
    zb: float
    ia: float  # spins (the sign, if any, is the parity)
    ib: float
    q: float  # Q-value, eV
    pnt: int  # 1 compute the penetrability, -1 do not, 0 default by MT
    shf: int  # 1 compute the shift factor
    mt: int  # reaction the pair stands for

    @property
    def is_gamma(self) -> bool:
        return self.mt == 102


@dataclass(frozen=True)
class Channel:
    """One channel of a spin group: a particle pair with orbital momentum l
    and channel spin s, and its two radii (10^-12 cm)."""

    pair: int  # index into RMLRange.pairs (0-based)
    l: int
    s: float
    bnd: float  # boundary condition (only used with the shift factor)
    ape: float  # effective radius: hard-sphere phase shift
    apt: float  # true radius: penetrability and shift


@dataclass(frozen=True)
class RMLResonance:
    energy: float
    widths: tuple[float, ...]  # one per channel of the spin group (signed)


@dataclass(frozen=True)
class RMLSpinGroup:
    """All levels of one J and parity, with the channels they decay into."""

    j: float
    parity: int  # +1 or -1
    channels: tuple[Channel, ...]
    resonances: tuple[RMLResonance, ...]


@dataclass(frozen=True)
class RMLRange:
    """A resolved range in the R-Matrix Limited format (LRF=7).

    Unlike LRF=1-3, every channel is listed explicitly: its particle pair,
    l, channel spin s and radii. Only J and parity are conserved, so a spin
    group may hold several neutron channels (both channel spins of a
    target with spin, or l = 0 and 2 for the same J)."""

    el: float
    eh: float
    ifg: int  # 0: widths in eV; 1: reduced-width amplitudes in eV^(1/2)
    krm: int  # R-matrix approximation: 3 = Reich-Moore
    krl: int  # 1: relativistic kinematics
    pairs: tuple[ParticlePair, ...]
    groups: tuple[RMLSpinGroup, ...]
    lrf: int = 7
    naps: int = 0  # unused: LRF=7 gives both radii per channel

    @property
    def formalism(self) -> str:
        return FORMALISMS[7]

    @property
    def ap(self) -> float:
        """Scattering radius of the s-wave neutron channel (APE), the LRF=7
        counterpart of AP (sigma_p = 4 pi AP^2 at low energy)."""
        n = self.neutron_pair
        for sg in self.groups:
            for ch in sg.channels:
                if ch.pair == n and ch.l == 0:
                    return ch.ape
        return 0.0

    @property
    def neutron_pair(self) -> int:
        """The entrance pair: a neutron (MA = 1) with Q = 0 or MT = 2."""
        for i, p in enumerate(self.pairs):
            if p.ma == 1.0 and (p.q == 0.0 or p.mt == 2):
                return i
        raise ValueError("LRF=7 range without a neutron + target pair")


@dataclass(frozen=True)
class UnresolvedJ:
    """Average parameters for one (l, J) sequence, tabulated in energy (LRF=2)."""

    j: float
    interpolation: int
    dof: tuple[
        float, float, float, float
    ]  # degrees of freedom: competitive, neutron, gamma, fission
    energies: tuple[float, ...]
    spacing: tuple[float, ...]  # average level spacing D
    gx: tuple[float, ...]  # average competitive width
    gn: tuple[float, ...]  # average reduced neutron width
    gg: tuple[float, ...]
    gf: tuple[float, ...]


@dataclass(frozen=True)
class UnresolvedRange:
    """Resonances too dense to resolve: only statistical averages are given."""

    el: float
    eh: float
    lrf: int
    naps: int
    spi: float
    ap: float
    lssf: int  # 1: MF3 already holds the infinitely dilute cross section here
    sequences: tuple[
        tuple[int, float, tuple[UnresolvedJ, ...]], ...
    ]  # (l, awri, per-J)


@dataclass(frozen=True)
class EmptyRange:
    """LRU = 0: no resonances, only a scattering radius."""

    el: float
    eh: float
    spi: float
    ap: float


EnergyRange = ResolvedRange | RMLRange | UnresolvedRange | EmptyRange


@dataclass(frozen=True)
class Isotope:
    zai: int
    abundance: float
    lfw: int  # 1 if unresolved fission widths are given
    ranges: tuple[EnergyRange, ...] = field(default_factory=tuple)

    @property
    def resolved(self) -> ResolvedRange | RMLRange | None:
        found = [r for r in self.ranges if isinstance(r, (ResolvedRange, RMLRange))]
        if len(found) > 1:
            raise ValueError("more than one resolved range")
        return found[0] if found else None

    @property
    def unresolved(self) -> UnresolvedRange | None:
        found = [r for r in self.ranges if isinstance(r, UnresolvedRange)]
        return found[0] if found else None


def read_mf2(cursor: Cursor) -> list[Isotope]:
    """Read an MF2/MT151 section into isotopes and their energy ranges."""
    head = cursor.head()
    isotopes = []
    for _ in range(head.n1):
        iso = cursor.cont()
        ranges = tuple(_read_range(cursor, iso.l2) for _ in range(iso.n1))
        isotopes.append(
            Isotope(zai=round(iso.c1), abundance=iso.c2, lfw=iso.l2, ranges=ranges)
        )
    return isotopes


def _read_range(c: Cursor, lfw: int) -> EnergyRange:
    rng = c.cont()
    el, eh, lru, lrf, nro, naps = rng
    ap_table = None
    if nro and lru != 0:
        _, tab = c.tab1()
        ap_table = Tabulated1D(tab.x, tab.y, tab.breaks, tab.laws)

    if lru == 0:
        s = c.cont()
        return EmptyRange(el, eh, spi=s.c1, ap=s.c2)
    if lru == 1 and lrf in (1, 2, 3):
        return _read_bw_or_rm(c, el, eh, lrf, naps, ap_table)
    if lru == 1 and lrf == 7:
        if ap_table is not None:
            raise NotImplementedError("LRF=7 with an energy-dependent radius")
        return _read_rml(c, el, eh)
    if lru == 2 and lrf == 2 and lfw in (0, 1):
        return _read_unresolved_lrf2(c, el, eh, naps)
    raise NotImplementedError(
        f"resonance range {el:g}-{eh:g} eV with LRU={lru}, LRF={lrf} is not supported yet"
    )


def _read_bw_or_rm(c: Cursor, el, eh, lrf, naps, ap_table) -> ResolvedRange:
    s = c.cont()  # SPI, AP, LAD (RM), 0, NLS, NLSC (RM)
    groups = []
    for _ in range(s.n1):
        lst, v = c.list_record()  # AWRI, APL (RM) or QX (BW), L, LRX, 6*NRS, NRS
        rows = [v[i : i + 6] for i in range(0, len(v), 6)]
        if lrf == 3:  # ER, AJ, GN, GG, GFA, GFB
            res = tuple(
                Resonance(e, j, gn, gg, gfa, gfb) for e, j, gn, gg, gfa, gfb in rows
            )
            apl = lst.c2
        else:  # ER, AJ, GT, GN, GG, GF   (GT = total width, redundant)
            res = tuple(
                Resonance(e, j, gn, gg, gf, 0.0) for e, j, _gt, gn, gg, gf in rows
            )
            apl = 0.0
        lrx = 0 if lrf == 3 else lst.l2
        groups.append(LGroup(l=lst.l1, awri=lst.c1, apl=apl, resonances=res, lrx=lrx))
    rm = lrf == 3
    return ResolvedRange(
        el=el,
        eh=eh,
        lrf=lrf,
        naps=naps,
        spi=s.c1,
        ap=s.c2,
        lad=s.l1 if rm else 0,
        nlsc=s.n2 if rm else 0,
        groups=tuple(groups),
        ap_table=ap_table,
    )


def _read_rml(c: Cursor, el: float, eh: float) -> RMLRange:
    """R-Matrix Limited (ENDF-102 section 2.2.1.6)::

    CONT 0, 0, IFG, KRM, NJS, KRL
    LIST 0, 0, NPP, 0, 12*NPP, 2*NPP       particle pairs, 12 numbers each
    per spin group:
      LIST AJ, PJ, KBK, KPS, 6*NCH, NCH    channels: PPI, L, SCH, BND, APE, APT
      LIST 0, 0, 0, NRS, 6*NX, NX          levels: ER, GAM_1 .. GAM_NCH,
                                           each level padded to whole lines
    """
    head = c.cont()
    ifg, krm, njs, krl = head.l1, head.l2, head.n1, head.n2
    _, v = c.list_record()
    # fmt: off
    pairs = tuple(
        ParticlePair(
            ma=v[i], mb=v[i + 1], za=v[i + 2], zb=v[i + 3], ia=v[i + 4],
            ib=v[i + 5], q=v[i + 6], pnt=round(v[i + 7]), shf=round(v[i + 8]),
            mt=round(v[i + 9]),
        )
        for i in range(0, len(v), 12)
    )
    # fmt: on
    groups = []
    for _ in range(njs):
        lst, v = c.list_record()
        aj, pj, kbk, kps, nch = lst.c1, lst.c2, lst.l1, lst.l2, lst.n2
        if kbk or kps:
            raise NotImplementedError("LRF=7 background R-matrix or phase shifts")
        # fmt: off
        channels = tuple(
            Channel(
                pair=round(v[i]) - 1, l=round(v[i + 1]), s=v[i + 2],
                bnd=v[i + 3], ape=v[i + 4], apt=v[i + 5],
            )
            for i in range(0, 6 * nch, 6)
        )
        # fmt: on
        lst, v = c.list_record()
        nrs = lst.l2
        stride = 6 * (nch // 6 + 1)  # every level starts on a new line
        levels = tuple(
            RMLResonance(v[i], tuple(v[i + 1 : i + 1 + nch]))
            for i in range(0, stride * nrs, stride)
        )
        if aj != 0.0:
            parity = 1 if aj > 0 else -1
        else:
            parity = -1 if pj < 0 else 1
        groups.append(RMLSpinGroup(abs(aj), parity, channels, levels))
    return RMLRange(el, eh, ifg, krm, krl, pairs, tuple(groups))


def _read_unresolved_lrf2(c: Cursor, el, eh, naps) -> UnresolvedRange:
    s = c.cont()  # SPI, AP, LSSF, 0, NLS, 0
    sequences = []
    for _ in range(s.n1):
        lc = c.cont()  # AWRI, 0, L, 0, NJS, 0
        js = []
        for _ in range(lc.n1):
            lst, v = c.list_record()  # AJ, 0, INT, 0, 6*NE+6, NE
            # first six numbers: 0, 0, AMUX, AMUN, AMUG, AMUF (degrees of freedom)
            dof = (v[2], v[3], v[4], v[5])
            rows = [v[i : i + 6] for i in range(6, len(v), 6)]  # ES, D, GX, GN, GG, GF
            js.append(
                UnresolvedJ(
                    j=lst.c1,
                    interpolation=lst.l1,
                    dof=dof,
                    energies=tuple(r[0] for r in rows),
                    spacing=tuple(r[1] for r in rows),
                    gx=tuple(r[2] for r in rows),
                    gn=tuple(r[3] for r in rows),
                    gg=tuple(r[4] for r in rows),
                    gf=tuple(r[5] for r in rows),
                )
            )
        sequences.append((lc.l1, lc.c1, tuple(js)))
    return UnresolvedRange(
        el=el,
        eh=eh,
        lrf=2,
        naps=naps,
        spi=s.c1,
        ap=s.c2,
        lssf=s.l1,
        sequences=tuple(sequences),
    )
