"""MF7: thermal scattering laws S(alpha, beta) [ENDF-102 §7].

Below a few eV the neutron's energy is comparable with the thermal motion
of the moderator atoms, and for hydrogen in water with the energies of
molecular rotation and vibration. Scattering then depends on how the atom is
bound, which the scattering law S(alpha, beta, T) describes. Per principal
atom, integrated over the azimuth,

    d2 sigma / dE' dmu = sigma_b / (2 kT) sqrt(E'/E) exp(-beta/2) S(alpha, beta)

with the energy transfer beta = (E' - E)/kT and the momentum transfer
alpha = (E + E' - 2 mu sqrt(E E')) / (A0 kT). sigma_b = sigma_free
((A0+1)/A0)^2 is the bound-atom cross section.

The file tabulates S for the principal atom (H in H2O) at a list of
temperatures; other atoms of the molecule are described by constants
(oxygen in water: a free gas). With LAT = 1 the tabulated alpha and beta are
computed with kT0 = 0.0253 eV instead of the actual kT (only the lookup is
rescaled; the prefactor and exp(-beta/2) use the actual temperature, as in
NJOY's THERMR).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..evaluation import Evaluation
from ..interpolation import Tabulated1D

KT0 = 0.0253  # eV, the reference temperature of LAT = 1
LN_S_FLOOR = -225.0  # ln S used for S <= 0 ("negligible"), as in NJOY


@dataclass(frozen=True)
class OtherAtom:
    """A non-principal atom of the molecule (B(7..12) etc.)."""

    kind: int  # 0: short-collision-time only, 1: free gas, 2: diffusive
    sigma_free: float  # M_n sigma_fn, all atoms of this type (b)
    awr: float
    count: int


@dataclass(frozen=True)
class ScatteringLaw:
    """S(alpha, beta) of the principal atom at one temperature."""

    za: int
    temperature: float  # K
    lat: int
    lasym: int
    sigma_free: float  # per principal atom (b)
    a0: float  # mass ratio used in alpha
    emax: float  # upper energy limit of the data (eV)
    count: int  # principal atoms per molecule (2 for H in H2O)
    others: tuple[OtherAtom, ...]
    alpha: np.ndarray  # increasing
    beta: np.ndarray  # increasing (>= 0 when LASYM = 0)
    ln_s: np.ndarray  # ln S, shape (len(beta), len(alpha))
    t_effective: float  # short-collision-time temperature (K)

    @property
    def sigma_bound(self) -> float:
        """Bound-atom cross section per principal atom (b)."""
        return self.sigma_free * ((self.a0 + 1.0) / self.a0) ** 2


def available_temperatures(ev: Evaluation) -> list[float]:
    cur = ev.cursor(7, 4)
    cur.cont()
    cur.list_record()
    cur.tab2()
    c, _ = cur.tab1()
    temps = [c.c1]
    for _ in range(c.l1):
        cl, _ = cur.list_record()
        temps.append(cl.c1)
    return temps


def scattering_law(ev: Evaluation, temperature: float) -> ScatteringLaw:
    """S(alpha, beta) at one of the file's temperatures (MF7 MT4).

    ENDF-102 §7.5.2 advises against interpolating S in temperature, so only
    tabulated temperatures (within 0.5 K) are accepted.
    """
    cur = ev.cursor(7, 4)
    head = cur.cont()
    za, lat, lasym = int(head.c1), head.l2, head.n1
    lb, b = cur.list_record()
    lln, ns = lb.l1, lb.n2
    if b[0] == 0.0:
        raise NotImplementedError("analytic-only scattering laws (B(1) = 0)")
    t2, _, beta_laws = cur.tab2()
    betas, rows, alpha = [], [], None
    index = None
    alpha_laws: list[int] = []
    temps: list[float] = []
    for _ in range(t2.n2):
        c, t = cur.tab1()
        temps = [c.c1]
        values = [t.y]
        for _ in range(c.l1):
            cl, v = cur.list_record()
            temps.append(cl.c1)
            values.append(v)
        if index is None:
            matches = [i for i, tt in enumerate(temps) if abs(tt - temperature) <= 0.5]
            if not matches:
                raise KeyError(
                    f"no S(alpha, beta) at {temperature} K; available: {temps}"
                )
            index = matches[0]
        if alpha is None:
            alpha, alpha_laws = t.x, t.laws
        betas.append(c.c2)
        rows.append(values[index])
    assert alpha is not None and index is not None
    if set(beta_laws) - {2, 4} or set(alpha_laws) - {2, 4}:
        raise NotImplementedError(
            "S(alpha, beta) interpolation other than lin-lin/log-lin"
        )
    s = np.array(rows, dtype=float)
    ln_s = s if lln else np.where(s > 0, np.log(np.where(s > 0, s, 1.0)), LN_S_FLOOR)
    # effective temperatures for the short-collision-time approximation
    c, teff = cur.tab1()
    t_effective = Tabulated1D(teff.x, teff.y, teff.breaks, teff.laws)(temps[index])
    others = tuple(
        OtherAtom(int(b[6 * n]), b[6 * n + 1], b[6 * n + 2], int(b[6 * n + 5]))
        for n in range(1, ns + 1)
    )
    count = int(b[5])
    return ScatteringLaw(
        za=za, temperature=temps[index], lat=lat, lasym=lasym,
        sigma_free=b[0] / count, a0=b[2], emax=b[3], count=count, others=others,
        alpha=np.asarray(alpha, dtype=float), beta=np.asarray(betas, dtype=float),
        ln_s=ln_s, t_effective=t_effective,
    )  # fmt: skip
