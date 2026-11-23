"""R-Matrix Limited (LRF=7) cross sections, Reich-Moore approximation (KRM=3).

ENDF-102 Appendix D.1.7. LRF=7 is the general container of R-matrix
evaluations. Each spin group (one J and parity) lists its channels
explicitly: particle pair, orbital momentum l, channel spin s and radii.
With KRM=3 the R-matrix is the Reich-Moore one, the same physics as LRF=3
(`reich_moore.py`), only written without the shortcuts of that format:

    R_cc'(E) = sum_lambda gamma_lc gamma_lc' / (E_l - E - i G_lg / 2)      (D.44)

over the levels lambda of the group. The capture channels are *eliminated*:
they do not appear as rows of R, their width G_lg only makes the level
energies complex. Widths are given at the level energy (IFG = 0), so the
reduced amplitude is gamma_lc = sign(GAM) sqrt(|GAM| / (2 P_c(|E_l|))) and
the energy-dependent width amplitude sqrt(2 P_c(E)) gamma_lc is

    a_lc(E) = sign(GAM) sqrt(|GAM| P_c(E) / P_c(|E_l|)).

Without shift factors (SHF = 0), L = i P and the matrix of D.42 is the one
of reich_moore.py:

    K = i P^1/2 R P^1/2 = (i/2) sum_l a_l a_l^T / (E_l - E - i G_lg/2),
    U_cc' = exp(-i phi_c) (2 [(1 - K)^-1]_cc' - delta_cc') exp(-i phi_c').

Then, for an incident neutron channel c, summed with g_J over the groups:

    total    (2 pi / k^2) g_J (1 - Re U_cc)                     (optical theorem)
    elastic  (pi / k^2)   g_J sum_c' |delta_cc' - U_cc'|^2       (c' neutron channels)
    capture  (pi / k^2)   g_J sum_l G_lg |A_l|^2,  A_l = sum_c' a_lc' Y_c'c / (E_l - E - i G_lg/2)

with Y = (1 - K)^-1, as in reich_moore.py: capture is the flux the levels
lose through their eliminated gamma widths, computed directly.

Potential scattering needs no extra terms: every (l, s, J) the evaluator
wants is listed as a channel, with or without levels, and a channel
without levels gives U = exp(-2 i phi), |1 - U|^2 = 4 sin^2 phi.

What is not done (errors, not silent): particle pairs other than the
gamma and the neutron + target ground state (inelastic or fission
channels, e.g. Fe-57's 14.4 keV level), shift factors, IFG = 1, other KRM.
"""

from __future__ import annotations

import cmath
import math
from dataclasses import dataclass

from .channel import K_CONSTANT, penetrability_shift_phase
from .reich_moore import _solve
from .resonances import RMLRange, RMLSpinGroup
from .slbw import ResonanceXS


@dataclass(frozen=True)
class PreparedGroup:
    """A spin group reduced to what the energy loop needs."""

    g: float  # spin statistical factor g_J
    mass: float  # target mass of the neutron pair (neutron masses)
    channels: tuple[tuple[int, float, float], ...]  # (l, APT, APE) per neutron channel
    levels: tuple[tuple[float, float, tuple[float, ...], tuple[float, ...]], ...]
    # (E_l, G_lg, |GAM| per neutron channel, P_c(|E_l|) per neutron channel)
    signs: tuple[tuple[float, ...], ...]  # sign of GAM per level and channel


def _k(energy: float, mass: float) -> float:
    """Centre-of-mass wave number (D.39) for a neutron on a target of `mass`
    neutron masses, in 1/(10^-12 cm)."""
    return K_CONSTANT * mass / (mass + 1.0) * math.sqrt(abs(energy))


def prepare(rr: RMLRange) -> list[PreparedGroup]:
    """Check what this implementation supports and precompute per level."""
    if rr.krm != 3:
        raise NotImplementedError(f"LRF=7 with KRM={rr.krm} (only Reich-Moore, 3)")
    if rr.ifg != 0:
        raise NotImplementedError("LRF=7 with reduced-width amplitudes (IFG=1)")
    if rr.krl != 0:
        raise NotImplementedError("LRF=7 with relativistic kinematics")
    n = rr.neutron_pair
    pair = rr.pairs[n]
    for p in rr.pairs:
        if not p.is_gamma and p is not pair:
            raise NotImplementedError(f"LRF=7 particle pair with MT={p.mt}")
    if pair.shf == 1:
        raise NotImplementedError("LRF=7 with shift factors")
    out = []
    for sg in rr.groups:
        neutron = [i for i, c in enumerate(sg.channels) if c.pair == n]
        gammas = [i for i, c in enumerate(sg.channels) if rr.pairs[c.pair].is_gamma]
        if len(neutron) + len(gammas) != len(sg.channels):
            raise NotImplementedError("LRF=7 channel of an unsupported pair")
        chans = tuple(
            (sg.channels[i].l, sg.channels[i].apt, sg.channels[i].ape) for i in neutron
        )
        levels, signs = [], []
        for r in sg.resonances:
            kr = _k(r.energy, pair.mb)
            pr = tuple(penetrability_shift_phase(l, kr * apt, kr * ape)[0]
                       for l, apt, ape in chans)  # fmt: skip
            gg = sum(abs(r.widths[i]) for i in gammas)
            levels.append((r.energy, gg, tuple(abs(r.widths[i]) for i in neutron), pr))
            signs.append(tuple(math.copysign(1.0, r.widths[i]) for i in neutron))
        g = (2.0 * sg.j + 1.0) / (
            (2.0 * abs(pair.ia) + 1.0) * (2.0 * abs(pair.ib) + 1.0)
        )
        out.append(PreparedGroup(g, pair.mb, chans, tuple(levels), tuple(signs)))
    return out


@dataclass(frozen=True)
class GroupSolution:
    k: float
    u: list[list[complex]]  # collision matrix between the neutron channels
    gamma_flux: list[float]  # per incident channel: sum_l G_lg |A_l|^2


def solve_group(pg: PreparedGroup, energy: float) -> GroupSolution:
    """Build K for one spin group, invert (1 - K) column by column."""
    k = _k(energy, pg.mass)
    nch = len(pg.channels)
    p, phi = [], []
    for l, apt, ape in pg.channels:
        pc, _, ph = penetrability_shift_phase(l, k * apt, k * ape)
        p.append(pc)
        phi.append(ph)
    kmat = [[0j] * nch for _ in range(nch)]
    levels = []  # (amplitudes a_lc(E), complex denominator, G_lg)
    for (er, gg, widths, pr), signs in zip(pg.levels, pg.signs, strict=True):
        amp = [s * math.sqrt(w * pc / pe) if pe > 0 else 0.0
               for s, w, pc, pe in zip(signs, widths, p, pr, strict=True)]  # fmt: skip
        denom = complex(er - energy, -0.5 * gg)
        levels.append((amp, denom, gg))
        scale = 0.5j / denom
        for c in range(nch):
            for c2 in range(nch):
                kmat[c][c2] += scale * amp[c] * amp[c2]
    one_minus_k = [[(1.0 if c == c2 else 0.0) - kmat[c][c2] for c2 in range(nch)]
                   for c in range(nch)]  # fmt: skip
    # Y = (1 - K)^-1 is symmetric: column c = solution of (1 - K) y = e_c
    cols = [_solve(one_minus_k, [1.0 + 0j if i == c else 0j for i in range(nch)])
            for c in range(nch)]  # fmt: skip
    omega = [cmath.exp(-1j * ph) for ph in phi]
    u = [[omega[c] * (2.0 * cols[c][c2] - (1.0 if c == c2 else 0.0)) * omega[c2]
          for c2 in range(nch)] for c in range(nch)]  # fmt: skip
    flux = []
    for c in range(nch):
        f = 0.0
        for amp, denom, gg in levels:
            a = sum(ac * y for ac, y in zip(amp, cols[c], strict=True)) / denom
            f += gg * abs(a) ** 2
        flux.append(f)
    return GroupSolution(k, u, flux)


class RMatrixLimited:
    """Callable: energy -> ResonanceXS (elastic, capture; no fission)."""

    def __init__(self, rr: RMLRange):
        self.groups = prepare(rr)

    def __call__(self, energy: float) -> ResonanceXS:
        elastic = capture = 0.0
        for pg in self.groups:
            sol = solve_group(pg, energy)
            pk = math.pi / (sol.k * sol.k) * pg.g
            nch = len(sol.u)
            for c in range(nch):
                elastic += pk * sum(abs((1.0 if c == c2 else 0.0) - sol.u[c][c2]) ** 2
                                    for c2 in range(nch))  # fmt: skip
                capture += pk * sol.gamma_flux[c]
        return ResonanceXS(elastic=elastic, capture=capture, fission=0.0)

    def total(self, energy: float) -> float:
        """Optical theorem, (2 pi / k^2) g_J sum_c (1 - Re U_cc): an
        independent route to elastic + capture."""
        total = 0.0
        for pg in self.groups:
            sol = solve_group(pg, energy)
            total += (
                2.0
                * math.pi
                / (sol.k * sol.k)
                * pg.g
                * sum(1.0 - sol.u[c][c].real for c in range(len(sol.u)))
            )
        return total


def r_matrix_limited(rr: RMLRange, energy: float) -> ResonanceXS:
    """Elastic and capture (barns) at one energy (rebuilds the groups; for
    many energies keep an RMatrixLimited)."""
    return RMatrixLimited(rr)(energy)


def levels_of(sg: RMLSpinGroup) -> list[tuple[float, float]]:
    """(energy, total width) of every level, for seeding a reconstruction."""
    return [(r.energy, sum(abs(w) for w in r.widths)) for r in sg.resonances]
