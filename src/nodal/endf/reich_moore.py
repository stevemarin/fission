"""Reich-Moore (LRF=3) cross sections from resolved resonance parameters.

R-matrix theory, the short version (ENDF-102, Appendix D.1.3):

Neutrons, fission fragments and gammas leave the compound nucleus through
*channels*. Resonances with the same orbital momentum l and spin J (a "spin
group") share the same channels, so their amplitudes add before squaring: they
interfere. For each spin group build the channel matrix

    K_cc'(E) = (i/2) sum_r  sqrt(G_rc) sqrt(G_rc') / (E_r - E - i G_rg / 2)

over the open channels c, c' (neutron, fission a, fission b). sqrt(G_rc) is the
signed amplitude of resonance r in channel c; the neutron width carries its
penetrability energy dependence G_rn(E) = G_rn(|E_r|) P_l(E)/P_l(|E_r|).

The Reich-Moore trick is in the denominator: the thousands of gamma channels
are not carried explicitly, their total width G_rg just makes each level energy
complex. That keeps the matrix at most 3x3 while staying very accurate.

From the collision matrix element U_nn = exp(-2 i phi_l) (2 [(1 - K)^-1]_nn - 1):

    total    sigma_t = (2 pi/k^2) g_J (1 - Re U_nn)
    elastic  sigma_n = (pi/k^2)   g_J |1 - U_nn|^2
    fission  sigma_f = (4 pi/k^2) g_J sum_c in fission |[(1 - K)^-1]_nc|^2
    capture  sigma_g = (pi/k^2)   g_J sum_r G_rg |A_r|^2   (flux lost through gammas)

with x = (1 - K)^-1 e_n and level amplitudes A_r = sum_c sqrt(G_rc) x_c / (E_r - E - i G_rg/2).
Capture could also be obtained as sigma_t - sigma_n - sigma_f, but that
difference cancels catastrophically far from resonances.

Unitarity |U_nn| <= 1 guarantees sigma_n >= 0, which SLBW cannot promise.
Spin values J of an l that have no resonances still scatter off the hard
sphere: they add (4 pi/k^2) g_J sin^2(phi_l), and l values up to NLSC - 1
contribute such potential scattering even without any resonances.
"""

from __future__ import annotations

import cmath
import math
from collections import defaultdict
from dataclasses import dataclass

from .channel import channel_radius, penetrability_shift_phase, spin_factor, wave_number
from .resonances import ResolvedRange, Resonance
from .slbw import ResonanceXS


@dataclass(frozen=True)
class SpinGroup:
    """Resonances sharing orbital momentum l and spin J, with their radii."""

    l: int
    j: float
    awri: float
    a: float  # channel radius (penetrability)
    ap: float  # scattering radius (phase shift)
    resonances: tuple[Resonance, ...]

    @property
    def has_fission(self) -> bool:
        return any(r.gfa != 0 or r.gfb != 0 for r in self.resonances)


def spin_groups(rr: ResolvedRange) -> list[SpinGroup]:
    """Split a Reich-Moore range into (l, J) spin groups."""
    if rr.lrf != 3:
        raise ValueError(f"not a Reich-Moore range (LRF={rr.lrf})")
    groups = []
    for lg in rr.groups:
        if lg.l > 0 and rr.spi > 0 and any(r.j < 0 for r in lg.resonances):
            # a negative J selects one of two channel spins; that needs two
            # neutron channels per J, which our nuclides never require
            raise NotImplementedError("channel-spin selection via negative J")
        ap = lg.apl if lg.apl > 0 else rr.ap
        a = channel_radius(lg.awri, ap, rr.naps)
        by_j: dict[float, list[Resonance]] = defaultdict(list)
        for r in lg.resonances:
            by_j[abs(r.j)].append(r)
        for j in sorted(by_j):
            groups.append(SpinGroup(lg.l, j, lg.awri, a, ap, tuple(by_j[j])))
    return groups


def _solve(matrix: list[list[complex]], rhs: list[complex]) -> list[complex]:
    """Solve a small complex linear system by Gaussian elimination with pivoting."""
    n = len(rhs)
    m = [row[:] + [rhs[i]] for i, row in enumerate(matrix)]
    for col in range(n):
        pivot = max(range(col, n), key=lambda r: abs(m[r][col]))
        m[col], m[pivot] = m[pivot], m[col]
        for r in range(col + 1, n):
            f = m[r][col] / m[col][col]
            for c in range(col, n + 1):
                m[r][c] -= f * m[col][c]
    x = [0j] * n
    for r in reversed(range(n)):
        x[r] = (m[r][n] - sum(m[r][c] * x[c] for c in range(r + 1, n))) / m[r][r]
    return x


@dataclass(frozen=True)
class SpinGroupSolution:
    """One spin group solved at one energy."""

    u_nn: complex  # elastic collision-matrix element
    row: list[complex]  # x = (1 - K)^-1 e_n: [x_n, x_fa, x_fb]
    k: float  # wave number
    phi: float  # hard-sphere phase shift
    gamma_flux: float  # sum_r G_rg |A_r|^2, the flux absorbed through gamma emission


def solve_spin_group(sg: SpinGroup, energy: float) -> SpinGroupSolution:
    """Build K for one spin group, solve (1 - K) x = e_n and collect the results."""
    k = wave_number(energy, sg.awri)
    p, _, phi = penetrability_shift_phase(sg.l, k * sg.a, k * sg.ap)
    nch = 3 if sg.has_fission else 1
    kmat = [[0j] * nch for _ in range(nch)]
    levels = []  # (channel amplitudes, complex energy denominator, gamma width)
    for r in sg.resonances:
        kr = wave_number(r.energy, sg.awri)
        pr, _, _ = penetrability_shift_phase(sg.l, kr * sg.a, kr * sg.ap)
        amp = [math.sqrt(r.gn * p / pr)]
        if nch == 3:
            amp += [math.copysign(math.sqrt(abs(f)), f) for f in (r.gfa, r.gfb)]
        denom = complex(r.energy - energy, -0.5 * r.gg)
        levels.append((amp, denom, r.gg))
        scale = 0.5j / denom
        for c in range(nch):
            for c2 in range(nch):
                kmat[c][c2] += scale * amp[c] * amp[c2]
    one_minus_k = [
        [(1.0 if c == c2 else 0.0) - kmat[c][c2] for c2 in range(nch)]
        for c in range(nch)
    ]
    # (1 - K) is symmetric, so row n of its inverse = solution of (1 - K) x = e_n
    row = _solve(one_minus_k, [1.0 + 0j] + [0j] * (nch - 1))
    u_nn = cmath.exp(-2j * phi) * (2.0 * row[0] - 1.0)

    # Capture directly, as the flux each level loses through its gamma width:
    # level amplitude A_r = sum_c sqrt(G_rc) x_c / (E_r - E - i G_rg/2).
    # (Computing it as total - elastic - fission would cancel catastrophically
    # far from resonances, where capture is 1e-6 of elastic.)
    gamma_flux = 0.0
    for amp, denom, gg in levels:
        a_r = sum(a * x for a, x in zip(amp, row, strict=True)) / denom
        gamma_flux += gg * abs(a_r) ** 2
    return SpinGroupSolution(u_nn, row, k, phi, gamma_flux)


@dataclass(frozen=True)
class PotentialTerm:
    """Hard-sphere scattering of (l, J) combinations that have no resonances.

    `weight` is the spin weight not covered by resonant spin groups: summed
    over all J, g_J adds up to 2l+1, so weight = (2l+1) - sum of resonant g_J.
    """

    l: int
    weight: float
    awri: float
    a: float
    ap: float

    def __call__(self, energy: float) -> float:
        k = wave_number(energy, self.awri)
        _, _, phi = penetrability_shift_phase(self.l, k * self.a, k * self.ap)
        return 4.0 * math.pi / (k * k) * self.weight * math.sin(phi) ** 2


def potential_terms(rr: ResolvedRange) -> list[PotentialTerm]:
    """Potential scattering not carried by any spin group, including l values
    listed only so that potential scattering converges (l < NLSC)."""
    g_used: dict[int, float] = defaultdict(float)
    for sg in spin_groups(rr):
        g_used[sg.l] += spin_factor(sg.j, rr.spi)
    awri = rr.groups[0].awri
    radii = {lg.l: (lg.apl if lg.apl > 0 else rr.ap) for lg in rr.groups}
    terms = []
    for l in range(max(rr.nlsc, len(rr.groups))):
        missing = (2 * l + 1) - g_used[l]
        if missing > 1e-12:
            ap = radii.get(l, rr.ap)
            terms.append(
                PotentialTerm(l, missing, awri, channel_radius(awri, ap, rr.naps), ap)
            )
    return terms


def reich_moore(rr: ResolvedRange, energy: float) -> ResonanceXS:
    """Reich-Moore elastic, capture and fission cross sections (barns) at one energy."""
    elastic = capture = fission = 0.0
    for sg in spin_groups(rr):
        sol = solve_spin_group(sg, energy)
        g = spin_factor(sg.j, rr.spi)
        pk = math.pi / (sol.k * sol.k)
        elastic += pk * g * abs(1.0 - sol.u_nn) ** 2
        fission += 4.0 * pk * g * sum(abs(x) ** 2 for x in sol.row[1:])
        capture += pk * g * sol.gamma_flux
    elastic += sum(term(energy) for term in potential_terms(rr))
    return ResonanceXS(elastic=elastic, capture=capture, fission=fission)


def total_from_optical_theorem(rr: ResolvedRange, energy: float) -> float:
    """sigma_t = (2 pi/k^2) sum g_J (1 - Re U_nn) + missing-J potential scattering.

    An independent route to the total; equals elastic + capture + fission
    when the collision matrix is unitary (used as a consistency check).
    """
    total = 0.0
    for sg in spin_groups(rr):
        sol = solve_spin_group(sg, energy)
        g = spin_factor(sg.j, rr.spi)
        total += 2.0 * math.pi / (sol.k * sol.k) * g * (1.0 - sol.u_nn.real)
    return total + sum(term(energy) for term in potential_terms(rr))


def main(argv: list[str] | None = None) -> None:
    """Compare Reich-Moore and SLBW across one resonance.

    uv run python -m nodal.endf.reich_moore U238 6.674
    """
    import argparse

    from ..data import EndfLibrary
    from .evaluation import Evaluation
    from .slbw import single_resonance, slbw

    ap = argparse.ArgumentParser(
        prog="python -m nodal.endf.reich_moore", description=main.__doc__
    )
    ap.add_argument("nuclide")
    ap.add_argument("energy", type=float, help="approximate resonance energy (eV)")
    ap.add_argument(
        "--widths", type=float, default=6.0, help="half-range in total widths"
    )
    args = ap.parse_args(argv)

    rr = (
        Evaluation.from_file(EndfLibrary().neutron(args.nuclide)).resonances[0].resolved
    )
    if rr is None:
        raise SystemExit(f"{args.nuclide} has no resolved resonances")
    group, res = min(
        ((g, r) for g in rr.groups for r in g.resonances),
        key=lambda gr: abs(gr[1].energy - args.energy),
    )
    width = single_resonance(rr, group, res, res.energy).width
    print(
        f"{args.nuclide} {rr.formalism}, resonance at {res.energy:g} eV (l={group.l}, J={res.j})"
    )
    print(
        f"  {'E (eV)':>10} | {'RM capture':>11} {'RM elastic':>11} {'RM fission':>11} |"
        f" {'SLBW capt.':>11} {'SLBW elas.':>11}"
    )
    for i in range(-10, 11):
        e = res.energy + args.widths * width * i / 10
        rm, bw = reich_moore(rr, e), slbw(rr, e)
        print(
            f"  {e:10.5f} | {rm.capture:11.5g} {rm.elastic:11.5g} {rm.fission:11.5g} |"
            f" {bw.capture:11.5g} {bw.elastic:11.5g}"
        )


if __name__ == "__main__":
    main()
