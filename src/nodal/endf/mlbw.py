"""Multi-Level Breit-Wigner (MLBW, LRF=2) cross sections.

MLBW keeps SLBW's single-level capture and fission, but lets levels with the
same spin J interfere in elastic scattering (ENDF-102, Appendix D.1.2). For
each spin group (l, J) the elastic collision-matrix element is

    U_nn = exp(-2 i phi_l) [1 + i sum_r Gn_r / (E'_r - E - i G_r / 2)]

and sigma_n = (pi/k^2) g_J |1 - U_nn|^2, which can never be negative.
(For a single level this reproduces the SLBW elastic formula exactly.)
Spin weight not covered by resonant J values scatters off the hard sphere,
as in Reich-Moore (`potential_terms`).

MLBW is adequate when fission is absent or weak and resonances do not
overlap strongly: here Zr isotopes and U-234.
"""

from __future__ import annotations

import cmath
import math

from .channel import penetrability_shift_phase, wave_number
from .reich_moore import potential_terms, spin_groups
from .resonances import LGroup, ResolvedRange
from .slbw import ResonanceXS, single_resonance


def mlbw(rr: ResolvedRange, energy: float) -> ResonanceXS:
    """MLBW elastic, capture and fission cross sections (barns) at one energy."""
    if any(g.lrx for g in rr.groups):
        raise NotImplementedError("MLBW with a competitive width (LRX=1)")
    elastic = capture = fission = 0.0
    for sg in spin_groups(rr):
        lgroup = LGroup(l=sg.l, awri=sg.awri, apl=sg.ap, resonances=sg.resonances)
        k = wave_number(energy, sg.awri)
        _, _, phi = penetrability_shift_phase(sg.l, k * sg.a, k * sg.ap)
        amplitude = 1.0 + 0j
        g = 0.0
        for res in sg.resonances:
            c = single_resonance(rr, lgroup, res, energy)
            capture += c.capture  # no interference in capture or fission
            fission += c.fission
            amplitude += 1j * c.gn / complex(c.shifted_energy - energy, -0.5 * c.width)
            g = c.g
        u_nn = cmath.exp(-2j * phi) * amplitude
        elastic += math.pi / (k * k) * g * abs(1.0 - u_nn) ** 2
    elastic += sum(term(energy) for term in potential_terms(rr))
    return ResonanceXS(elastic=elastic, capture=capture, fission=fission)
