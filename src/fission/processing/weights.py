"""Weighting spectra for group averaging.

A group constant is a flux-weighted average, sigma_g = Int sigma phi / Int phi
over the group. The true flux phi(E) is what the transport calculation will
compute, so the library is built with a *typical* spectrum w(E) instead;
within each group only its shape matters.

`ThermalSlowingDownFission` is NJOY's analytic weight (GROUPR iwt = 4):

    w(E) = a_b E exp(-E/kT)        E <= E_b   thermal Maxwellian
         = 1/E                     E_b < E <= E_c   slowing down
         = a_c sqrt(E) exp(-E/T_f) E > E_c    fission spectrum

with a_b, a_c chosen to make w continuous. NJOY leaves the breaks E_b, E_c to
the user. Here they are where the logarithmic slopes also match, so that w is
smooth: d ln w / d ln E = 1 - E/kT = -1 gives E_b = 2 kT, and
1/2 - E/T_f = -1 gives E_c = 1.5 T_f (NJOY's iwt = 6 uses the same 2 kT).
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from .doppler import BOLTZMANN_EV_PER_K

# Fission temperature from the U-235 thermal prompt spectrum of
# ENDF/B-VIII.1: <E'> = 2.0002 MeV (chapter 9) = 1.5 T_f for a Maxwellian.
FISSION_TEMPERATURE_U235 = 2.0002e6 / 1.5


@dataclass(frozen=True)
class Constant:
    breakpoints: tuple[float, ...] = ()

    def __call__(self, energy: float) -> float:
        return 1.0

    def describe(self) -> str:
        return "constant"


@dataclass(frozen=True)
class InverseE:
    breakpoints: tuple[float, ...] = ()

    def __call__(self, energy: float) -> float:
        return 1.0 / energy

    def describe(self) -> str:
        return "1/E"


@dataclass(frozen=True)
class ThermalSlowingDownFission:
    """Maxwellian -- 1/E -- fission spectrum (NJOY iwt = 4), slope-matched breaks."""

    kt: float  # thermal temperature (eV)
    tf: float = FISSION_TEMPERATURE_U235  # fission temperature (eV)

    @classmethod
    def at(cls, temperature: float, tf: float = FISSION_TEMPERATURE_U235):
        """Thermal part at a moderator temperature in kelvin."""
        return cls(BOLTZMANN_EV_PER_K * temperature, tf)

    @property
    def eb(self) -> float:
        return 2.0 * self.kt

    @property
    def ec(self) -> float:
        return 1.5 * self.tf

    @property
    def ab(self) -> float:
        return 1.0 / (math.exp(-self.eb / self.kt) * self.eb**2)

    @property
    def ac(self) -> float:
        return 1.0 / (math.exp(-self.ec / self.tf) * self.ec**1.5)

    @property
    def breakpoints(self) -> tuple[float, ...]:
        return (self.eb, self.ec)

    def __call__(self, energy: float) -> float:
        if energy <= self.eb:
            return self.ab * energy * math.exp(-energy / self.kt)
        if energy <= self.ec:
            return 1.0 / energy
        return self.ac * math.sqrt(energy) * math.exp(-energy / self.tf)

    def describe(self) -> str:
        return f"thermal({self.kt:.5g} eV)-1/E-fission({self.tf:.5g} eV)"


type Weight = Constant | InverseE | ThermalSlowingDownFission
