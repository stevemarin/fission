"""Steady-state temperatures of a fuel rod.

Heat generated in the pellet at q' (W per metre of rod) flows out through
the pellet, the gas gap and the clad into the coolant. In cylindrical
geometry, with uniform heat generation:

* pellet: int_{T_s}^{T_c} k(T) dT = q' / (4 pi) (centre T_c, surface T_s);
* gap: T_s - T_ci = q' / (2 pi r_p h_gap);
* clad: T_ci - T_co = q' ln(r_co / r_ci) / (2 pi k_clad).

The properties are the simplified ones of the OECD/NEA PWR MOX/UO2 core
transient benchmark specification (2003), section 7, after Finnemann
(1991) (free: https://www.oecd-nea.org/upload/docs/application/pdf/2020-10/
mox_benchmark_specifications_2003.pdf):

    k_UO2 = 1.05 + 2150 / (T - 73.15) W/(m K)                       Eq. (3)
    k_Zr = 7.51 + 2.09e-2 T - 1.45e-5 T^2 + 7.67e-9 T^3 W/(m K)     Eq. (4)
    h_gap = 1e4 W/(m^2 K)                                          Eq. (8)
    T_Doppler = 0.3 T_c + 0.7 T_s                                  Eq. (2), alpha = 0.7

From the clad to the water, a film coefficient h (``convection``):
T_co - T_water = q' / (2 pi r_co h); without h the clad surface is at the
water temperature.
"""

from __future__ import annotations

import math
from dataclasses import dataclass


def k_fuel(t: float) -> float:
    return 1.05 + 2150.0 / (t - 73.15)


def k_fuel_integral(t: float) -> float:
    """int k dT (an antiderivative), W/m."""
    return 1.05 * t + 2150.0 * math.log(t - 73.15)


def k_clad(t: float) -> float:
    return 7.51 + 2.09e-2 * t - 1.45e-5 * t**2 + 7.67e-9 * t**3


@dataclass
class RodTemperatures:
    centre: float
    surface: float  # pellet surface
    clad_inner: float
    clad_outer: float

    @property
    def doppler(self) -> float:
        return 0.3 * self.centre + 0.7 * self.surface


def rod(q_linear: float, coolant: float, pellet_radius: float, clad_inner: float,
        clad_outer: float, gap_conductance: float = 1e4,
        film: float | None = None) -> RodTemperatures:  # fmt: skip
    """Temperatures (K) for a linear heat rate q' (W/m) and coolant
    temperature (K); radii in metres; film: clad-to-water heat transfer
    coefficient (W/(m^2 K)), None for none."""
    t_co = coolant
    if film is not None:
        t_co += q_linear / (2 * math.pi * clad_outer * film)
    t_ci = t_co
    for _ in range(50):  # clad conductivity at the mean clad temperature
        new = t_co + q_linear * math.log(clad_outer / clad_inner) / (
            2 * math.pi * k_clad(0.5 * (t_ci + t_co))
        )
        if abs(new - t_ci) < 1e-9:
            break
        t_ci = new
    t_s = t_ci + q_linear / (2 * math.pi * pellet_radius * gap_conductance)
    target = k_fuel_integral(t_s) + q_linear / (4 * math.pi)
    t_c = t_s + q_linear / (4 * math.pi * k_fuel(t_s))  # start: constant k
    for _ in range(50):  # Newton on int_{T_s}^{T_c} k dT = q' / (4 pi)
        step = (k_fuel_integral(t_c) - target) / k_fuel(t_c)
        t_c -= step
        if abs(step) < 1e-9:
            break
    return RodTemperatures(t_c, t_s, t_ci, t_co)
