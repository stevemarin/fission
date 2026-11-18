"""Channel quantities shared by every resonance formalism.

A neutron of energy E meeting a nucleus is described in the centre-of-mass
frame by its wave number k (inverse reduced wavelength). With orbital angular
momentum l, it must tunnel through a centrifugal barrier to reach the nuclear
surface at the channel radius a. Outside the nucleus the neutron moves freely,
so everything is expressed through the dimensionless rho = k a and the
"hard-sphere" functions of rho:

    P_l(rho)    penetrability: how easily the neutron reaches the surface.
                For small rho, P_l ~ rho^(2l+1): s-waves (l=0) get in freely,
                p-waves (l=1) are suppressed by rho^2 at low energy.
    S_l(rho)    shift factor: moves the resonance energy slightly with E.
    phi_l(rho)  hard-sphere phase shift: scattering off an impenetrable ball,
                the source of the smooth "potential scattering" 4 pi a^2.

Units follow ENDF: energies in eV, radii in 10^-12 cm, so pi/k^2 comes out
directly in barns (1 b = 10^-24 cm^2 = (10^-12 cm)^2).
"""

from __future__ import annotations

import math

# CODATA 2018
NEUTRON_MASS_EV = 939.56542052e6  # m_n c^2 (eV)
HBAR_C_EV_CM = 197.3269804e6 * 1e-13  # hbar c (eV cm): 197.327 MeV fm

# k = sqrt(2 m_n E) / hbar, expressed in 1/(10^-12 cm):  k = K_CONSTANT * sqrt(E [eV])
K_CONSTANT = math.sqrt(2.0 * NEUTRON_MASS_EV) / HBAR_C_EV_CM * 1e-12  # ~ 2.19680e-3


def wave_number(energy: float, awri: float) -> float:
    """Centre-of-mass neutron wave number k in 1/(10^-12 cm).

    awri / (awri + 1) converts the lab-frame energy to the centre of mass
    (the target recoils, so the neutron has slightly less relative momentum).
    """
    return K_CONSTANT * awri / (awri + 1.0) * math.sqrt(abs(energy))


def channel_radius(awri: float, ap: float, naps: int) -> float:
    """Channel radius a (10^-12 cm) used in penetrability and shift factors.

    NAPS = 0: a from the nuclear-size formula a = 0.123 AWRI^(1/3) + 0.08
    NAPS = 1: a = AP, the scattering radius
    The phase shift always uses the scattering radius AP (or APL).
    """
    if naps == 0:
        return 0.123 * awri ** (1.0 / 3.0) + 0.08
    if naps == 1:
        return ap
    raise ValueError(f"NAPS={naps} is not defined for resolved resonances")


def penetrability_shift_phase(
    l: int, rho: float, rho_hat: float
) -> tuple[float, float, float]:
    """(P_l(rho), S_l(rho), phi_l(rho_hat)) for a hard sphere, l = 0 ... 3.

    rho = k a (channel radius), rho_hat = k AP (scattering radius).
    """
    r2 = rho * rho
    h2 = rho_hat * rho_hat
    if l == 0:
        return rho, 0.0, rho_hat
    if l == 1:
        return r2 * rho / (1 + r2), -1.0 / (1 + r2), rho_hat - math.atan(rho_hat)
    if l == 2:
        d = 9 + 3 * r2 + r2 * r2
        return (
            r2 * r2 * rho / d,
            -(18 + 3 * r2) / d,
            rho_hat - math.atan2(3 * rho_hat, 3 - h2),
        )
    if l == 3:
        d = 225 + 45 * r2 + 6 * r2 * r2 + r2 * r2 * r2
        return (
            r2 * r2 * r2 * rho / d,
            -(675 + 90 * r2 + 6 * r2 * r2) / d,
            rho_hat - math.atan2(rho_hat * (15 - h2), 15 - 6 * h2),
        )
    raise NotImplementedError(f"l = {l} not implemented")


def spin_factor(j: float, spi: float) -> float:
    """Statistical weight g_J = (2J+1) / ((2i+1)(2I+1)), neutron spin i = 1/2.

    The fraction of neutron-target spin orientations that can form a
    compound state of spin J.
    """
    return (2.0 * abs(j) + 1.0) / (2.0 * (2.0 * spi + 1.0))
