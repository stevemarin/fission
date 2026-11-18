"""Kinematics of neutron emission, and the Kalbach-Mann angular systematics.

Masses are in units of the neutron mass (ENDF's AWR). Non-relativistic.

Two-body reactions  n + A -> b + B  (elastic: b = n, B = A, Q = 0; discrete
inelastic: Q = -E_level). In the centre-of-mass frame the outgoing particle
has a fixed energy

    E'_cm = (M_B / (M_B + m_b)) * (A/(A+1) E + Q),

and its direction is described by the CM cosine mu_cm. Adding back the
motion of the centre of mass gives the laboratory values [ENDF-102 eq. 6.7]:

    E'_lab = E'_cm + m_b E / (A+1)^2 + 2 sqrt(m_b E E'_cm) mu_cm / (A+1)
    mu_lab = sqrt(E'_cm / E'_lab) mu_cm + sqrt(m_b E / E'_lab) / (A+1)

(for a neutron projectile, m_a = 1). For elastic scattering this is the
classic E' = E (A^2 + 2 A mu_cm + 1) / (A+1)^2.
"""

from __future__ import annotations

import math


def cm_energy(energy: float, awr: float, q: float = 0.0, awp: float = 1.0) -> float:
    """Energy of the emitted particle in the CM frame (eV); <= 0 below threshold."""
    m_residual = awr + 1.0 - awp
    return m_residual / (m_residual + awp) * (awr / (awr + 1.0) * energy + q)


def cm_to_lab(
    energy: float, e_cm: float, mu_cm: float, awr: float, awp: float = 1.0
) -> tuple[float, float]:
    """(E'_lab, mu_lab) of a particle emitted with E'_cm and mu_cm."""
    a1 = awr + 1.0
    e_lab = (
        e_cm + awp * energy / a1**2 + 2.0 * math.sqrt(awp * energy * e_cm) * mu_cm / a1
    )
    if e_lab <= 0.0:
        return 0.0, 1.0
    mu_lab = math.sqrt(e_cm / e_lab) * mu_cm + math.sqrt(awp * energy / e_lab) / a1
    return e_lab, max(-1.0, min(1.0, mu_lab))


def two_body(
    energy: float, mu_cm: float, awr: float, q: float = 0.0, awp: float = 1.0
) -> tuple[float, float]:
    """(E'_lab, mu_lab) for a two-body reaction at CM cosine mu_cm."""
    return cm_to_lab(energy, max(cm_energy(energy, awr, q, awp), 0.0), mu_cm, awr, awp)


def elastic_alpha(awr: float) -> float:
    """alpha = ((A-1)/(A+1))^2: the smallest E'/E in elastic scattering."""
    return ((awr - 1.0) / (awr + 1.0)) ** 2


def mean_lethargy_gain(awr: float) -> float:
    """xi = <ln(E/E')> for isotropic CM elastic scattering."""
    alpha = elastic_alpha(awr)
    if alpha == 0.0:
        return 1.0
    return 1.0 + alpha * math.log(alpha) / (1.0 - alpha)


# ------------------------------------------------------------ Kalbach-Mann
# Parameters of Kalbach's systematics as given in ENDF-102 §6.2.3.2 (Table V
# of Kalbach, Phys. Rev. C 37, 2350 (1988)); energies in MeV.
_C1, _C2, _C3 = 0.04, 1.8e-6, 6.7e-7
_ET1, _ET3 = 130.0, 41.0
# Only neutron emission after neutron capture is needed here: projectile and
# emitted particle are both neutrons, with M_a = 1, m_b = 1/2 and no internal
# binding (I = 0). Other particles would need ENDF-102 Appendix H.
_M_NEUTRON_IN, _M_NEUTRON_OUT = 1.0, 0.5


def separation_energy(za_compound: int, za_residual: int) -> float:
    """Kalbach's separation energy (MeV) of a neutron from the compound nucleus.

    A liquid-drop mass formula without pairing and shell terms, as ENDF-102
    prescribes for this purpose. za_* are 1000 Z + A.
    """

    def nza(za: int) -> tuple[int, int, int]:
        z, a = divmod(za, 1000)
        return a - z, z, a

    nc, zc, ac = nza(za_compound)
    nb, zb, ab = nza(za_residual)
    return (
        15.68 * (ac - ab)
        - 28.07 * ((nc - zc) ** 2 / ac - (nb - zb) ** 2 / ab)
        - 18.56 * (ac ** (2 / 3) - ab ** (2 / 3))
        + 33.22 * ((nc - zc) ** 2 / ac ** (4 / 3) - (nb - zb) ** 2 / ab ** (4 / 3))
        - 0.717 * (zc**2 / ac ** (1 / 3) - zb**2 / ab ** (1 / 3))
        + 1.211 * (zc**2 / ac - zb**2 / ab)
    )


def kalbach_slope(
    energy: float, e_out_cm: float, za_target: int, awr_target: float
) -> float:
    """Kalbach's slope a(E, E'_cm) for a neutron emitted after neutron capture.

    energy and e_out_cm in eV. a grows with the emission energy: fast
    pre-compound particles are strongly forward peaked.
    """
    # Kalbach's single-emission picture: every neutron leaves the same compound
    # nucleus (A+1), so the residual is the target and S_a = S_b = S_n.
    s_n = separation_energy(za_target + 1, za_target)
    eps_a = energy * awr_target / (awr_target + 1.0) / 1e6
    eps_b = e_out_cm * (awr_target + 1.0) / awr_target / 1e6
    e_a, e_b = eps_a + s_n, eps_b + s_n
    x1 = min(e_a, _ET1) * e_b / e_a
    x3 = min(e_a, _ET3) * e_b / e_a
    return _C1 * x1 + _C2 * x1**3 + _C3 * _M_NEUTRON_IN * _M_NEUTRON_OUT * x3**4


def kalbach_mann(mu: float, a: float, r: float) -> float:
    """Normalized angular shape a / (2 sinh a) [cosh(a mu) + r sinh(a mu)].

    r is the pre-compound fraction: r = 0 gives a distribution symmetric
    about 90 degrees (compound-nucleus evaporation), r = 1 forward peaked.
    """
    if a < 1e-8:
        return 0.5 * (1.0 + r * a * mu)
    return a / (2.0 * math.sinh(a)) * (math.cosh(a * mu) + r * math.sinh(a * mu))


def kalbach_mann_mean_cosine(a: float, r: float) -> float:
    """<mu> = r (coth a - 1/a): only the pre-compound part is forward peaked."""
    if a < 1e-6:
        return r * a / 3.0
    return r * (1.0 / math.tanh(a) - 1.0 / a)
