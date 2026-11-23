"""Properties of compressed liquid water: IAPWS-IF97, region 1.

IAPWS R7-97(2012), "Revised Release on the IAPWS Industrial Formulation
1997 for the Thermodynamic Properties of Water and Steam", section 5.1
(free: https://iapws.org/technical-guidance/release/IF97-Rev.download).
Region 1 covers liquid water from 273.15 to 623.15 K at pressures above
saturation, which includes a PWR core (15.5 MPa, 565-600 K).

The release gives the dimensionless Gibbs free energy (Eq. 7)

    gamma(pi, tau) = g / (R T) = sum_i n_i (7.1 - pi)^I_i (tau - 1.222)^J_i,

with pi = p / 16.53 MPa, tau = 1386 K / T and the 34 coefficients of
Table 2. Every property follows from its derivatives (Table 3):

    specific volume   v = pi gamma_pi R T / p,
    enthalpy          h = tau gamma_tau R T,
    heat capacity     c_p = -tau^2 gamma_tautau R.

The temperature of a given enthalpy is found by Newton's method on h(T)
(slope c_p).
"""

from __future__ import annotations

import math

import numpy as np

R = 0.461526e3  # J/(kg K), Eq. (1)
P_STAR = 16.53e6  # Pa
T_STAR = 1386.0  # K

# Table 2: I_i, J_i, n_i
_TABLE2 = (
    (0, -2, 0.14632971213167),
    (0, -1, -0.84548187169114),
    (0, 0, -0.37563603672040e1),
    (0, 1, 0.33855169168385e1),
    (0, 2, -0.95791963387872),
    (0, 3, 0.15772038513228),
    (0, 4, -0.16616417199501e-1),
    (0, 5, 0.81214629983568e-3),
    (1, -9, 0.28319080123804e-3),
    (1, -7, -0.60706301565874e-3),
    (1, -1, -0.18990068218419e-1),
    (1, 0, -0.32529748770505e-1),
    (1, 1, -0.21841717175414e-1),
    (1, 3, -0.52838357969930e-4),
    (2, -3, -0.47184321073267e-3),
    (2, 0, -0.30001780793026e-3),
    (2, 1, 0.47661393906987e-4),
    (2, 3, -0.44141845330846e-5),
    (2, 17, -0.72694996297594e-15),
    (3, -4, -0.31679644845054e-4),
    (3, 0, -0.28270797985312e-5),
    (3, 6, -0.85205128120103e-9),
    (4, -5, -0.22425281908000e-5),
    (4, -2, -0.65171222895601e-6),
    (4, 10, -0.14341729937924e-12),
    (5, -8, -0.40516996860117e-6),
    (8, -11, -0.12734301741641e-8),
    (8, -6, -0.17424871230634e-9),
    (21, -29, -0.68762131295531e-18),
    (23, -31, 0.14478307828521e-19),
    (29, -38, 0.26335781662795e-22),
    (30, -39, -0.11947622640071e-22),
    (31, -40, 0.18228094581404e-23),
    (32, -41, -0.93537087292458e-25),
)
I = np.array([t[0] for t in _TABLE2], dtype=float)
J = np.array([t[1] for t in _TABLE2], dtype=float)
N = np.array([t[2] for t in _TABLE2])


def _terms(p: float, t: float) -> tuple[float, float]:
    pi, tau = p / P_STAR, T_STAR / t
    return 7.1 - pi, tau - 1.222


def gamma_pi(p: float, t: float) -> float:
    a, b = _terms(p, t)
    return float(np.sum(-N * I * a ** (I - 1) * b**J))


def gamma_tau(p: float, t: float) -> float:
    a, b = _terms(p, t)
    return float(np.sum(N * a**I * J * b ** (J - 1)))


def gamma_tautau(p: float, t: float) -> float:
    a, b = _terms(p, t)
    return float(np.sum(N * a**I * J * (J - 1) * b ** (J - 2)))


def specific_volume(p: float, t: float) -> float:
    """m^3/kg at pressure p (Pa) and temperature t (K)."""
    return (p / P_STAR) * gamma_pi(p, t) * R * t / p


def density(p: float, t: float) -> float:
    """kg/m^3."""
    return 1.0 / specific_volume(p, t)


def enthalpy(p: float, t: float) -> float:
    """J/kg."""
    return (T_STAR / t) * gamma_tau(p, t) * R * t


def heat_capacity(p: float, t: float) -> float:
    """Isobaric, J/(kg K)."""
    tau = T_STAR / t
    return -(tau**2) * gamma_tautau(p, t) * R


def temperature(p: float, h: float, guess: float = 580.0) -> float:
    """The temperature (K) of liquid water with enthalpy h (J/kg) at p (Pa)."""
    t = guess
    for _ in range(50):
        step = (enthalpy(p, t) - h) / heat_capacity(p, t)
        t -= step
        if abs(step) < 1e-10:
            return t
    raise ValueError(f"no temperature for h = {h:g} J/kg at {p:g} Pa")


def gamma_pipi(p: float, t: float) -> float:
    a, b = _terms(p, t)
    return float(np.sum(N * I * (I - 1) * a ** (I - 2) * b**J))


def gamma_pitau(p: float, t: float) -> float:
    a, b = _terms(p, t)
    return float(np.sum(-N * I * a ** (I - 1) * J * b ** (J - 1)))


def heat_capacity_v(p: float, t: float) -> float:
    """Isochoric, J/(kg K) (IF97 Table 3):
    c_v / R = -tau^2 gamma_tautau + (gamma_pi - tau gamma_pitau)^2 / gamma_pipi."""
    tau = T_STAR / t
    g_pi, g_pipi = gamma_pi(p, t), gamma_pipi(p, t)
    return R * (-(tau**2) * gamma_tautau(p, t)
                + (g_pi - tau * gamma_pitau(p, t)) ** 2 / g_pipi)  # fmt: skip


def drho_dp(p: float, t: float) -> float:
    """(d rho / d p)_T in kg/(m^3 Pa): v = gamma_pi R T / p*, so
    dv/dp = gamma_pipi R T / p*^2 and d rho / dp = -rho^2 dv/dp."""
    rho = density(p, t)
    return -(rho**2) * gamma_pipi(p, t) * R * t / P_STAR**2


# ---------------------------------------------------------------- transport properties
# IAPWS R12-08, "Release on the IAPWS Formulation 2008 for the Viscosity of
# Ordinary Water Substance" (free: https://iapws.org/technical-guidance/
# release/viscosity.download), and IAPWS R15-11, "... 2011 for the Thermal
# Conductivity of Ordinary Water Substance" (free: https://iapws.org/
# technical-guidance/release/ThCond.download), both in their form for
# industrial use (section 3): density and heat capacities from IF97, and no
# critical enhancement of the viscosity (mu_2 = 1).
T_REF = 647.096  # K
RHO_REF = 322.0  # kg/m^3
P_REF = 22.064e6  # Pa
R_COND = 461.51805  # J/(kg K), Eq. (6) of R15-11 (used for c_p bar there)

# R12-08 Table 1 (H_i) and Table 2 (H_ij, the others are zero)
_H0 = (1.67752, 2.20462, 0.6366564, -0.241605)
_H1 = {
    (0, 0): 5.20094e-1, (1, 0): 8.50895e-2, (2, 0): -1.08374, (3, 0): -2.89555e-1,
    (0, 1): 2.22531e-1, (1, 1): 9.99115e-1, (2, 1): 1.88797, (3, 1): 1.26613,
    (5, 1): 1.20573e-1, (0, 2): -2.81378e-1, (1, 2): -9.06851e-1,
    (2, 2): -7.72479e-1, (3, 2): -4.89837e-1, (4, 2): -2.57040e-1,
    (0, 3): 1.61913e-1, (1, 3): 2.57399e-1, (0, 4): -3.25372e-2, (3, 4): 6.98452e-2,
    (4, 5): 8.72102e-3, (3, 6): -4.35673e-3, (5, 6): -5.93264e-4,
}  # fmt: skip
# R15-11 Table 1 (L_k), Table 2 (L_ij: rows i = 0..4, columns j = 0..5)
_L0 = (2.443221e-3, 1.323095e-2, 6.770357e-3, -3.454586e-3, 4.096266e-4)
_L1 = (
    (1.60397357, -0.646013523, 0.111443906, 0.102997357, -0.0504123634, 0.00609859258),
    (2.33771842, -2.78843778, 1.53616167, -0.463045512, 0.0832827019, -0.00719201245),
    (2.19650529, -4.54580785, 3.55777244, -1.40944978, 0.275418278, -0.0205938816),
    (-1.21051378, 1.60812989, -0.621178141, 0.0716373224, 0.0, 0.0),
    (-2.7203370, 4.57586331, -3.18369245, 1.1168348, -0.19268305, 0.012913842),
)
# R15-11 Table 3 (critical region) and Table 6 (A_ij of zeta(T_R, rho), Eq. 25)
_LAMBDA, _QD, _NU, _GAMMA, _XI0, _GAMMA0, _TR = (
    177.8514,
    1 / 0.40,
    0.630,
    1.239,
    0.13,
    0.06,
    1.5,
)
_A = (
    (6.53786807199516, 6.52717759281799, 5.35500529896124, 1.55225959906681, 1.11999926419994),
    (-5.61149954923348, -6.30816983387575, -3.96415689925446, 0.464621290821181, 0.595748562571649),
    (3.39624167361325, 8.08379285492595, 8.91990208918795, 8.93237374861479, 9.88952565078920),
    (-2.27492629730878, -9.82240510197603, -12.0338729505790, -11.0321960061126, -10.3255051147040),
    (10.2631854662709, 12.1358413791395, 9.19494865194302, 6.16780999933360, 4.66861294457414),
    (1.97815050331519, -5.54349664571295, -2.16866274479712, -0.965458722086812, -0.503243546373828),
)  # fmt: skip
_A_EDGES = (0.310559006, 0.776397516, 1.242236025, 1.863354037)


def viscosity(t: float, rho: float) -> float:
    """Pa s, R12-08 Eq. (10) with mu_2 = 1: mu = mu_0(T) mu_1(T, rho) mu*."""
    tb, rb = t / T_REF, rho / RHO_REF
    mu0 = 100.0 * math.sqrt(tb) / sum(h / tb**i for i, h in enumerate(_H0))
    s = sum(h * (1 / tb - 1) ** i * (rb - 1) ** j for (i, j), h in _H1.items())
    return 1e-6 * mu0 * math.exp(rb * s)


def _zeta_ref(rb: float) -> float:
    """zeta(T_R, rho) of R15-11 Eq. (25), dimensionless."""
    j = sum(rb > edge for edge in _A_EDGES)
    return 1.0 / sum(_A[i][j] * rb**i for i in range(6))


def conductivity(p: float, t: float) -> float:
    """W/(m K), R15-11 industrial formulation: lambda = lambda_0 lambda_1 +
    lambda_2, with the critical enhancement lambda_2 from IF97 properties."""
    rho = density(p, t)
    tb, rb = t / T_REF, rho / RHO_REF
    lam0 = math.sqrt(tb) / sum(c / tb**k for k, c in enumerate(_L0))
    s = sum(
        _L1[i][j] * (1 / tb - 1) ** i * (rb - 1) ** j
        for i in range(5)
        for j in range(6)
    )
    lam1 = math.exp(rb * s)
    # critical enhancement, Eqs. (18)-(25)
    cp, cv = heat_capacity(p, t), heat_capacity_v(p, t)
    zeta = drho_dp(p, t) * P_REF / RHO_REF
    dchi = max(rb * (zeta - _zeta_ref(rb) * _TR / tb), 0.0)
    xi = _XI0 * (dchi / _GAMMA0) ** (_NU / _GAMMA)  # nm
    y = _QD * xi
    lam2 = 0.0
    if y >= 1.2e-7:
        kappa = cp / cv
        z = 2 / (math.pi * y) * (((1 - 1 / kappa) * math.atan(y) + y / kappa)
                                 - (1 - math.exp(-1 / (1 / y + y * y / (3 * rb * rb)))))  # fmt: skip
        mub = viscosity(t, rho) / 1e-6
        lam2 = _LAMBDA * rb * (cp / R_COND) * tb / mub * z
    return 1e-3 * (lam0 * lam1 + lam2)
