"""B1 fundamental mode of a homogenized lattice (step C3).

The lattice calculation assumes an infinite array of identical cells: no
leakage. In a core the same assembly leaks, and leakage hardens the
spectrum (thermal neutrons travel less far but fast ones leak more easily
per collision). The usual correction [Hebert, Applied Reactor Physics] is to
solve the homogenized medium in its fundamental mode,

    phi(x, E, mu) = psi(E, mu) exp(i B x),

the shape of the flux in a large bare reactor, and to find the critical
buckling B^2 that makes it exactly critical (k = 1). With P1 scattering,
expanding psi in Legendre moments phi_g (flux) and i J_g (current) and
closing the angular equation exactly (not with the diffusion
approximation) gives the B1 equations for every group g:

    Sigma_t,g phi_g + B J_g = sum_g' Sigma_s0,g'->g phi_g' + chi_g F / k,
    3 alpha_g Sigma_t,g J_g - B phi_g = 3 sum_g' Sigma_s1,g'->g J_g',

with F = sum_g nu Sigma_f,g phi_g and, for x = B / Sigma_t,g,

    alpha(x) = x^2 arctan(x) / (3 (x - arctan(x)))   (B^2 > 0),
    alpha(x) = x^2 artanh(x) / (3 (artanh(x) - x))   (B^2 < 0, x = |B| / Sigma_t),

which tends to 1 as B -> 0. With alpha = 1 these are the P1 equations, and
J_g = D_g B phi_g with D_g = 1 / (3 (Sigma_t,g - Sigma_s1,g->g)) when the
P1 matrix is diagonal: B1 *computes* the diffusion coefficients,
D_g = J_g / (B phi_g), including the transport across groups.

Eliminating J: J = B M phi with M = (3 diag(alpha Sigma_t) - 3 S1^T)^-1,
then (diag(Sigma_t) - S0^T + B^2 M) phi = chi F / k, so that
k(B^2) = nu Sigma_f . A(B^2)^-1 chi. The critical buckling solves k = 1.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .homogenize import FineConstants


def alpha(b2: float, sigma_t: np.ndarray) -> np.ndarray:
    """The B1 correction factor of every group (1 in the P1 approximation)."""
    if b2 == 0.0:
        return np.ones_like(sigma_t)
    x = np.sqrt(abs(b2)) / sigma_t
    out = np.empty_like(x)
    small = x < 1e-3
    sign = 1.0 if b2 > 0 else -1.0
    # series: alpha = 1 + 4 x^2 / 15 (B^2 > 0), 1 - 4 x^2 / 15 (B^2 < 0)
    out[small] = 1.0 + sign * 4.0 * x[small] ** 2 / 15.0
    xs = x[~small]
    if b2 > 0:
        out[~small] = xs**2 * np.arctan(xs) / (3.0 * (xs - np.arctan(xs)))
    else:
        if np.any(xs >= 1.0):
            raise ValueError("negative B^2 larger than Sigma_t^2: no real mode")
        out[~small] = xs**2 * np.arctanh(xs) / (3.0 * (np.arctanh(xs) - xs))
    return out


@dataclass
class Mode:
    """The B1 fundamental mode at one buckling."""

    b2: float  # 1/cm^2
    k: float
    flux: np.ndarray  # phi_g, normalized to F = 1
    current: np.ndarray  # J_g / B (cm), so that D_g = current_g / flux_g

    @property
    def diffusion(self) -> np.ndarray:
        return self.current / self.flux


def mode(c: FineConstants, b2: float) -> Mode:
    """Solve the B1 equations at a given buckling."""
    s0, s1 = c.scatter[0], c.scatter[1]
    a = alpha(b2, c.total)
    m = np.linalg.inv(3.0 * np.diag(a * c.total) - 3.0 * s1.T)
    operator = np.diag(c.total) - s0.T + b2 * m
    shape = np.linalg.solve(operator, c.chi)  # the flux for a source chi F / k = chi
    k = float(c.nu_fission @ shape)  # F of that flux: k = F / (F / k)
    flux = shape / k  # normalized to F = 1
    return Mode(b2, k, flux, m @ flux)


def critical(c: FineConstants, tol: float = 1e-12, max_iterations: int = 100) -> Mode:
    """The buckling that makes the fundamental mode critical (k = 1).

    k(B^2) decreases with B^2. The secant method on k - 1 starts from the
    one-group estimate B^2 = (k_inf - 1) / M^2, with the migration area
    M^2 = <D> / <Sigma_a> averaged over the infinite-medium spectrum."""
    m0 = mode(c, 0.0)
    migration = float(m0.current.sum()) / float(c.absorption @ m0.flux)
    x0, f0 = 0.0, m0.k - 1.0
    x1 = (m0.k - 1.0) / migration
    f1 = mode(c, x1).k - 1.0
    for _ in range(max_iterations):
        if abs(f1) < tol:
            return mode(c, x1)
        x2 = x1 - f1 * (x1 - x0) / (f1 - f0)
        x0, f0 = x1, f1
        x1, f1 = x2, mode(c, x2).k - 1.0
    raise RuntimeError("critical buckling search did not converge")
