"""Condensation to few groups and the constants handed to the nodal code.

With a spectrum phi_g (and currents J_g / B, the diffusion "current per unit
buckling"), every few-group constant is the rate-preserving average over its
fine groups g in G:

    phi_G = sum_g phi_g,          Sigma_x,G = sum_g Sigma_x,g phi_g / phi_G,
    Sigma_s,G->G' = sum_{g in G} sum_{g' in G'} Sigma_s0,g->g' phi_g / phi_G,
    chi_G = sum_g chi_g,          D_G = sum_g (J_g / B) / phi_G.

D_G is current-weighted, so the leakage D_G B^2 phi_G of the few-group
diffusion equation is exactly the B1 leakage sum_g B J_g. With the critical
B1 spectrum the few-group constants therefore give k = 1 at the critical
buckling, and with the infinite-medium spectrum (B = 0) they reproduce the
lattice k-infinity: both are tests of the bookkeeping.
"""

from __future__ import annotations

import itertools
import json
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from .b1 import Mode
from .homogenize import FineConstants

THERMAL_CUT = 0.625  # eV: the usual two-group boundary for LWRs


@dataclass
class FewGroup:
    """Homogenized few-group constants (1/cm, D in cm, kappa Sigma_f in eV/cm)."""

    name: str
    boundaries: tuple[float, ...]  # eV, descending
    flux: np.ndarray  # phi_G of the weighting spectrum (normalized to F = 1)
    D: np.ndarray
    total: np.ndarray
    sigma_a: np.ndarray
    sigma_s: np.ndarray  # [G, G'] from G to G' (P0), diagonal included
    nu_sigma_f: np.ndarray
    sigma_f: np.ndarray
    kappa_sigma_f: np.ndarray
    chi: np.ndarray
    buckling: float  # B^2 of the weighting spectrum (0 = infinite medium)
    adf: np.ndarray = field(default_factory=lambda: np.ones(0))

    def __post_init__(self):
        if self.adf.size == 0:
            self.adf = np.ones(len(self.D))

    def k(self, b2: float = 0.0) -> float:
        """Eigenvalue of the homogeneous few-group diffusion problem at B^2."""
        a = np.diag(self.total + self.D * b2) - self.sigma_s.T
        return float(self.nu_sigma_f @ np.linalg.solve(a, self.chi))

    def material(self) -> dict:
        """Keyword arguments of ``nodal.Material`` (branch claude-nodal)."""
        return {
            "name": self.name,
            "D": self.D.tolist(),
            "sigma_a": self.sigma_a.tolist(),
            "nu_sigma_f": self.nu_sigma_f.tolist(),
            "sigma_s": self.sigma_s.tolist(),
            "chi": self.chi.tolist(),
            "kappa_sigma_f": self.kappa_sigma_f.tolist(),
            "adf": self.adf.tolist(),
        }

    def write(self, path: Path, **metadata) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        data = {
            "material": self.material(),
            "boundaries_eV": list(self.boundaries),
            "buckling_cm-2": self.buckling,
            "k_infinity": self.k(0.0),
            "flux": self.flux.tolist(),
            "sigma_f": self.sigma_f.tolist(),
            "units": "1/cm; D in cm; kappa_sigma_f in eV/cm",
            **metadata,
        }
        path.write_text(json.dumps(data, indent=1))


def condense(
    c: FineConstants, weight: Mode, boundaries: tuple[float, ...] = (THERMAL_CUT,)
) -> FewGroup:
    """Few-group constants of `c` with the spectrum and currents of `weight`.

    boundaries: the inner few-group boundaries (eV); each must coincide with
    a fine-group boundary."""
    fine = np.asarray(c.groups.boundaries)
    edges = [fine[0], *sorted(boundaries, reverse=True), fine[-1]]
    cuts = []
    for e in edges:
        i = int(np.argmin(np.abs(fine - e)))
        if not np.isclose(fine[i], e, rtol=1e-6):
            raise ValueError(f"{e} eV is not a boundary of {c.groups.name}")
        cuts.append(i)
    owner = np.zeros(c.groups.size, dtype=int)
    for big, (lo, hi) in enumerate(itertools.pairwise(cuts)):
        owner[lo:hi] = big
    size = len(cuts) - 1
    phi = weight.flux
    flux = np.bincount(owner, phi, size)

    def average(x: np.ndarray) -> np.ndarray:
        return np.bincount(owner, x * phi, size) / flux

    sigma_s = np.zeros((size, size))
    for g in range(c.groups.size):
        sigma_s[owner[g]] += np.bincount(owner, c.scatter[0][g], size) * phi[g]
    sigma_s /= flux[:, None]
    return FewGroup(
        c.name,
        tuple(float(fine[i]) for i in cuts),
        flux,
        np.bincount(owner, weight.current, size) / flux,
        average(c.total),
        average(c.absorption),
        sigma_s,
        average(c.nu_fission),
        average(c.fission),
        average(c.kappa_fission),
        np.bincount(owner, c.chi, size),
        weight.b2,
    )
