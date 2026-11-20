"""Flux-volume homogenization of a solved lattice (step C3).

A nodal code sees each assembly (here: each pin cell) as one homogeneous
node. Its cross sections must reproduce the reaction rates of the
heterogeneous transport solution: for every reaction x and fine group g,

    Sigma_x,g = sum_r V_r Sigma_x,r,g phi_r,g / sum_r V_r phi_r,g,

summed over the flat-source regions r of the cell. The scattering matrices
are weighted with the scalar flux of the departure group, also the P1
matrix (in an infinite lattice the net current is zero, so there is no
better weight). chi is weighted with each region's fission production, so
that the homogeneous cell emits the same fission neutrons.

The totals are the plain ones, without the transport correction used by
MOC: B1 (``b1.py``) treats the P1 scattering explicitly. The removal
Sigma_t - sum_g' Sigma_s0,g->g' is the same with or without the correction,
so the homogeneous infinite-medium eigenvalue equals the MOC k-infinity.

**Discontinuity factors.** The homogeneous node has a flat flux; the real
one dips in the fuel and peaks in the water. The ratio of the heterogeneous
flux on the node surface to the homogeneous one (the average),
f_g = phi_surface,g / phi_g, lets the nodal code reproduce the interface
currents between different assemblies [Hebert]. The surface flux comes from
the boundary angular fluxes of MOC (``moc.surface_flux``), averaged over
the four reflective sides.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from fission.processing.groups import GroupStructure
from fission.transport.moc import surface_flux

from .kappa import recoverable_energy


@dataclass
class FineConstants:
    """Homogenized cell constants in the library's groups (1/cm)."""

    name: str
    groups: GroupStructure
    flux: np.ndarray  # cell-average scalar flux phi_g (MOC normalization)
    total: np.ndarray
    scatter: np.ndarray  # (2, G, G): P0, P1 from g to g'
    nu_fission: np.ndarray
    fission: np.ndarray
    kappa_fission: np.ndarray  # eV / cm
    chi: np.ndarray
    surface: np.ndarray | None = None  # side-averaged surface flux phi_s,g

    @property
    def absorption(self) -> np.ndarray:
        """Removal minus out-scattering: Sigma_t - sum_g' Sigma_s0,g->g'."""
        return self.total - self.scatter[0].sum(axis=1)

    def k_infinity(self) -> float:
        """Eigenvalue of the infinite homogeneous medium with these constants."""
        a = np.diag(self.total) - self.scatter[0].T
        return float(self.nu_fission @ np.linalg.solve(a, self.chi))


def homogenize(result, name: str = "cell") -> FineConstants:
    """Constants of the whole cell of a ``heterogeneous.Result``."""
    sol = result.solution
    problem = sol.problem
    labels = np.array(problem.geometry.materials())
    names = sorted(set(labels.tolist()))  # MOC materials, in Problem.xs order
    groups = next(iter(next(iter(result.libraries.values())).values())).groups
    size = groups.size
    upper, lower = np.asarray(groups.boundaries[:-1]), np.asarray(groups.boundaries[1:])
    midpoint = np.sqrt(upper * lower)
    volume = problem.areas.sum()
    flux_total = (problem.areas[:, None] * sol.flux).sum(axis=0)  # sum_r V_r phi_r
    rates = {k: np.zeros(size) for k in ("total", "nu", "fission", "kappa")}
    scatter = np.zeros((2, size, size))
    chi = np.zeros(size)
    production = 0.0
    for index, label in enumerate(names):
        mask = labels == label
        flux = (problem.areas[mask, None] * sol.flux[mask]).sum(axis=0)
        region_production = 0.0
        for n, lib in result.libraries[label].items():
            d = result.densities[label][n]
            rates["total"] += d * lib.sigma[1] * flux
            s = lib.scatter
            for order in range(min(2, s.shape[0])):
                scatter[order] += d * s[order] * flux[:, None]
            if 18 in lib.sigma:
                kappa = np.array([recoverable_energy(n)(e) for e in midpoint])
                rates["fission"] += d * lib.sigma[18] * flux
                rates["kappa"] += d * kappa * lib.sigma[18] * flux
            if lib.nu_fission is not None:
                nu = d * lib.nu_fission * flux
                rates["nu"] += nu
                region_production += float(nu.sum())
        chi += problem.xs[index].chi * region_production
        production += region_production
    with np.errstate(invalid="ignore", divide="ignore"):
        per = {
            k: np.where(flux_total > 0, v / flux_total, 0.0) for k, v in rates.items()
        }
        scatter = np.where(flux_total[None, :, None] > 0,
                           scatter / flux_total[None, :, None], 0.0)  # fmt: skip
    surface = None
    if sol.psi is not None:
        sides = surface_flux(problem, sol.psi)
        if sides:
            surface = np.mean(list(sides.values()), axis=0)
    return FineConstants(
        name, groups, flux_total / volume, per["total"], scatter, per["nu"],
        per["fission"], per["kappa"], chi / production, surface,
    )  # fmt: skip
