"""The neutron spectrum and k-infinity of a homogeneous mixture.

In an infinite homogeneous medium there is no leakage and no spatial
dependence, so the multigroup transport equation reduces to a balance per
group g (group 0 the fastest):

    Sigma_t,g phi_g = sum_g' Sigma_s,g'->g phi_g' + chi_g / k  sum_g' nu Sigma_f,g' phi_g'

or (T - S^T) phi = (1/k) F phi with T = diag(Sigma_t) and F = chi (nu
Sigma_f)^T. k-infinity is the largest eigenvalue of (T - S^T)^-1 F. With a
few hundred groups this is solved directly; the transport solver of the
pin cell will need iterations instead. Each fissionable nuclide brings its
own chi, so F is a sum of such outer products.

Only P0 scattering matters here (no directions in an infinite medium); the
inelastic matrix carries the (n,xn) multiplicities.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from fission.model import Component
from fission.processing.groupr import GroupLibrary

THERMAL_CUTOFF = 0.625  # eV, the cadmium-like boundary of two-group theory


@dataclass
class Mixture:
    """Macroscopic multigroup data (1/cm) of a homogeneous mixture."""

    components: list[Component]
    libraries: list[GroupLibrary]

    @property
    def groups(self):
        return self.libraries[0].groups

    def _sum(self, f) -> np.ndarray:
        return sum(
            c.density * f(lib)
            for c, lib in zip(self.components, self.libraries, strict=True)
        )

    @property
    def total(self) -> np.ndarray:
        return self._sum(lambda lib: lib.sigma[1])

    @property
    def scatter(self) -> np.ndarray:
        """P0 transfer matrix S[g, g'] from g to g', including (n,xn) yields."""
        return self._sum(lambda lib: lib.scatter[0])

    def reaction(self, mt: int) -> np.ndarray:
        return self._sum(lambda lib: lib.sigma.get(mt, np.zeros(lib.groups.size)))

    @property
    def nu_fission(self) -> np.ndarray:
        return self._sum(
            lambda lib: lib.nu_fission if lib.nu_fission is not None else 0.0
        )

    @property
    def removal(self) -> np.ndarray:
        """Net loss per collision: Sigma_t minus all neutrons re-emitted by scattering
        (absorption, minus the extra neutrons of (n,xn))."""
        return self.total - self.scatter.sum(axis=1)

    def fission_matrix(self) -> np.ndarray:
        """F = sum over fissionable components of chi_i (nu Sigma_f,i)^T."""
        g = self.groups.size
        f = np.zeros((g, g))
        for c, lib in zip(self.components, self.libraries, strict=True):
            if lib.nu_fission is not None and lib.chi is not None:
                f += np.outer(lib.chi, c.density * lib.nu_fission)
        return f


@dataclass
class Solution:
    k: float
    flux: (
        np.ndarray
    )  # group fluxes (integrated over each group), normalized: sum nu F phi = 1
    mixture: Mixture

    def rates(self) -> dict[str, float]:
        """Reaction rates per fission neutron, split at the thermal cutoff."""
        m, phi = self.mixture, self.flux
        th = np.array(
            [
                m.groups.bounds(g)[1] <= THERMAL_CUTOFF * (1 + 1e-9)
                for g in range(m.groups.size)
            ]
        )
        fuel = [c.region == "fuel" for c in m.components]
        absorption_fuel = sum(
            c.density * (lib.sigma.get(102, 0) + lib.sigma.get(18, 0))
            for c, lib, in_fuel in zip(m.components, m.libraries, fuel, strict=True)
            if in_fuel
        )
        nf, rem = m.nu_fission * phi, m.removal * phi
        return {
            "production": float(nf.sum()),
            "production_thermal": float(nf[th].sum()),
            "removal": float(rem.sum()),
            "removal_thermal": float(rem[th].sum()),
            "fuel_absorption_thermal": float((absorption_fuel * phi)[th].sum()),
        }

    def four_factors(self) -> dict[str, float]:
        """k-infinity = eta f p epsilon, with the factors defined from the rates:

        epsilon = all fission neutrons / those from thermal fission (fast fission factor)
        p       = thermal removals / all removals (resonance escape probability)
        f       = thermal absorption in fuel / all thermal removals (thermal utilization)
        eta     = thermal fission neutrons / thermal absorption in fuel
        """
        r = self.rates()
        eps = r["production"] / r["production_thermal"]
        p = r["removal_thermal"] / r["removal"]
        f = r["fuel_absorption_thermal"] / r["removal_thermal"]
        eta = r["production_thermal"] / r["fuel_absorption_thermal"]
        return {
            "epsilon": eps,
            "p": p,
            "f": f,
            "eta": eta,
            "product": eps * p * f * eta,
        }


def solve(mixture: Mixture, use_rust: bool = True) -> Solution:
    """k-infinity and the fundamental-mode spectrum.

    Python: a direct dense eigen-solve (numpy). Rust: LU factorization once,
    then power iteration, the pattern the transport solver will use.
    """
    f = mixture.fission_matrix()
    if use_rust:
        from fission import _rust

        k, phi, _ = _rust.infinite_medium(
            mixture.total.tolist(), mixture.scatter.ravel().tolist(), f.ravel().tolist()
        )
        phi = np.abs(phi)
        phi /= float(mixture.nu_fission @ phi)
        return Solution(float(k), phi, mixture)
    a = np.diag(mixture.total) - mixture.scatter.T
    m = np.linalg.solve(a, f)
    values, vectors = np.linalg.eig(m)
    i = int(np.argmax(values.real))
    k = float(values[i].real)
    phi = np.abs(vectors[:, i].real)
    phi /= float(mixture.nu_fission @ phi)
    return Solution(k, phi, mixture)
