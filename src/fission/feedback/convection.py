"""Heat transfer from the clad to the water.

Turbulent forced convection in single-phase water, with the Dittus-Boelter
correlation in McAdams' form for a wall hotter than the fluid,

    Nu = h D_h / lambda = 0.023 Re^0.8 Pr^0.4,
    Re = G D_h / mu,   Pr = c_p mu / lambda,

valid for 1e4 <= Re <= 1.2e5 and 0.7 <= Pr <= 160, and Presser's factor
for a square rod bundle, Nu = psi Nu_tube with

    psi = 0.9217 + 0.1478 P/D - 0.1130 exp(-7 (P/D - 1))   (1.05 <= P/D <= 1.9),

as documented in the MOOSE subchannel and Navier-Stokes modules (free:
https://github.com/idaholab/moose, modules/subchannel/doc/content/source/
scmclosures/SCMHTCDittusBoelter.md and modules/navier_stokes/doc/content/
source/functormaterials/DittusBoelterFunctorMaterial.md; after Dittus and
Boelter (1930), McAdams (1954), Presser (1967)). The water properties are
those of IAPWS (``water``), at the bulk water temperature.
"""

from __future__ import annotations

import math

from . import water


def presser(pitch_over_diameter: float) -> float:
    """Bundle correction psi for a square array."""
    x = pitch_over_diameter
    return 0.9217 + 0.1478 * x - 0.1130 * math.exp(-7.0 * (x - 1.0))


def film_coefficient(p: float, t: float, mass_flux: float, diameter: float,
                     pitch_over_diameter: float | None = None) -> float:  # fmt: skip
    """h in W/(m^2 K) for water at p (Pa), t (K), mass flux G (kg/(m^2 s))
    in a channel of hydraulic diameter D_h (m)."""
    rho = water.density(p, t)
    mu = water.viscosity(t, rho)
    k = water.conductivity(p, t)
    re = mass_flux * diameter / mu
    pr = water.heat_capacity(p, t) * mu / k
    nu = 0.023 * re**0.8 * pr**0.4
    if pitch_over_diameter is not None:
        nu *= presser(pitch_over_diameter)
    return nu * k / diameter


def reynolds_prandtl(p: float, t: float, mass_flux: float, diameter: float,
                     ) -> tuple[float, float]:  # fmt: skip
    rho = water.density(p, t)
    mu = water.viscosity(t, rho)
    return mass_flux * diameter / mu, water.heat_capacity(
        p, t
    ) * mu / water.conductivity(p, t)
