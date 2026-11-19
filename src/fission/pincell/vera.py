"""VERA Core Physics Benchmark Progression Problem 1: a PWR pin cell.

Source: A. T. Godfrey, "VERA Core Physics Benchmark Progression Problem
Specifications", CASL-U-2012-0131-004 (Rev. 4, 2014): Table 1 (rod
geometry), Table P1-1 (states), Table P1-4 (isotopic number densities) and
Table P1-5 (CE KENO-VI reference eigenvalues, ENDF/B-VII.0).
https://corephysics.com/docs/CASL-U-2012-0131-004.pdf

Number densities are in atoms per barn-cm. Isotope names follow the ENDF
library ('U235', 'Zr90', ...).
"""

from __future__ import annotations

import math
from dataclasses import dataclass

# Table 1: fuel rod geometry (cm)
PELLET_RADIUS = 0.4096
CLAD_INNER_RADIUS = 0.418
CLAD_OUTER_RADIUS = 0.475
PITCH = 1.26

# Table P1-4: reference mixing table
FUEL = {
    "U234": 6.11864e-06, "U235": 7.18132e-04, "U236": 3.29861e-06,
    "U238": 2.21546e-02, "O16": 4.57642e-02,
}  # fmt: skip
GAP = {"He4": 2.68714e-05}
CLAD = {
    "Zr90": 2.18865e-02, "Zr91": 4.77292e-03, "Zr92": 7.29551e-03,
    "Zr94": 7.39335e-03, "Zr96": 1.19110e-03,
    "Sn112": 4.68066e-06, "Sn114": 3.18478e-06, "Sn115": 1.64064e-06,
    "Sn116": 7.01616e-05, "Sn117": 3.70592e-05, "Sn118": 1.16872e-04,
    "Sn119": 4.14504e-05, "Sn120": 1.57212e-04, "Sn122": 2.23417e-05,
    "Sn124": 2.79392e-05,
    "Fe54": 8.68307e-06, "Fe56": 1.36306e-04, "Fe57": 3.14789e-06, "Fe58": 4.18926e-07,
    "Cr50": 3.30121e-06, "Cr52": 6.36606e-05, "Cr53": 7.21860e-06, "Cr54": 1.79686e-06,
    "Hf174": 3.54138e-09, "Hf176": 1.16423e-07, "Hf177": 4.11686e-07,
    "Hf178": 6.03806e-07, "Hf179": 3.01460e-07, "Hf180": 7.76449e-07,
}  # fmt: skip
MODERATOR_1A = {
    "O16": 2.48112e-02,
    "H1": 4.96224e-02,
    "B10": 1.07070e-05,
    "B11": 4.30971e-05,
}
MODERATOR_1B = {
    "O16": 2.20729e-02,
    "H1": 4.41459e-02,
    "B10": 9.52537e-06,
    "B11": 3.83408e-05,
}

# Nuclides whose ENDF/B-VIII.1 resonance data this code cannot process yet:
# Fe-54 and Fe-57 (R-matrix limited format, LRF = 7), Fe-58 (an unresolved-
# range form, LRF = 1), Cr-53 (Reich-Moore with negative J selecting the
# channel spin). Together 1.9e-5 /b-cm of the cladding (0.04% of its atoms);
# they are left out.
UNSUPPORTED = frozenset({"Fe54", "Fe57", "Fe58", "Cr53"})


@dataclass(frozen=True)
class Case:
    name: str
    fuel_temperature: float  # K
    moderator_temperature: float  # K (the clad and gap take this temperature)
    moderator: dict[str, float]
    keff: float  # CE KENO-VI reference (heterogeneous pin cell, ENDF/B-VII.0)
    keff_sigma: float


# Table P1-1 and Table P1-5
CASES = {
    "1A": Case("1A", 565.0, 565.0, MODERATOR_1A, 1.187038, 0.000054),
    "1B": Case("1B", 600.0, 600.0, MODERATOR_1B, 1.182149, 0.000068),
    "1C": Case("1C", 900.0, 600.0, MODERATOR_1B, 1.171722, 0.000072),
    "1D": Case("1D", 1200.0, 600.0, MODERATOR_1B, 1.162603, 0.000071),
}


def region_areas() -> dict[str, float]:
    """Cross-sectional areas (cm^2) of the four regions of the square cell."""
    fuel = math.pi * PELLET_RADIUS**2
    gap = math.pi * CLAD_INNER_RADIUS**2 - fuel
    clad = math.pi * (CLAD_OUTER_RADIUS**2 - CLAD_INNER_RADIUS**2)
    moderator = PITCH**2 - math.pi * CLAD_OUTER_RADIUS**2
    return {"fuel": fuel, "gap": gap, "clad": clad, "moderator": moderator}


@dataclass(frozen=True)
class Component:
    """One nuclide of the homogenized cell, from one region."""

    nuclide: str
    region: str
    temperature: float
    density: float  # cell-averaged atoms per barn-cm
    thermal: str  # "free-gas" or an S(alpha, beta) material


def homogenized(case: Case, include_unsupported: bool = False) -> list[Component]:
    """Volume-weighted number densities of every nuclide in every region.

    The same nuclide in two regions (O-16 in fuel and water) stays two
    components: they have different temperatures and scattering laws.
    """
    areas = region_areas()
    cell = PITCH**2
    regions = {
        "fuel": (FUEL, case.fuel_temperature),
        "gap": (GAP, case.moderator_temperature),
        "clad": (CLAD, case.moderator_temperature),
        "moderator": (case.moderator, case.moderator_temperature),
    }
    out = []
    for region, (composition, temperature) in regions.items():
        for nuclide, density in composition.items():
            if nuclide in UNSUPPORTED and not include_unsupported:
                continue
            thermal = (
                "HinH2O" if (nuclide == "H1" and region == "moderator") else "free-gas"
            )
            out.append(
                Component(
                    nuclide,
                    region,
                    temperature,
                    density * areas[region] / cell,
                    thermal,
                )
            )
    return out
