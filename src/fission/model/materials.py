"""Number densities of a material at a state.

A compound of mass density rho (g/cm3) whose molecule contains n_e atoms of
element e has

    N_e = rho * N_A / M * n_e,        M = sum_e n_e * m_e,

atoms per cm3, where m_e is the mean atomic mass of the element. Isotopics
are given by atom fraction a_i (then m_e = sum_i a_i m_i) or by weight
fraction w_i (then 1/m_e = sum_i w_i / m_i and a_i = m_e w_i / m_i).
Multiplying by 1e-24 cm2/barn gives atoms per barn-cm, the unit of
the cross sections in barns.

Soluble boron is given in ppm by weight of the solution: a fraction
w_B = ppm * 1e-6 of the mass is boron, 1 - w_B the host (the convention of
the VERA specification, CASL-U-2012-0131-004 Problem 1, "Materials").

Atomic masses are AWR * m_n, with AWR the mass ratio to the neutron from the
nuclide's ENDF evaluation (MF1/MT451) and m_n the CODATA 2018 neutron mass.
"""

from __future__ import annotations

import re
from collections.abc import Callable
from functools import cache

from .spec import Isotopics, Material, State

NEUTRON_MASS_U = 1.00866491595  # CODATA 2018, unified atomic mass units
AVOGADRO = 6.02214076e23  # /mol, exact (SI 2019)
BARN = 1e-24  # cm^2

type MassFunction = Callable[[str], float]


@cache
def atomic_mass(nuclide: str) -> float:
    """Atomic mass (u) of a nuclide from its ENDF/B-VIII.1 evaluation."""
    from fission.data import EndfLibrary
    from fission.endf import Evaluation

    return (
        Evaluation.from_file(EndfLibrary().neutron(nuclide)).header.awr * NEUTRON_MASS_U
    )


def element_of(nuclide: str) -> str:
    """'U235' -> 'U', 'He4' -> 'He'."""
    m = re.fullmatch(r"([A-Z][a-z]?)\d+", nuclide)
    if m is None:
        raise ValueError(f"not a nuclide name: {nuclide!r}")
    return m.group(1)


def atom_fractions(
    iso: Isotopics, mass: MassFunction
) -> tuple[dict[str, float], float]:
    """(atom fraction of every isotope, mean atomic mass of the element)."""
    total = sum(iso.percent.values())
    frac = {n: p / total for n, p in iso.percent.items()}
    if iso.basis == "atom":
        mean = sum(a * mass(n) for n, a in frac.items())
        return frac, mean
    mean = 1.0 / sum(w / mass(n) for n, w in frac.items())
    return {n: mean * w / mass(n) for n, w in frac.items()}, mean


def number_densities(
    material: Material, state: State | None = None, mass: MassFunction = atomic_mass
) -> dict[str, float]:
    """Atoms per barn-cm of every nuclide of the material at the state."""
    if material.atom_densities is not None:
        return dict(material.atom_densities)
    assert material.formula is not None and material.isotopics is not None
    if material.density == "moderator":
        if state is None:
            raise ValueError(
                f"material {material.name!r} takes its density from a state"
            )
        rho = state.moderator_density
    else:
        assert material.density is not None
        rho = material.density
    w_boron = 0.0
    if material.soluble_boron is not None:
        if state is None:
            raise ValueError(f"material {material.name!r} takes its boron from a state")
        w_boron = state.boron_ppm * 1e-6

    out: dict[str, float] = {}

    def add(nuclide: str, n: float) -> None:
        out[nuclide] = out.get(nuclide, 0.0) + n

    elements = {
        e: atom_fractions(material.isotopics[e], mass) for e in material.formula
    }
    molar_mass = sum(material.formula[e] * elements[e][1] for e in material.formula)
    molecules = rho * (1.0 - w_boron) * AVOGADRO / molar_mass * BARN
    for e, count in material.formula.items():
        for nuclide, a in elements[e][0].items():
            add(nuclide, molecules * count * a)
    if material.soluble_boron is not None and w_boron > 0.0:
        frac, mean = atom_fractions(material.soluble_boron, mass)
        atoms = rho * w_boron * AVOGADRO / mean * BARN
        for nuclide, a in frac.items():
            add(nuclide, atoms * a)
    return out
