"""Smear a 2D lattice into one homogeneous mixture (volume weighting).

A lattice is infinitely tall, so volume fractions are area fractions: the
homogenized density of a nuclide is sum_r N_r A_r / A over the regions r
that contain it. This conserves the number of atoms in the slice. The same
nuclide in two materials (oxygen in fuel and water) stays two components,
because they have different temperatures and scattering laws.
"""

from __future__ import annotations

from dataclasses import dataclass

from .materials import MassFunction, atomic_mass, number_densities
from .spec import Lattice, Model, State


@dataclass(frozen=True)
class Component:
    """One nuclide of a homogenized mixture, from one material."""

    nuclide: str
    region: str  # the material it came from
    temperature: float  # K
    density: float  # slice-averaged atoms per barn-cm
    thermal: str  # "free-gas" or an S(alpha, beta) law


def material_areas(model: Model, lattice: Lattice) -> dict[str, float]:
    """Area (cm^2) covered by each material in one slice of the lattice."""
    if lattice.material is not None:
        return {lattice.material: 1.0}
    assert lattice.pins is not None and lattice.pitch is not None
    areas: dict[str, float] = {}
    for row in lattice.pins:
        for name in row:
            for material, area in model.pins[name].areas(lattice.pitch):
                areas[material] = areas.get(material, 0.0) + area
    return areas


def homogenize(
    model: Model, lattice: str, state: str, mass: MassFunction = atomic_mass
) -> list[Component]:
    """Area-weighted components of a lattice at a state, in material order."""
    lat = model.lattices[lattice]
    st: State = model.states[state]
    areas = material_areas(model, lat)
    total = sum(areas.values())
    out = []
    for name, area in areas.items():
        material = model.materials[name]
        temperature = st.temperature(material.temperature)
        for nuclide, n in number_densities(material, st, mass).items():
            out.append(
                Component(
                    nuclide,
                    name,
                    temperature,
                    n * area / total,
                    material.thermal_law(nuclide),
                )
            )
    return out
