"""Declarative reactor models: materials -> pins -> lattices -> assemblies -> core.

from fission.model import load_model, homogenize
model = load_model("vera-p1")            # models/vera-p1.toml
homogenize(model, "cell", "1A")          # area-weighted components
"""

from .homogenize import Component, homogenize, material_areas
from .load import MODELS_DIR, ModelError, load_model, parse
from .materials import atomic_mass, number_densities
from .multigroup import MacroXS, load_multigroup
from .spec import (
    BRANCH_AXES,
    Assembly,
    AxialSegment,
    Branches,
    Core,
    CoreThermal,
    Isotopics,
    Lattice,
    Material,
    Model,
    Pin,
    Reference,
    State,
)

__all__ = [
    "BRANCH_AXES",
    "MODELS_DIR",
    "Assembly",
    "AxialSegment",
    "Branches",
    "Component",
    "Core",
    "CoreThermal",
    "Isotopics",
    "Lattice",
    "MacroXS",
    "Material",
    "Model",
    "ModelError",
    "Pin",
    "Reference",
    "State",
    "atomic_mass",
    "homogenize",
    "load_model",
    "load_multigroup",
    "material_areas",
    "number_densities",
    "parse",
]
