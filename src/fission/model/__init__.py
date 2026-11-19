"""Declarative reactor models: materials -> pins -> lattices -> assemblies -> core.

from fission.model import load_model, homogenize
model = load_model("vera-p1")            # models/vera-p1.toml
homogenize(model, "cell", "1A")          # area-weighted components
"""

from .homogenize import Component, homogenize, material_areas
from .load import MODELS_DIR, ModelError, load_model, parse
from .materials import atomic_mass, number_densities
from .spec import (
    Assembly,
    AxialSegment,
    Core,
    Isotopics,
    Lattice,
    Material,
    Model,
    Pin,
    Reference,
    State,
)

__all__ = [
    "MODELS_DIR",
    "Assembly",
    "AxialSegment",
    "Component",
    "Core",
    "Isotopics",
    "Lattice",
    "Material",
    "Model",
    "ModelError",
    "Pin",
    "Reference",
    "State",
    "atomic_mass",
    "homogenize",
    "load_model",
    "material_areas",
    "number_densities",
    "parse",
]
