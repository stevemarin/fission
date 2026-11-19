"""Macroscopic multigroup cross sections given directly as data.

Benchmarks such as C5G7 specify materials by their macroscopic group
constants (1/cm), not by nuclides. A model points at such a library with
``[model] multigroup = "file.toml"``, and a material uses one entry with
``macroscopic = "name"``. The library file holds, per material,

    transport (or total), absorption, scatter[g][g'] (g -> g', g = 0 fastest),
    and for fissile materials fission, nu and chi.

The transport solver uses ``transport`` as its total cross section when it is
given (transport-corrected, isotropic scattering), otherwise ``total``.
"""

from __future__ import annotations

import tomllib
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np


class LibraryError(ValueError):
    pass


@dataclass(frozen=True)
class MacroXS:
    """Group constants of one material (1/cm); group 0 is the fastest."""

    name: str
    total: np.ndarray  # the total used in transport (transport-corrected if given)
    absorption: np.ndarray
    scatter: np.ndarray  # [g, g'] from g to g'
    fission: np.ndarray
    nu_fission: np.ndarray
    chi: np.ndarray

    @property
    def groups(self) -> int:
        return len(self.total)

    @property
    def fissile(self) -> bool:
        return bool(np.any(self.nu_fission > 0))


def _vector(m: dict[str, Any], key: str, where: str, groups: int) -> np.ndarray:
    v = m.get(key)
    if not isinstance(v, list) or len(v) != groups:
        raise LibraryError(f"{where}.{key}: expected {groups} numbers")
    a = np.array(v, dtype=float)
    if np.any(a < 0):
        raise LibraryError(f"{where}.{key}: negative value")
    return a


def parse_multigroup(data: dict[str, Any], source: str = "") -> dict[str, MacroXS]:
    head = data.get("library", {})
    if head.get("format") != 1:
        raise LibraryError(f"{source}library.format: expected 1")
    groups = head.get("groups")
    if not isinstance(groups, int) or groups < 1:
        raise LibraryError(f"{source}library.groups: expected a positive integer")
    out = {}
    for name, m in data.get("materials", {}).items():
        where = f"{source}materials.{name}"
        total = _vector(m, "transport" if "transport" in m else "total", where, groups)
        absorption = _vector(m, "absorption", where, groups)
        scatter = np.array(m.get("scatter", []), dtype=float)
        if scatter.shape != (groups, groups) or np.any(scatter < 0):
            raise LibraryError(
                f"{where}.scatter: expected {groups} x {groups} non-negative"
            )
        if "fission" in m:
            fission = _vector(m, "fission", where, groups)
            nu_fission = fission * _vector(m, "nu", where, groups)
            chi = _vector(m, "chi", where, groups)
            if abs(chi.sum() - 1) > 1e-3:
                raise LibraryError(f"{where}.chi: sums to {chi.sum():g}, not 1")
            chi = chi / chi.sum()
        else:
            fission = nu_fission = chi = np.zeros(groups)
        out[name] = MacroXS(name, total, absorption, scatter, fission, nu_fission, chi)
    return out


def load_multigroup(path: Path) -> dict[str, MacroXS]:
    with open(path, "rb") as f:
        return parse_multigroup(tomllib.load(f), f"{path.name}: ")
