"""Macroscopic multigroup cross sections given directly as data.

Benchmarks such as C5G7 specify materials by their macroscopic group
constants (1/cm), not by nuclides. A model points at such a library with
``[model] multigroup = "file.toml"``, and a material uses one entry with
``macroscopic = "name"``. The library file holds, per material,

    transport (or total), absorption, scatter[g][g'] (g -> g', g = 0 fastest),
    and for fissile materials fission, nu and chi (or nu_fission and chi when
    a benchmark gives only the product). Diffusion benchmarks also give the
    diffusion coefficient D per group (``diffusion``, cm); without it, a
    diffusion solver uses D = 1 / (3 Sigma_tr).

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
    diffusion: np.ndarray | None = None  # D (cm), if the library gives it
    kappa_fission: np.ndarray | None = None  # recoverable energy x Sigma_f (eV/cm)
    adf: np.ndarray | None = None  # assembly discontinuity factors (1 if None)

    @property
    def diffusion_coefficient(self) -> np.ndarray:
        """D per group: as given, or 1 / (3 Sigma_tr) from the transport total."""
        if self.diffusion is not None:
            return self.diffusion
        return 1.0 / (3.0 * self.total)

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
        if "fission" in m or "nu_fission" in m:
            if "nu_fission" in m:  # only the product is given
                nu_fission = _vector(m, "nu_fission", where, groups)
                fission = (
                    _vector(m, "fission", where, groups) if "fission" in m
                    else np.zeros(groups)
                )  # fmt: skip
            else:
                fission = _vector(m, "fission", where, groups)
                nu_fission = fission * _vector(m, "nu", where, groups)
            chi = _vector(m, "chi", where, groups)
            if abs(chi.sum() - 1) > 1e-3:
                raise LibraryError(f"{where}.chi: sums to {chi.sum():g}, not 1")
            chi = chi / chi.sum()
        else:
            fission = nu_fission = chi = np.zeros(groups)
        diffusion = _vector(m, "diffusion", where, groups) if "diffusion" in m else None
        if diffusion is not None and np.any(diffusion <= 0):
            raise LibraryError(f"{where}.diffusion: must be positive")
        out[name] = MacroXS(name, total, absorption, scatter, fission, nu_fission, chi,
                            diffusion)  # fmt: skip
    return out


def load_multigroup(path: Path) -> dict[str, MacroXS]:
    with open(path, "rb") as f:
        return parse_multigroup(tomllib.load(f), f"{path.name}: ")
