"""Energy group structures.

A multigroup calculation replaces sigma(E) by one value per energy group.
Groups are numbered from the *highest* energy down, the usual convention of
reactor physics (group 0 = fastest):

    group g covers [E_(g+1), E_g],   E_0 > E_1 > ... > E_G.

Standard structures are copied from NJOY2016 (src/groupr.f90, subroutine
gengpn), which carries the library structures in common use.
"""

from __future__ import annotations

import itertools
import math
from bisect import bisect_right
from dataclasses import dataclass

# XMAS 172-group structure (NJOY ign=18), boundaries in eV, ascending.
# Designed for light-water reactor lattice calculations (JEF project).
_XMAS_172 = (
    1.00001e-05, 0.003, 0.005, 0.0069, 0.01, 0.015,
    0.02, 0.025, 0.03, 0.035, 0.042, 0.05,
    0.058, 0.067, 0.077, 0.08, 0.095, 0.100001,
    0.115, 0.134, 0.14, 0.16, 0.18, 0.189,
    0.22, 0.248, 0.28, 0.3, 0.3145, 0.32,
    0.35, 0.391, 0.4, 0.433, 0.485, 0.5,
    0.54, 0.625, 0.705, 0.78, 0.79, 0.85,
    0.86, 0.91, 0.93, 0.95, 0.972, 0.986,
    0.996, 1.02, 1.035, 1.045, 1.071, 1.097,
    1.11, 1.12535, 1.15, 1.17, 1.235, 1.3,
    1.3375, 1.37, 1.44498, 1.475, 1.5, 1.59,
    1.67, 1.755, 1.84, 1.93, 2.02, 2.1,
    2.13, 2.36, 2.55, 2.6, 2.72, 2.76792,
    3.3, 3.38075, 4.0, 4.12925, 5.04348, 5.34643,
    6.16012, 7.52398, 8.31529, 9.18981, 9.90555, 11.2245,
    13.7096, 15.9283, 19.4548, 22.6033, 24.9805, 27.6077,
    30.5113, 33.7201, 37.2665, 40.169, 45.5174, 48.2516,
    51.578, 55.5951, 67.9041, 75.6736, 91.6609, 136.742,
    148.625, 203.995, 304.325, 371.703, 453.999, 677.287,
    748.518, 914.242, 1010.39, 1234.1, 1433.82, 1507.33,
    2034.68, 2248.67, 3354.63, 3526.62, 5004.51, 5530.84,
    7465.86, 9118.82, 11137.8, 15034.4, 16615.6, 24787.5,
    27394.4, 29283.0, 36978.6, 40867.7, 55165.6, 67379.5,
    82297.5, 111090.0, 122773.0, 183156.0, 247235.0, 273237.0,
    301974.0, 407622.0, 450492.0, 497871.0, 550232.0, 608101.0,
    820850.0, 907180.0, 1002590.0, 1108030.0, 1224560.0, 1353350.0,
    1652990.0, 2018970.0, 2231300.0, 2465970.0, 3011940.0, 3678790.0,
    4493290.0, 5488120.0, 6065310.0, 6703200.0, 8187310.0, 10000000.0,
    11618300.0, 13840300.0, 14918200.0, 17332500.0, 19640300.0,
)  # fmt: skip


@dataclass(frozen=True)
class GroupStructure:
    """Group boundaries in eV, descending: boundaries[g] is the top of group g."""

    name: str
    boundaries: tuple[float, ...]

    def __post_init__(self):
        b = self.boundaries
        if len(b) < 2 or any(x <= y for x, y in itertools.pairwise(b)):
            raise ValueError("group boundaries must be strictly decreasing")

    @classmethod
    def from_ascending(cls, name: str, edges) -> GroupStructure:
        return cls(name, tuple(sorted(edges, reverse=True)))

    @property
    def size(self) -> int:
        return len(self.boundaries) - 1

    def bounds(self, g: int) -> tuple[float, float]:
        """(E_low, E_high) of group g."""
        return self.boundaries[g + 1], self.boundaries[g]

    def lethargy_width(self, g: int) -> float:
        lo, hi = self.bounds(g)
        return math.log(hi / lo)

    def find(self, energy: float) -> int | None:
        """The group containing `energy`, or None outside the structure."""
        b = self.boundaries
        if not b[-1] <= energy <= b[0]:
            return None
        ascending = b[::-1]
        i = min(bisect_right(ascending, energy) - 1, self.size - 1)
        return self.size - 1 - i


XMAS_172 = GroupStructure.from_ascending("XMAS-172", _XMAS_172)


def equal_lethargy(
    n: int, lo: float, hi: float, name: str | None = None
) -> GroupStructure:
    """n groups of equal lethargy width between lo and hi (eV)."""
    edges = [hi * (lo / hi) ** (k / n) for k in range(n + 1)]
    edges[-1] = lo
    return GroupStructure(name or f"{n} equal-lethargy groups", tuple(edges))
