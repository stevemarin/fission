"""Few-group constants for the nodal code (step C3): homogenize, B1, condense."""

from .b1 import Mode, critical, mode
from .condense import FewGroup, condense
from .homogenize import FineConstants, homogenize

__all__ = [
    "FewGroup",
    "FineConstants",
    "Mode",
    "condense",
    "critical",
    "homogenize",
    "mode",
]
