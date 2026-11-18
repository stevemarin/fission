"""Tabulated functions and the five ENDF interpolation laws.

ENDF stores a function y(x) as points (x_i, y_i) plus a rule for what happens
between points. The rule can change along the table: the table is split into
*regions*, region r ends at point NBT[r] (1-based) and uses law INT[r]:

    INT = 1  histogram    y = y_i                            (constant)
    INT = 2  lin-lin      y linear in x
    INT = 3  lin-log      y linear in ln(x)
    INT = 4  log-lin      ln(y) linear in x
    INT = 5  log-log      ln(y) linear in ln(x)              (power law)

Each law is exact for one family of functions, which is why evaluators pick
them: lin-lin is exact for straight lines, log-log for power laws such as the
1/v capture cross section (sigma ~ E^-1/2), log-lin for exponentials.

The evaluated data is defined *by* its interpolation law: interpolating an
ENDF table any other way gives a slightly different cross section.
"""

from __future__ import annotations

import math
from bisect import bisect_left, bisect_right
from collections.abc import Sequence
from itertools import pairwise

HISTOGRAM, LIN_LIN, LIN_LOG, LOG_LIN, LOG_LOG = 1, 2, 3, 4, 5


def interpolate(law: int, x: float, x1: float, x2: float, y1: float, y2: float) -> float:
    """Value at x between (x1, y1) and (x2, y2) using one ENDF interpolation law."""
    if law == HISTOGRAM or x2 == x1:
        return y1
    # logarithms need positive values; fall back to linear where they aren't
    log_x = law in (LIN_LOG, LOG_LOG) and x1 > 0 and x > 0
    log_y = law in (LOG_LIN, LOG_LOG) and y1 > 0 and y2 > 0
    if law not in (LIN_LIN, LIN_LOG, LOG_LIN, LOG_LOG):
        raise ValueError(f"unsupported interpolation law {law}")
    t = math.log(x / x1) / math.log(x2 / x1) if log_x else (x - x1) / (x2 - x1)
    if log_y:
        return y1 * math.exp(t * math.log(y2 / y1))
    return y1 + t * (y2 - y1)


class Tabulated1D:
    """A function y(x) given as ENDF points and interpolation regions.

    Outside [x_0, x_last] the function is 0, the ENDF convention for cross
    sections (a threshold reaction's table starts at its threshold).
    A repeated x value marks a discontinuity; exactly at it the value from the
    upper side is returned.
    """

    def __init__(self, x: Sequence[float], y: Sequence[float],
                 breaks: Sequence[int] | None = None, laws: Sequence[int] | None = None):
        if len(x) != len(y) or not x:
            raise ValueError("x and y must be non-empty and the same length")
        if any(b < a for a, b in pairwise(x)):
            raise ValueError("x must be non-decreasing")
        self.x = list(x)
        self.y = list(y)
        self.breaks = list(breaks) if breaks else [len(x)]
        self.laws = list(laws) if laws else [LIN_LIN]
        if len(self.breaks) != len(self.laws) or self.breaks[-1] != len(x):
            raise ValueError("interpolation regions must cover the whole table")

    def law(self, i: int) -> int:
        """Interpolation law of the interval between points i and i+1 (0-based)."""
        # the interval's upper point is i+2 in ENDF's 1-based numbering
        return self.laws[bisect_left(self.breaks, i + 2)]

    def __call__(self, x: float) -> float:
        xs = self.x
        if x < xs[0] or x > xs[-1]:
            return 0.0
        if x == xs[-1]:
            return self.y[-1]
        i = bisect_right(xs, x) - 1
        return interpolate(self.law(i), x, xs[i], xs[i + 1], self.y[i], self.y[i + 1])

    def __len__(self) -> int:
        return len(self.x)

    def __repr__(self) -> str:
        return (f"Tabulated1D({len(self.x)} points, x in [{self.x[0]:.4g}, {self.x[-1]:.4g}], "
                f"laws {self.laws})")
