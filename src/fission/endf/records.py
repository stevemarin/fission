"""The ENDF-6 record layer: turning 80-column text lines into numbers.

Every line of an ENDF-6 file has the same fixed layout::

    columns  1-66   six data fields, 11 characters each
    columns 67-70   MAT   material number      (e.g. 9237 = U-238)
    columns 71-72   MF    file number          (kind of data: 1 header, 3 cross sections, ...)
    columns 73-75   MT    section number       (reaction: 1 total, 2 elastic, 102 capture, ...)
    columns 76-80   NS    line sequence number (optional, ignored here)

Lines are grouped into a handful of *record types*. Every section is a sequence
of these records, and the format manual describes each section in their terms:

    TEXT   66 characters of free text
    CONT   C1 C2 L1 L2 N1 N2   two floats then four integers
    HEAD   same layout as CONT, starting a section (C1 = ZA, C2 = AWR)
    LIST   a CONT whose N1 = NPL, followed by NPL floats, six per line
    TAB1   a CONT whose N1 = NR, N2 = NP, followed by NR interpolation
           pairs (NBT, INT) and NP (x, y) pairs: a tabulated function y(x)
    TAB2   a CONT with NR interpolation pairs for a function of two variables
           (the data for each value of the second variable follows as records)

Numbers are written in a compact Fortran style without the "E": 9.223800+4
means 9.2238e+4 and 1.0-5 means 1.0e-5.
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import NamedTuple

FIELD = 11  # characters per data field
FIELDS = 6  # data fields per line
DATA_COLS = 66  # columns holding data


def parse_float(field: str) -> float:
    """Parse an ENDF floating-point field: '9.223800+4' -> 92238.0, blank -> 0.0."""
    s = field.replace(" ", "")
    if not s:
        return 0.0
    try:
        return float(s)  # handles '1.0E-5', '1.0e-5', '0.5', ...
    except ValueError:
        pass
    # compact form: the exponent starts at the last sign that isn't the leading one
    i = max(s.rfind("+"), s.rfind("-"))
    if i > 0:
        return float(f"{s[:i]}e{s[i:]}")
    raise ValueError(f"not an ENDF float: {field!r}")


def parse_int(field: str) -> int:
    s = field.strip()
    return int(s) if s else 0


def split_fields(line: str) -> list[str]:
    """The six 11-character data fields of a line (short lines are padded)."""
    data = line[:DATA_COLS].ljust(DATA_COLS)
    return [data[i : i + FIELD] for i in range(0, DATA_COLS, FIELD)]


def control(line: str) -> tuple[int, int, int]:
    """(MAT, MF, MT) from columns 67-75."""
    return parse_int(line[66:70]), parse_int(line[70:72]), parse_int(line[72:75])


class Cont(NamedTuple):
    """A CONT (or HEAD) record: two floats and four integers."""

    c1: float
    c2: float
    l1: int
    l2: int
    n1: int
    n2: int


class Tab1Data(NamedTuple):
    """Raw content of a TAB1 record: interpolation regions and (x, y) points."""

    breaks: list[int]  # NBT: 1-based index of the last point of each region
    laws: list[int]  # INT: interpolation law used in each region
    x: list[float]
    y: list[float]


class Cursor:
    """Reads records one after another from the lines of a section."""

    def __init__(self, lines: Sequence[str]):
        self.lines = lines
        self.pos = 0

    def _line(self) -> str:
        if self.pos >= len(self.lines):
            raise EOFError("read past the end of the section")
        line = self.lines[self.pos]
        self.pos += 1
        return line

    def text(self) -> str:
        return self._line()[:DATA_COLS]

    def cont(self) -> Cont:
        f = split_fields(self._line())
        return Cont(
            parse_float(f[0]),
            parse_float(f[1]),
            parse_int(f[2]),
            parse_int(f[3]),
            parse_int(f[4]),
            parse_int(f[5]),
        )

    head = cont  # a HEAD record has the same layout as a CONT record

    def _numbers(self, n: int, convert=parse_float) -> list:
        """The next n numbers, read six per line."""
        out: list = []
        while len(out) < n:
            out.extend(convert(f) for f in split_fields(self._line()))
        return out[:n]

    def list_record(self) -> tuple[Cont, list[float]]:
        c = self.cont()
        return c, self._numbers(c.n1)

    def _interpolation(self, nr: int) -> tuple[list[int], list[int]]:
        pairs = self._numbers(2 * nr, parse_int)
        return pairs[0::2], pairs[1::2]

    def tab1(self) -> tuple[Cont, Tab1Data]:
        c = self.cont()
        breaks, laws = self._interpolation(c.n1)
        xy = self._numbers(2 * c.n2)
        return c, Tab1Data(breaks, laws, xy[0::2], xy[1::2])

    def tab2(self) -> tuple[Cont, list[int], list[int]]:
        c = self.cont()
        breaks, laws = self._interpolation(c.n1)
        return c, breaks, laws
