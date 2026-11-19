"""Cyclic tracks: the straight lines along which MOC integrates.

The method of characteristics follows neutrons along straight lines (the
characteristics) crossing the geometry in a set of directions. For each
azimuthal angle phi, parallel tracks a distance delta apart cover the
rectangle W x H. A track stands for the strip of width delta around it, so
sum over the tracks of delta * (length inside a region) estimates the
region's area, and fluxes integrated along the tracks give region averages.

**Reflective boundaries.** A neutron leaving through a reflective side comes
back with its direction mirrored. For the angular flux to be passed on
exactly, the mirrored track must start where the old one ends. This holds
if the tracks of phi and pi - phi hit every side at the same evenly spaced
points. With n_x tracks starting on the bottom side and n_y on the left,
they are spaced W/n_x along x and H/n_y along y, which fixes the angle:

    tan phi = (H / n_y) / (W / n_x)  ->  phi = atan(H n_x / (W n_y)),
    delta = (W / n_x) sin phi.

The requested angle and spacing are therefore adjusted: n_x and n_y are
the counts needed for the requested spacing, rounded up, and the angle and
spacing then follow from them. Following the reflections, every track then
returns to its start after a finite number of segments: the tracks are
**cyclic**.

Angles are kept in (0, pi): the first quarter phi_a and its mirror
pi - phi_a. Each track is swept both ways, so directions in (pi, 2 pi) are
the backward sweeps. A track's forward exit and backward exit each link to
(next track, forward?) or to nothing (vacuum).

Azimuthal weights: angle a stands for the arc between the midpoints of its
neighbours (0 and pi/2 at the ends of the quarter), as a fraction of pi.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

LEFT, RIGHT, BOTTOM, TOP = 0, 1, 2, 3
EDGES = ("left", "right", "bottom", "top")


@dataclass(frozen=True)
class Angle:
    phi: float  # radians, in (0, pi)
    spacing: float  # cm, perpendicular distance between tracks
    weight: float  # fraction of the half circle; the weights sum to 1
    nx: int  # tracks entering through the bottom side
    ny: int  # tracks entering through the left (phi < pi/2) or right side


@dataclass
class TrackSet:
    width: float
    height: float
    angles: list[Angle]
    boundary: dict[str, str]  # edge -> "reflective" | "vacuum"
    angle: np.ndarray  # angle index of every track
    x0: np.ndarray  # start point
    y0: np.ndarray
    x1: np.ndarray  # end point
    y1: np.ndarray
    # where the forward (and backward) sweep continues: track index or -1
    # (vacuum), and whether it is swept forward there
    next_fwd: np.ndarray = field(repr=False)
    fwd_forward: np.ndarray = field(repr=False)
    next_bwd: np.ndarray = field(repr=False)
    bwd_forward: np.ndarray = field(repr=False)

    def __len__(self) -> int:
        return len(self.x0)

    @property
    def weights(self) -> np.ndarray:
        """Azimuthal weight of every angle (they sum to 1)."""
        return np.array([a.weight for a in self.angles])

    @property
    def lengths(self) -> np.ndarray:
        return np.sqrt((self.x1 - self.x0) ** 2 + (self.y1 - self.y0) ** 2)

    def cycle(self, track: int, forward: bool = True) -> list[tuple[int, bool]]:
        """The sweeps visited from (track, forward) until it comes back, or
        until it leaves through a vacuum side."""
        out = [(track, forward)]
        t, f = track, forward
        while True:
            t, f = (
                (int(self.next_fwd[t]), bool(self.fwd_forward[t]))
                if f
                else (int(self.next_bwd[t]), bool(self.bwd_forward[t]))
            )
            if t < 0 or (t, f) == (track, forward):
                return out
            out.append((t, f))
            if len(out) > 4 * len(self):
                raise RuntimeError("track links do not close")


def _angles(width: float, height: float, n_azim: int, spacing: float) -> list[Angle]:
    if n_azim < 4 or n_azim % 4:
        raise ValueError(
            "n_azim (directions in [0, 2 pi)) must be a positive multiple of 4"
        )
    quarter = n_azim // 4
    first = []
    for a in range(quarter):
        phi = math.pi / 2 * (a + 0.5) / quarter  # requested angle, evenly spread
        nx = int(width / spacing * math.sin(phi)) + 1
        ny = int(height / spacing * math.cos(phi)) + 1
        phi = math.atan(height * nx / (width * ny))  # cyclic angle
        first.append((phi, (width / nx) * math.sin(phi), nx, ny))
    phis = [p for p, _, _, _ in first]
    out = []
    for a, (phi, delta, nx, ny) in enumerate(first):
        lo = 0.0 if a == 0 else 0.5 * (phis[a - 1] + phi)
        hi = math.pi / 2 if a == quarter - 1 else 0.5 * (phi + phis[a + 1])
        out.append(Angle(phi, delta, (hi - lo) / math.pi, nx, ny))
    # the mirror angles pi - phi, in the same order
    out += [Angle(math.pi - a.phi, a.spacing, a.weight, a.nx, a.ny) for a in out]
    return out


def _exit(
    x: float, y: float, c: float, s: float, w: float, h: float
) -> tuple[float, float, int]:
    """Where the ray from (x, y) in direction (c, s) leaves the box, and through which side."""
    tx = (w - x) / c if c > 0 else (-x / c if c < 0 else math.inf)
    ty = (h - y) / s if s > 0 else (-y / s if s < 0 else math.inf)
    if tx < ty:
        return x + tx * c, y + tx * s, RIGHT if c > 0 else LEFT
    return x + ty * c, y + ty * s, TOP if s > 0 else BOTTOM


def cyclic_tracks(
    width: float,
    height: float,
    n_azim: int,
    spacing: float,
    boundary: str | dict[str, str] = "reflective",
) -> TrackSet:
    """Cyclic tracks over a W x H box for n_azim directions in [0, 2 pi).

    boundary: one condition for all sides, or a dict per side ("left",
    "right", "bottom", "top"), each "reflective" or "vacuum".
    """
    bc = {e: boundary for e in EDGES} if isinstance(boundary, str) else dict(boundary)
    if set(bc) != set(EDGES) or not set(bc.values()) <= {"reflective", "vacuum"}:
        raise ValueError(f"boundary needs {EDGES}, each 'reflective' or 'vacuum'")
    angles = _angles(width, height, n_azim, spacing)

    # Track starts and ends, and a lookup (angle, start/end, side, k) -> track.
    # Along the bottom and top, points sit at (k + 1/2) W / nx; along the
    # sides at (k + 1/2) H / ny.
    rows: list[tuple[int, float, float, float, float]] = []
    where: dict[tuple[int, str, int, int], int] = {}

    def slot(side: int, x: float, y: float, a: Angle) -> int:
        pos, step = (x, width / a.nx) if side in (BOTTOM, TOP) else (y, height / a.ny)
        k = round(pos / step - 0.5)
        if abs(pos - (k + 0.5) * step) > 1e-9 * max(width, height):
            raise RuntimeError("track end is not on the cyclic grid")
        return k

    for ia, a in enumerate(angles):
        c, s = math.cos(a.phi), math.sin(a.phi)
        side_in = LEFT if c > 0 else RIGHT
        starts = [((k + 0.5) * width / a.nx, 0.0, BOTTOM) for k in range(a.nx)]
        x_side = 0.0 if side_in == LEFT else width
        starts += [(x_side, (k + 0.5) * height / a.ny, side_in) for k in range(a.ny)]
        for sx, sy, side0 in starts:
            ex, ey, side1 = _exit(sx, sy, c, s, width, height)
            t = len(rows)
            rows.append((ia, sx, sy, ex, ey))
            where[(ia, "start", side0, slot(side0, sx, sy, a))] = t
            where[(ia, "end", side1, slot(side1, ex, ey, a))] = t

    arr = np.array(rows)
    angle = arr[:, 0].astype(np.int64)
    x0, y0, x1, y1 = (np.ascontiguousarray(arr[:, i]) for i in range(1, 5))
    n = len(rows)
    links = {name: np.full(n, -1, dtype=np.int64) for name in ("nf", "nb")}
    forward = {name: np.ones(n, dtype=bool) for name in ("ff", "bf")}
    half = len(angles) // 2
    mirror = [(i + half) % len(angles) for i in range(len(angles))]  # phi <-> pi - phi

    def link(t: int, x: float, y: float, dx: float, dy: float) -> tuple[int, bool]:
        """The sweep that continues one leaving (x, y) in direction (dx, dy)."""
        side = _side(x, y, width, height)
        if bc[EDGES[side]] == "vacuum":
            return -1, True
        if side in (LEFT, RIGHT):
            dx = -dx
        else:
            dy = -dy
        m = mirror[int(angle[t])]  # the mirrored direction lies along angle m
        a = angles[m]
        k = slot(side, x, y, a)
        if dy > 0:  # forward sweeps of angles in (0, pi) go up: start here
            return where[(m, "start", side, k)], True
        return where[(m, "end", side, k)], False  # backward sweep from its end

    for t in range(n):
        phi = angles[int(angle[t])].phi
        c, s = math.cos(phi), math.sin(phi)
        links["nf"][t], forward["ff"][t] = link(t, x1[t], y1[t], c, s)
        links["nb"][t], forward["bf"][t] = link(t, x0[t], y0[t], -c, -s)
    return TrackSet(
        width, height, angles, bc, angle, x0, y0, x1, y1,
        links["nf"], forward["ff"], links["nb"], forward["bf"],
    )  # fmt: skip


def _side(x: float, y: float, width: float, height: float) -> int:
    """The side of the box a boundary point lies on."""
    d = [abs(x), abs(x - width), abs(y), abs(y - height)]
    return int(np.argmin(d))
