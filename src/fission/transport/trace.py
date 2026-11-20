"""Ray tracing: cut every track into segments, one per region crossed.

Along a track p(t) = p0 + t d (|d| = 1, 0 <= t <= L) the region changes
only where the track crosses a surface. There are three kinds of surface:

* the grid lines between cells, x = g + i p and y = g + j p (g: the gap
  around the lattice, usually 0): t = (g + i p - x0) / d_x, and the same for y;
* the circles of a cell, |p(t) - c| = r:
  t^2 + 2 b t + q = 0 with b = d . (p0 - c), q = |p0 - c|^2 - r^2,
  so t = -b +- sqrt(b^2 - q) when b^2 > q;
* the sector lines, through the centre c in direction u:
  t = ((c - p0) x u) / (d x u), with a x b = a_x b_y - a_y b_x.

Collecting all crossings, sorting them, and asking which region contains the
midpoint of each piece gives the segments. That is slower than walking from
surface to surface, but it is simple and cannot get lost. Pieces shorter
than 1e-10 cm (tangent crossings, lines through a corner) are dropped, and
neighbouring pieces in the same region are merged.

With every segment known, sum_tracks delta * (length in region i) estimates
the area of region i. This is the first check of a tracking: it must
approach the exact area as the spacing delta shrinks.
"""

from __future__ import annotations

import itertools
import math
from dataclasses import dataclass

import numpy as np

from .geometry import TWO_PI, Geometry
from .tracks import TrackSet

MIN_LENGTH = 1e-10  # cm


@dataclass
class Segments:
    """All segments, track after track: those of track i are start[i]:start[i+1]."""

    region: np.ndarray  # int64
    length: np.ndarray  # cm
    start: np.ndarray  # int64, len(tracks) + 1

    def of(self, track: int) -> tuple[np.ndarray, np.ndarray]:
        a, b = self.start[track], self.start[track + 1]
        return self.region[a:b], self.length[a:b]


def trace_track(
    geometry: Geometry, x0: float, y0: float, x1: float, y1: float
) -> list[tuple[int, float]]:
    """(region, length) of every segment of the track from (x0, y0) to (x1, y1)."""
    p, g = geometry.pitch, geometry.gap
    length = math.sqrt((x1 - x0) * (x1 - x0) + (y1 - y0) * (y1 - y0))
    dx, dy = (x1 - x0) / length, (y1 - y0) / length
    # crossings of the grid lines between cells
    cuts = [0.0, length]
    if dx != 0.0:
        for i in range(1, geometry.nx):
            t = (g + i * p - x0) / dx
            if 0.0 < t < length:
                cuts.append(t)
    if dy != 0.0:
        for j in range(1, geometry.ny):
            t = (g + j * p - y0) / dy
            if 0.0 < t < length:
                cuts.append(t)
    cuts.sort()
    offsets = geometry.offsets
    out: list[tuple[int, float]] = []
    for a, b in itertools.pairwise(cuts):
        if b - a < MIN_LENGTH:
            continue
        tm = 0.5 * (a + b)
        ix, iy = geometry.cell_at(x0 + tm * dx, y0 + tm * dy)
        cell = iy * geometry.nx + ix
        mesh = geometry.cells[cell]
        cx, cy = geometry.centre(ix, iy)
        # crossings inside this cell: circles, then sector lines
        local = [a, b]
        ox, oy = x0 - cx, y0 - cy
        bb = dx * ox + dy * oy
        for r in mesh.radii:
            disc = bb * bb - (ox * ox + oy * oy - r * r)
            if disc > 0.0:
                root = math.sqrt(disc)
                for t in (-bb - root, -bb + root):
                    if a < t < b:
                        local.append(t)
        if mesh.sectors > 1:
            for k in range(mesh.sectors):
                theta = TWO_PI * k / mesh.sectors
                ux, uy = math.cos(theta), math.sin(theta)
                den = dx * uy - dy * ux
                if den != 0.0:
                    t = ((cx - x0) * uy - (cy - y0) * ux) / den
                    if a < t < b:
                        local.append(t)
        local.sort()
        for s, e in itertools.pairwise(local):
            if e - s < MIN_LENGTH:
                continue
            tm = 0.5 * (s + e)
            region = offsets[cell] + mesh.region(x0 + tm * dx - cx, y0 + tm * dy - cy)
            if out and out[-1][0] == region:
                out[-1] = (region, out[-1][1] + (e - s))
            else:
                out.append((region, e - s))
    return out


def _cell_arrays(geometry: Geometry) -> tuple[list[list[float]], list[int], list[int]]:
    return (
        [list(c.radii) for c in geometry.cells],
        [c.sectors for c in geometry.cells],
        geometry.offsets,
    )


def trace(geometry: Geometry, tracks: TrackSet, use_rust: bool = True) -> Segments:
    """Segments of every track (threaded in Rust; use_rust=False is the reference)."""
    if (tracks.width, tracks.height) != (geometry.width, geometry.height):
        raise ValueError("tracks and geometry cover different rectangles")
    if use_rust:
        from fission import _rust

        region, length, start = _rust.trace_tracks(
            geometry.pitch, geometry.nx, geometry.ny, *_cell_arrays(geometry),
            tracks.x0, tracks.y0, tracks.x1, tracks.y1, geometry.gap,
        )  # fmt: skip
        return Segments(region, length, start)
    regions: list[int] = []
    lengths: list[float] = []
    starts = [0]
    for t in range(len(tracks)):
        segs = trace_track(
            geometry, tracks.x0[t], tracks.y0[t], tracks.x1[t], tracks.y1[t]
        )
        regions.extend(r for r, _ in segs)
        lengths.extend(s for _, s in segs)
        starts.append(len(regions))
    return Segments(
        np.array(regions, dtype=np.int64),
        np.array(lengths),
        np.array(starts, dtype=np.int64),
    )


def track_areas(geometry: Geometry, tracks: TrackSet, segments: Segments) -> np.ndarray:
    """Area of every region estimated by each angle: shape (angles, regions).

    Row a is sum over the tracks of angle a of spacing * segment length. The
    weighted average over angles is ``tracks.weights @ track_areas(...)``.
    """
    out = np.zeros((len(tracks.angles), geometry.regions))
    track_of_segment = np.repeat(np.arange(len(tracks)), np.diff(segments.start))
    angle = tracks.angle[track_of_segment]
    spacing = np.array([a.spacing for a in tracks.angles])[angle]
    np.add.at(out, (angle, segments.region), spacing * segments.length)
    return out


# ---------------------------------------------------------------- cell faces
@dataclass
class Crossings:
    """Where tracks cross the faces of the lattice cells (the CMFD mesh).

    Faces are numbered x-faces first (the vertical line x = i p in row iy is
    iy * (nx + 1) + i), then y-faces (y = j p in column ix is
    n_x_faces + j * nx + ix); the outer boundary is included. For track t,
    crossings start[t]:start[t+1] each give the position along the track
    (0 = its start, n = after its last segment; position k is just before
    segment k), the face, and the sign: +1 if the forward sweep crosses in
    the +x (or +y) direction.
    """

    start: np.ndarray  # int64, tracks + 1
    position: np.ndarray  # int64
    face: np.ndarray  # int64
    sign: np.ndarray  # float64, +1 or -1
    n_faces: int


def face_counts(geometry: Geometry) -> tuple[int, int]:
    """(number of x-faces, number of y-faces)."""
    return geometry.ny * (geometry.nx + 1), (geometry.ny + 1) * geometry.nx


def crossings(geometry: Geometry, tracks: TrackSet, segments: Segments) -> Crossings:
    """Face crossings of every track, from the cells of consecutive segments."""
    nx, ny = geometry.nx, geometry.ny
    n_xf, n_yf = face_counts(geometry)
    tol = 1e-9 * max(geometry.width, geometry.height)
    cell = geometry.cell_of_regions()[segments.region]
    rows: list[tuple[int, int, int, float]] = []  # (track, position, face, sign)

    def x_face(i: int, iy: int) -> int:
        return iy * (nx + 1) + i

    def y_face(ix: int, j: int) -> int:
        return n_xf + j * nx + ix

    for t in range(len(tracks)):
        a, b = int(segments.start[t]), int(segments.start[t + 1])
        cells = cell[a:b]
        # entering through the outer boundary (forward sweep)
        ix, iy = int(cells[0] % nx), int(cells[0] // nx)
        x0, y0 = tracks.x0[t], tracks.y0[t]
        if abs(x0) < tol:
            rows.append((t, 0, x_face(0, iy), 1.0))
        elif abs(x0 - geometry.width) < tol:
            rows.append((t, 0, x_face(nx, iy), -1.0))
        elif abs(y0) < tol:
            rows.append((t, 0, y_face(ix, 0), 1.0))
        else:
            rows.append((t, 0, y_face(ix, ny), -1.0))
        # between cells; a diagonal step (through a corner) counts as an
        # x-crossing followed by a y-crossing
        for k in range(1, b - a):
            c1, c2 = int(cells[k - 1]), int(cells[k])
            if c1 == c2:
                continue
            ix1, iy1, ix2, iy2 = c1 % nx, c1 // nx, c2 % nx, c2 // nx
            dx, dy = ix2 - ix1, iy2 - iy1
            if abs(dx) > 1 or abs(dy) > 1:
                raise RuntimeError(f"track {t} jumps from cell {c1} to {c2}")
            if dx:
                rows.append((t, k, x_face(max(ix1, ix2), iy1), float(dx)))
            if dy:
                rows.append((t, k, y_face(ix2, max(iy1, iy2)), float(dy)))
        # leaving through the outer boundary
        ix, iy = int(cells[-1] % nx), int(cells[-1] // nx)
        x1, y1 = tracks.x1[t], tracks.y1[t]
        n = b - a
        if abs(x1 - geometry.width) < tol:
            rows.append((t, n, x_face(nx, iy), 1.0))
        elif abs(x1) < tol:
            rows.append((t, n, x_face(0, iy), -1.0))
        elif abs(y1 - geometry.height) < tol:
            rows.append((t, n, y_face(ix, ny), 1.0))
        else:
            rows.append((t, n, y_face(ix, 0), -1.0))
    arr = np.array(rows)
    counts = np.bincount(arr[:, 0].astype(np.int64), minlength=len(tracks))
    return Crossings(
        np.concatenate([[0], np.cumsum(counts)]).astype(np.int64),
        np.ascontiguousarray(arr[:, 1].astype(np.int64)),
        np.ascontiguousarray(arr[:, 2].astype(np.int64)),
        np.ascontiguousarray(arr[:, 3]),
        n_xf + n_yf,
    )
