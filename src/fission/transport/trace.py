"""Ray tracing: cut every track into segments, one per region crossed.

Along a track p(t) = p0 + t d (|d| = 1, 0 <= t <= L) the region changes
only where the track crosses a surface. There are three kinds of surface:

* the grid lines between cells, x = i p and y = j p:
  t = (i p - x0) / d_x, and the same for y;
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
    p = geometry.pitch
    length = math.sqrt((x1 - x0) * (x1 - x0) + (y1 - y0) * (y1 - y0))
    dx, dy = (x1 - x0) / length, (y1 - y0) / length
    # crossings of the grid lines between cells
    cuts = [0.0, length]
    if dx != 0.0:
        for i in range(1, geometry.nx):
            t = (i * p - x0) / dx
            if 0.0 < t < length:
                cuts.append(t)
    if dy != 0.0:
        for j in range(1, geometry.ny):
            t = (j * p - y0) / dy
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
        cx, cy = (ix + 0.5) * p, (iy + 0.5) * p
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
            tracks.x0, tracks.y0, tracks.x1, tracks.y1,
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
