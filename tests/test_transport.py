"""Flat-source regions, cyclic tracks and ray tracing (step B1)."""

import math

import numpy as np
import pytest

from fission.model import parse
from fission.transport import (
    CellMesh,
    Geometry,
    cyclic_tracks,
    split_ring,
    trace,
    trace_track,
    track_areas,
)
from fission.transport.tracks import BOTTOM, EDGES, LEFT, RIGHT, TOP

PIN = CellMesh((0.4096, 0.418, 0.475), ("fuel", "gap", "clad", "water"), 8)


def _pin(sectors=8):
    return Geometry(1.26, 1, 1, (CellMesh(PIN.radii, PIN.materials, sectors),))


def _mixed():
    """A 3 x 2 lattice of three different cells (not square, not symmetric)."""
    small = CellMesh((0.3,), ("absorber", "water"), 4)
    empty = CellMesh((), ("water",), 1)
    cells = (PIN, small, empty, empty, PIN, small)
    return Geometry(1.26, 3, 2, cells)


# ---------------------------------------------------------------- regions
def test_equal_area_rings():
    radii = split_ring(0.1, 0.4096, 3)
    areas = np.diff([0.1**2, *(r * r for r in radii)])
    assert areas == pytest.approx(np.full(3, areas[0]), rel=1e-13)
    assert radii[-1] == pytest.approx(0.4096, rel=1e-15)


@pytest.mark.parametrize("sectors", [1, 2, 3, 4, 8])
def test_cell_areas_are_exact(sectors):
    mesh = CellMesh(PIN.radii, PIN.materials, sectors)
    areas = np.array(mesh.areas(1.26))
    assert areas.sum() == pytest.approx(1.26**2, rel=1e-14)
    assert np.all(areas > 0)
    # sample the cell on a square grid (shifted off the sector lines)
    n = 400
    x = (np.arange(n) + 0.37) / n * 1.26 - 0.63
    counts = np.zeros(mesh.regions)
    for xi in x:
        for yi in x:
            counts[mesh.region(xi, yi)] += 1
    sampled = counts / counts.sum() * 1.26**2
    big = areas > 0.02  # the thin gap needs a finer grid
    assert sampled[big] == pytest.approx(areas[big], rel=0.01)


def test_lattice_from_model_keeps_map_orientation():
    model = parse(
        {
            "model": {"format": 1, "name": "t"},
            "materials": {
                m: {"temperature": "moderator", "atom_densities": {"H1": 1.0}}
                for m in ("a", "b", "w")
            },
            "pins": {
                "A": {"radii": [0.4], "materials": ["a", "w"]},
                "B": {"radii": [0.3, 0.5], "materials": ["b", "w", "b"]},
            },
            "lattices": {
                "L": {"pitch": 1.26, "legend": {"A": "A", "B": "B"}, "map": "B A\nA A"}
            },
        }
    )
    g = Geometry.from_lattice(model, "L", rings={"a": 2}, sectors=4)
    assert (g.nx, g.ny) == (2, 2)
    # the map's top-left B is cell (ix 0, iy 1)
    assert g.cells[2].materials == ("b", "w", "b")
    assert g.cells[0].materials == ("a", "a", "w")
    assert g.cells[0].radii[1] == pytest.approx(0.4)
    assert g.regions == 3 * 12 + 3 * 4
    assert len(g.materials()) == g.regions
    assert g.areas().sum() == pytest.approx((2 * 1.26) ** 2, rel=1e-14)
    assert g.materials()[g.region(0.63, 1.89)] == "b"  # centre of B


# ---------------------------------------------------------------- tracks
@pytest.mark.parametrize(("width", "height"), [(1.26, 1.26), (3.78, 2.52)])
def test_cyclic_angles_and_spacing(width, height):
    tracks = cyclic_tracks(width, height, 32, 0.05)
    assert tracks.weights.sum() == pytest.approx(1.0, rel=1e-14)
    assert len(tracks.angles) == 16
    for a in tracks.angles[:8]:
        assert 0 < a.phi < math.pi / 2
        assert math.tan(a.phi) == pytest.approx(
            height * a.nx / (width * a.ny), rel=1e-13
        )
        assert a.spacing <= 0.05
        # the same spacing measured along y
        assert (height / a.ny) * math.cos(a.phi) == pytest.approx(a.spacing, rel=1e-12)
    assert len(tracks) == 2 * sum(a.nx + a.ny for a in tracks.angles[:8])
    # each angle's strips tile the box exactly
    for ia, a in enumerate(tracks.angles):
        mine = tracks.angle == ia
        assert (a.spacing * tracks.lengths[mine]).sum() == pytest.approx(
            width * height, rel=1e-12
        )


def _direction(tracks, t, forward):
    phi = tracks.angles[int(tracks.angle[t])].phi
    sign = 1.0 if forward else -1.0
    return sign * math.cos(phi), sign * math.sin(phi)


def _ends(tracks, t, forward):
    """(entry point, exit point) of a sweep."""
    p, q = (tracks.x0[t], tracks.y0[t]), (tracks.x1[t], tracks.y1[t])
    return (p, q) if forward else (q, p)


@pytest.mark.parametrize(("width", "height"), [(1.26, 1.26), (3.78, 2.52)])
def test_reflective_links_mirror_the_direction_at_the_same_point(width, height):
    tracks = cyclic_tracks(width, height, 16, 0.1)
    seen = set()
    for t in range(len(tracks)):
        for forward in (True, False):
            nxt = tracks.next_fwd[t] if forward else tracks.next_bwd[t]
            nf = tracks.fwd_forward[t] if forward else tracks.bwd_forward[t]
            seen.add((int(nxt), bool(nf)))
            out = _ends(tracks, t, forward)[1]
            into = _ends(tracks, int(nxt), bool(nf))[0]
            assert into == pytest.approx(out, abs=1e-12)
            (dx, dy), (ex, ey) = (_direction(tracks, t, forward),
                                  _direction(tracks, int(nxt), bool(nf)))  # fmt: skip
            on_side = abs(out[0]) < 1e-12 or abs(out[0] - width) < 1e-12
            mirrored = (-dx, dy) if on_side else (dx, -dy)
            assert (ex, ey) == pytest.approx(mirrored, abs=1e-12)
    # every sweep is entered from exactly one other: the links are a permutation
    assert len(seen) == 2 * len(tracks)
    lengths = [len(tracks.cycle(t)) for t in range(0, len(tracks), 7)]
    assert min(lengths) >= 2  # cycles close (cycle() raises otherwise)


def test_vacuum_sides_end_the_sweep():
    bc = {
        "left": "reflective",
        "bottom": "reflective",
        "right": "vacuum",
        "top": "vacuum",
    }
    tracks = cyclic_tracks(3.78, 2.52, 16, 0.1, bc)
    for t in range(len(tracks)):
        x, y = tracks.x1[t], tracks.y1[t]  # forward sweeps go up
        leaves_top = abs(y - 2.52) < 1e-12
        leaves_right = abs(x - 3.78) < 1e-12
        assert (tracks.next_fwd[t] < 0) == (leaves_top or leaves_right)
        x, y = tracks.x0[t], tracks.y0[t]
        assert (tracks.next_bwd[t] < 0) == (abs(x - 3.78) < 1e-12)
    path = tracks.cycle(0)
    assert len(path) >= 1  # ends at a vacuum side instead of closing
    assert set(EDGES) == {"left", "right", "bottom", "top"}
    assert (LEFT, RIGHT, BOTTOM, TOP) == (0, 1, 2, 3)


# ---------------------------------------------------------------- ray tracing
def test_chords_of_an_off_centre_track():
    y = 0.63 + 0.1
    segs = trace_track(_pin(), 0.0, y, 1.26, y)
    assert len(segs) == 10  # 7 ring pieces; sector lines cut the fuel chord in 4
    assert sum(s for _, s in segs) == pytest.approx(1.26, rel=1e-15)
    for ring, r in enumerate(PIN.radii):
        inside = sum(s for reg, s in segs if reg // 8 <= ring)
        assert inside == pytest.approx(2 * math.sqrt(r * r - 0.01), rel=1e-13)


def test_segments_cover_every_track():
    g = _mixed()
    tracks = cyclic_tracks(g.width, g.height, 16, 0.05)
    seg = trace(g, tracks, use_rust=False)
    per_track = np.add.reduceat(seg.length, seg.start[:-1])
    assert per_track == pytest.approx(tracks.lengths, rel=1e-13)
    assert seg.region.min() >= 0 and seg.region.max() < g.regions
    for t in range(0, len(tracks), 11):  # neighbours always differ
        regions, _ = seg.of(t)
        assert np.all(regions[1:] != regions[:-1])


@pytest.mark.parametrize("make", [_pin, _mixed])
def test_track_areas_converge_to_the_exact_areas(make):
    g = make()
    exact = g.areas()
    errors = []
    for spacing in (0.04, 0.01, 0.0025):
        tracks = cyclic_tracks(g.width, g.height, 32, spacing)
        areas = track_areas(g, tracks, trace(g, tracks))
        assert areas.sum(axis=1) == pytest.approx(g.width * g.height, rel=1e-12)
        errors.append(np.abs(tracks.weights @ areas - exact).max())
    assert errors[2] < errors[1] < errors[0]
    assert errors[2] < 2e-5  # cm^2; the cells are 1.59 cm^2


@pytest.mark.parametrize("make", [_pin, _mixed])
def test_rust_trace_matches_python(make):
    g = make()
    tracks = cyclic_tracks(g.width, g.height, 32, 0.02)
    py = trace(g, tracks, use_rust=False)
    rs = trace(g, tracks, use_rust=True)
    np.testing.assert_array_equal(rs.start, py.start)
    np.testing.assert_array_equal(rs.region, py.region)
    assert np.max(np.abs(rs.length - py.length)) <= 1e-12 * g.width


def test_outer_ring_only_refines_the_mesh():
    model = parse(
        {
            "model": {"format": 1, "name": "t"},
            "materials": {
                m: {"temperature": "moderator", "atom_densities": {"H1": 1.0}}
                for m in ("f", "w")
            },
            "pins": {
                "p": {"radii": [0.54], "materials": ["f", "w"]},
                "e": {"radii": [], "materials": ["w"]},
                "big": {"radii": [0.62], "materials": ["f", "w"]},
            },
            "lattices": {
                "L": {
                    "pitch": 1.26,
                    "legend": {"P": "p", "E": "e", "B": "big"},
                    "map": "P E\nB P",
                }
            },
        }
    )
    plain = Geometry.from_lattice(model, "L", sectors=4)
    ringed = Geometry.from_lattice(model, "L", sectors=4, outer_ring=0.6)
    pin, empty, big = (ringed.cells[i] for i in (2, 3, 0))  # map rows run from the top
    assert pin.radii == (0.54, 0.6) and pin.materials == ("f", "w", "w")
    assert empty.radii == (0.6,) and empty.materials == ("w", "w")
    assert big.radii == (0.62,)  # pin already reaches past the ring: unchanged

    # the same material areas, more regions
    def by_material(g):
        out = {}
        for m, a in zip(g.materials(), g.areas(), strict=True):
            out[m] = out.get(m, 0.0) + a
        return out

    assert by_material(ringed) == pytest.approx(by_material(plain), rel=1e-14)
    assert ringed.regions > plain.regions
