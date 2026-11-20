"""VERA problem 2: the 17x17 lattice model and its reference pin powers (C4)."""

import numpy as np
import pytest

from fission.model import load_model
from fission.pincell.assembly import reference_pin_powers


def test_vera_p2_lattice():
    m = load_model("vera-p2")
    lat = m.lattices["assembly"]
    flat = [p for row in lat.pins for p in row]
    assert flat.count("guide_tube") == 24
    assert flat.count("instrument_tube") == 1
    assert lat.pins[8][8] == "instrument_tube"
    assert lat.width == pytest.approx(21.50)  # Table 2: the assembly pitch
    # octant symmetry of the map
    grid = np.array(lat.pins)
    assert (grid == grid.T).all() and (grid == grid[::-1]).all()


@pytest.mark.parametrize(("state", "high", "low"),
                         [("2A", 1.0516, 0.9389), ("2D", 1.0474, 0.9404)])  # fmt: skip
def test_reference_pin_powers_unfold(state, high, low):
    ref = reference_pin_powers(state)
    m = load_model("vera-p2")
    fuel = np.array(m.lattices["assembly"].pins) == "fuel_rod"
    assert np.array_equal(~np.isnan(ref), fuel)  # values exactly at the rods
    assert np.nanmean(ref) == pytest.approx(1.0, abs=2e-4)  # normalized (4 digits)
    assert np.nanmax(ref) == high and np.nanmin(ref) == low  # Figures P2-3, P2-6
    np.testing.assert_array_equal(np.nan_to_num(ref), np.nan_to_num(ref.T))
