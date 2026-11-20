"""Cross-section tables (M3): multilinear interpolation and fitted polynomials."""

import itertools

import numpy as np
import pytest

from fission.model import BRANCH_AXES
from fission.tables import (
    XSTable,
    coordinate,
    fit,
    k_infinity,
    multilinear,
    press,
    select,
)

AXES = {
    "fuel_temperature": (565.0, 900.0, 1200.0),
    "moderator_temperature": (565.0, 600.0),
    "moderator_density": (0.66, 0.74),
    "boron_ppm": (0.0, 1000.0, 2000.0),
}


def _table(f, transforms=None) -> XSTable:
    """A table whose 'total' is f(coordinates) and whose other constants
    satisfy total = sigma_a + sum of sigma_s rows, at every node."""
    transforms = {} if transforms is None else transforms
    shape = tuple(len(AXES[a]) for a in BRANCH_AXES)
    total = np.zeros(shape + (2,))
    sigma_s = np.zeros(shape + (2, 2))
    for index in itertools.product(*(range(n) for n in shape)):
        x = [coordinate(a, AXES[a][i], transforms) for a, i in zip(BRANCH_AXES, index)]
        total[index] = [f(x), 2 * f(x) + 1]
        sigma_s[index] = [[0.5 * f(x), 0.1], [0.0, f(x)]]
    data = {
        "total": total,
        "sigma_s": sigma_s,
        "sigma_a": total - sigma_s.sum(axis=-1),
        "nu_sigma_f": np.broadcast_to([0.01, 0.2], shape + (2,)).copy(),
        "chi": np.broadcast_to([1.0, 0.0], shape + (2,)).copy(),
        "D": np.ones(shape + (2,)),
    }
    return XSTable(dict(AXES), data, transforms)


def multi(x):  # linear in each variable, cross terms included
    return (
        1
        + 1e-3 * x[0]
        + 2 * x[2]
        + 1e-4 * x[3]
        + 1e-6 * x[0] * x[3]
        + 1e-3 * x[1] * x[2]
    )


def test_multilinear_is_exact_for_multilinear_functions():
    t = _table(multi)
    for state in ({"fuel_temperature": 700.0, "moderator_temperature": 580.0,
                   "moderator_density": 0.70, "boron_ppm": 1500.0},
                  {"fuel_temperature": 1300.0, "moderator_temperature": 610.0,
                   "moderator_density": 0.65, "boron_ppm": -100.0}):  # fmt: skip
        x = [state[a] for a in BRANCH_AXES]
        assert t("total", state)[0] == pytest.approx(multi(x), rel=1e-13)
    # nodes are reproduced, and a single-node axis is a constant
    nodes = [np.array([1.0, 2.0, 4.0]), np.array([7.0])]
    values = np.array([[1.0], [3.0], [2.0]])
    assert multilinear(nodes, values, [2.0, 99.0]) == pytest.approx(3.0)
    assert multilinear(nodes, values, [3.0, 0.0]) == pytest.approx(2.5)
    assert multilinear(nodes, values, [5.0, 0.0]) == pytest.approx(1.5)  # extrapolated


def test_interpolation_keeps_linear_relations_between_constants():
    t = _table(lambda x: 1 + np.sin(x[0] / 300) * x[3] / 2000)  # not multilinear
    state = {"fuel_temperature": 777.0, "moderator_temperature": 590.0,
             "moderator_density": 0.71, "boron_ppm": 1234.0}  # fmt: skip
    c = t.constants(state)
    np.testing.assert_allclose(c["total"], c["sigma_a"] + c["sigma_s"].sum(axis=1),
                               rtol=1e-14)  # fmt: skip
    assert k_infinity(c) > 0


def test_sqrt_axis_is_exact_for_sqrt_t_dependence():
    def f(x):
        return 2 + 0.01 * x[0]  # x[0] = sqrt(T) when transformed

    t = _table(f, transforms={"fuel_temperature": "sqrt"})
    state = {"fuel_temperature": 700.0, "moderator_temperature": 565.0,
             "moderator_density": 0.66, "boron_ppm": 0.0}  # fmt: skip
    assert t("total", state)[0] == pytest.approx(2 + 0.01 * np.sqrt(700.0), rel=1e-13)


def test_press_residuals_are_leave_one_out_refits():
    rng = np.random.default_rng(2)
    x = np.column_stack([np.ones(12), rng.normal(size=12), rng.normal(size=12)])
    y = rng.normal(size=(12, 2))
    e = press(x, y)
    for i in range(12):
        keep = np.arange(12) != i
        coef, *_ = np.linalg.lstsq(x[keep], y[keep], rcond=None)
        np.testing.assert_allclose(e[i], y[i] - x[i] @ coef, rtol=1e-10)


def test_selection_finds_the_terms_of_a_polynomial():
    def poly(x):  # in the normalized variables of the grid
        u0 = (x[0] - 882.5) / 317.5
        u3 = (x[3] - 1000.0) / 1000.0
        return 3 + 0.4 * u0 + 0.2 * u3 * u0 + 0.05 * u3**2

    t = _table(poly)
    f, path = select(t, keys=("total",))
    assert set(f.terms) == {(0, 0, 0, 0), (1, 0, 0, 0), (1, 0, 0, 1), (0, 0, 0, 2)}
    state = {"fuel_temperature": 650.0, "moderator_temperature": 570.0,
             "moderator_density": 0.69, "boron_ppm": 300.0}  # fmt: skip
    x = [state[a] for a in BRANCH_AXES]
    assert f("total", state)[0] == pytest.approx(poly(x), rel=1e-12)
    assert path[-1][1] <= path[0][1]
    g = fit(t, f.terms)
    np.testing.assert_allclose(g("total", state), f("total", state), rtol=1e-14)


def test_table_file_round_trip(tmp_path):
    t = _table(multi, transforms={"fuel_temperature": "sqrt"})
    t.write(tmp_path / "t.json")
    back = XSTable.read(tmp_path / "t.json")
    state = {"fuel_temperature": 1000.0, "moderator_temperature": 575.0,
             "moderator_density": 0.7, "boron_ppm": 500.0}  # fmt: skip
    np.testing.assert_array_equal(back("sigma_s", state), t("sigma_s", state))
    one = t.without("fuel_temperature", 1)
    assert one.shape == (2, 2, 2, 3) and one.axes["fuel_temperature"] == (565.0, 1200.0)
