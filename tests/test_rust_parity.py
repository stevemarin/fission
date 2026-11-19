"""Rust kernels (fission._rust) must reproduce the Python references.

Tolerances are set by floating-point order of operations: the ports follow
the Python code step by step, so agreement is at the 1e-12 level or better.
"""

import math
import random

import numpy as np
import pytest

from fission import _rust
from fission.data import EndfLibrary
from fission.endf import Evaluation, Tabulated1D, fast
from fission.endf.channel import K_CONSTANT
from fission.endf.mlbw import mlbw
from fission.endf.reconstruct import reconstruct
from fission.endf.reich_moore import reich_moore
from fission.endf.resonances import LGroup, ResolvedRange, Resonance

REL = 1e-12


def test_wave_number_constant_identical():
    assert _rust.k_constant() == K_CONSTANT


@pytest.mark.parametrize("law", [1, 2, 3, 4, 5])
def test_interpolation_laws(law):
    table = Tabulated1D(
        [0.5, 1.0, 2.0, 2.0, 4.0, 9.0],
        [3.0, 2.0, 0.7, 1.4, 5.0, 2.5],
        breaks=[3, 6],
        laws=[law, 5],
    )
    rng = random.Random(law)
    xs = np.array(
        [0.4, 0.5, 2.0, 9.0, 9.5] + [rng.uniform(0.5, 9.0) for _ in range(500)]
    )
    rust = _rust.tab1_eval(table, xs)
    for x, r in zip(xs, rust, strict=True):
        assert r == pytest.approx(table(float(x)), rel=1e-13, abs=0)


def _compare(py, rs, scale):
    for a, b in zip((py.elastic, py.capture, py.fission), rs, strict=True):
        assert b == pytest.approx(a, rel=REL, abs=REL * scale)


def test_single_level_with_fission():
    res = Resonance(1.14, 4.0, 1.5e-5, 35e-3, 0.12, -0.05)
    group = LGroup(l=0, awri=233.0, apl=0.96, resonances=(res,))
    rr = ResolvedRange(1e-5, 2250.0, 3, 1, 3.5, 0.96, 0, 4, (group,))
    rm = fast.reich_moore(rr)
    for e in (1e-5, 0.0253, 1.0, 1.14, 1.2, 100.0):
        py = reich_moore(rr, e)
        _compare(py, rm(e), py.elastic + py.capture + py.fission)


@pytest.mark.parametrize(
    "res, spi, l",
    [
        (Resonance(292.0, 2.0, 0.6, 0.14, 0.0, 0.0), 2.5, 0),
        (Resonance(4008.0, 1.5, 0.34, 0.09, 0.0, 0.0), 0.0, 1),  # p-wave: shift factor
        (Resonance(-50.0, 0.5, 2.0, 0.03, 0.01, 0.0), 0.0, 0),  # bound level, fission
    ],
)
def test_mlbw_single_level(res, spi, l):
    group = LGroup(l=l, awri=90.0, apl=0.0, resonances=(res,))
    rr = ResolvedRange(1e-5, 2e5, 2, 0, spi, 0.72, 0, 0, (group,))
    model = fast.mlbw(rr)
    for e in (1e-5, 0.0253, 290.0, 292.0, 4000.0, 4008.0, 1e5):
        py = mlbw(rr, e)
        _compare(py, model(e), py.elastic + py.capture + py.fission)


# --------------------------------------------------- real ENDF/B-VIII.1 data
def _evaluation(nuclide):
    lib = EndfLibrary()
    if not (lib.root / "archives" / "neutrons-version.VIII.1.tar.gz").exists():
        pytest.skip("ENDF/B-VIII.1 neutron data not downloaded")
    return Evaluation.from_file(lib.neutron(nuclide))


@pytest.fixture(scope="module")
def u238():
    return _evaluation("U238")


@pytest.fixture(scope="module")
def u235():
    return _evaluation("U235")


@pytest.fixture(scope="module")
def zr91():
    return _evaluation("Zr91")


@pytest.fixture(scope="module")
def u234():
    return _evaluation("U234")


@pytest.mark.parametrize("name", ["zr91", "u234"])
def test_mlbw_point_values(name, request):
    rr = request.getfixturevalue(name).resonances[0].resolved
    model = fast.mlbw(rr)
    for i in range(150):
        e = 1e-5 * (rr.eh / 1e-5) ** (i / 149)
        py = mlbw(rr, e)
        _compare(py, model(e), py.elastic + py.capture + py.fission)


@pytest.mark.parametrize("name", ["u238", "u235"])
def test_point_values(name, request):
    ev = request.getfixturevalue(name)
    rr = ev.resonances[0].resolved
    rm = fast.reich_moore(rr)
    energies = [1e-5 * (rr.eh / 1e-5) ** (i / 120) for i in range(120)]
    energies += [6.55, 6.67428, 20.87, 0.29, 1.14]  # dips and peaks
    for e in energies:
        py = reich_moore(rr, e)
        _compare(py, rm(e), py.elastic + py.capture + py.fission)
    # vectorized entry point gives the same numbers as the scalar one
    el, ca, fi = rm.evaluate(np.array(energies))
    assert [rm(e) for e in energies] == list(zip(el, ca, fi, strict=True))


@pytest.mark.parametrize(
    "name, lo, hi",
    [
        ("u238", 6.3, 7.0),
        ("u235", 0.9, 1.4),
        ("zr91", 200.0, 2000.0),
        ("u234", 1.0, 30.0),
    ],
)
def test_reconstruction_window(name, lo, hi, request):
    """Same seeds and same refinement decisions: the grids should be identical."""
    ev = request.getfixturevalue(name)
    py, rs = reconstruct(ev, lo, hi), fast.reconstruct(ev, lo, hi)
    assert rs.energies == py.energies
    # errors are judged against the total at each point: inside interference
    # dips elastic is a small difference of ~11 b terms
    totals = [
        e + c + f for e, c, f in zip(py.elastic, py.capture, py.fission, strict=True)
    ]
    for reaction in ("elastic", "capture", "fission"):
        pairs = zip(py.reaction(reaction), rs.reaction(reaction), totals, strict=True)
        for a, b, total in pairs:
            assert b == pytest.approx(a, rel=REL, abs=REL * total)


@pytest.mark.slow
def test_u238_resonance_integral(u238):
    """Capture resonance integral of the resolved range (0.5 eV - 20 keV).

    Published total RI (0.5 eV to infinity) is ~275 b (Mughabghab: 275 +- 3 b);
    the part above 20 keV adds only a few barns, so the resolved part must sit
    just below it.
    """
    p = fast.reconstruct(u238)
    ri = p.integral_over_e("capture", 0.5, 2e4)
    assert 268.0 < ri < 277.0
    assert math.isfinite(p.integral_over_e("elastic", 0.5, 2e4))


# ---------------------------------------------------------------- slowing-down march
def _march_problem(regions):
    from fission.pincell.shielding import march  # noqa: F401  (the reference)

    e = np.geomspace(2e3 * math.exp(0.8), 1.0, 20000)
    u = np.log(e[0] / e)
    start = int(np.searchsorted(u, 0.8)) + 1
    res = 50.0 / (1 + ((e - 100.0) / 0.25) ** 2) + 30.0 / (1 + ((e - 36.0) / 0.1) ** 2)
    moderator = np.full_like(e, 1.0)
    absorber = 0.1 + 0.3 * res / 50.0  # some resonance scattering
    total = moderator + absorber + res
    scatter = [moderator, absorber]
    alphas = [0.5, (237.0 / 239.0) ** 2]
    rng = np.random.default_rng(3)
    p = 0.2 + 0.6 * total / (total + 1.0)
    if regions == 1:
        return {"energies": e, "total": total, "scatter": scatter, "alphas": alphas,
                    "escape": 1 - p, "history": np.ones(start), "start": start}  # fmt: skip
    # rings: a random row-stochastic-like collision matrix with escapes
    base = rng.uniform(0.1, 1.0, (regions, regions))
    collide = (base / base.sum(axis=1, keepdims=True))[None] * p[:, None, None]
    escape = 1.0 - collide.sum(axis=2)
    return {"energies": e, "total": total, "scatter": scatter, "alphas": alphas,
                "escape": escape, "history": np.ones((start, regions)), "start": start,
                "collide": collide, "volumes": np.array([1.0, 1.3, 0.7][:regions])}  # fmt: skip


@pytest.mark.parametrize("regions", [1, 3])
def test_slowing_down_march(regions):
    from fission.pincell.shielding import march

    kw = _march_problem(regions)
    outside = np.full(len(kw["energies"]), 1e-3)
    lv_py, lv_rs = np.zeros(len(outside)), np.zeros(len(outside))
    py = march(**kw, share=0.05, outside_absorption=outside, levels=lv_py)
    rs = march(
        **kw, share=0.05, outside_absorption=outside, levels=lv_rs, use_rust=True
    )
    np.testing.assert_allclose(rs, py, rtol=1e-12, atol=0)
    np.testing.assert_allclose(lv_rs, lv_py, rtol=1e-12, atol=0)
    assert np.min(py) < 0.5  # the resonances really dip the flux
