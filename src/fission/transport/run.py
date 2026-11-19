"""k-eigenvalue MOC calculation of a model's core (C5G7 by default).

    uv run python -m fission.transport.run c5g7 --azim 32 --spacing 0.05

For C5G7 the result is compared with the benchmark's multigroup MCNP
reference (NEA/NSC/DOC(2003)16, Table 3 eigenvalue, Table 4 extreme pin
powers, Table 5 assembly powers), with pin powers normalized to an average
of 1 over the 1056 fuel pins.
"""

from __future__ import annotations

import argparse
import time

import numpy as np

from fission.model import Model, load_model

from .geometry import Geometry
from .moc import Solution, build, solve
from .tracks import cyclic_tracks

# NEA/NSC/DOC(2003)16: reference multigroup MCNP solution of the 2D problem
C5G7_REFERENCE = {
    "k": 1.186550,  # Table 3 (68% confidence: 0.003%)
    "max_pin": 2.498,  # Table 4 (+-0.16%)
    "min_pin": 0.23,  # Table 4 (+-0.58%)
    "inner_uo2": 492.8,  # Table 5 (+-0.10%)
    "mox": 211.7,  # Table 5 (+-0.18%)
    "outer_uo2": 139.8,  # Table 5 (+-0.20%)
}


def fuel_pins(model: Model, geometry: Geometry) -> np.ndarray:
    """Mask of cells (row 0 at the top) whose pin contains fissile fuel,
    excluding fission chambers (fission rate tallied but not a fuel pin)."""
    mask = np.zeros(len(geometry.cells), dtype=bool)
    for i, cell in enumerate(geometry.cells):
        inner = model.materials[cell.materials[0]]
        name = inner.macroscopic or ""
        xs = model.multigroup.get(name)
        mask[i] = xs is not None and xs.fissile and "chamber" not in name
    return mask.reshape(geometry.ny, geometry.nx)[::-1]


def pin_powers(sol: Solution, model: Model) -> np.ndarray:
    """Fission rate per pin, normalized to an average of 1 over fuel pins
    (NaN elsewhere); row 0 at the top, as in the benchmark's Figure 3."""
    rates = sol.cell_rates("fission")
    fuel = fuel_pins(model, sol.problem.geometry)
    out = np.full(rates.shape, np.nan)
    out[fuel] = rates[fuel] / rates[fuel].mean()
    return out


def c5g7_report(sol: Solution, model: Model) -> str:
    ref = C5G7_REFERENCE
    p = pin_powers(sol, model)
    n = 17
    blocks = {
        "inner_uo2": p[:n, :n],
        "mox": p[:n, n : 2 * n],
        "outer_uo2": p[n : 2 * n, n : 2 * n],
    }
    k_pcm = (sol.k - ref["k"]) * 1e5
    lines = [
        f"  k = {sol.k:.6f}   reference {ref['k']:.6f}   difference {k_pcm:+.0f} pcm",
        (
            f"  max pin power {np.nanmax(p):.3f} (reference {ref['max_pin']:.3f}, "
            f"{100 * (np.nanmax(p) / ref['max_pin'] - 1):+.2f}%), "
            f"min {np.nanmin(p):.3f} (reference {ref['min_pin']:.2f})"
        ),
    ]
    for name, block in blocks.items():
        power = np.nansum(block)
        lines.append(
            f"  {name:9s} assembly power {power:6.1f} (reference {ref[name]:.1f}, "
            f"{100 * (power / ref[name] - 1):+.2f}%)"
        )
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(prog="python -m fission.transport.run", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)  # fmt: skip
    ap.add_argument("model", nargs="?", default="c5g7")
    ap.add_argument("--azim", type=int, default=32, help="directions in [0, 2 pi)")
    ap.add_argument("--spacing", type=float, default=0.05, help="track spacing (cm)")
    ap.add_argument("--polar", type=int, default=3, help="polar angles in (0, pi/2)")
    ap.add_argument("--rings", type=int, default=2, help="equal-area rings per pin")
    ap.add_argument("--sectors", type=int, default=8)
    ap.add_argument("--tol", type=float, default=1e-6, help="eigenvalue tolerance")
    args = ap.parse_args(argv)
    model = load_model(args.model)
    t0 = time.perf_counter()
    pin_materials = {
        m for p in model.pins.values() for m in p.materials[: len(p.radii)]
    }
    geometry = Geometry.from_core(
        model, dict.fromkeys(pin_materials, args.rings), args.sectors
    )
    assert model.core is not None
    tracks = cyclic_tracks(geometry.width, geometry.height, args.azim, args.spacing,
                           model.core.boundary)  # fmt: skip
    problem = build(model, geometry, tracks, args.polar)
    t1 = time.perf_counter()
    print(f"{model.name}: {geometry.regions} regions, {len(tracks)} tracks, "
          f"{len(problem.segments.length)} segments, {args.polar} polar angles "
          f"(set-up {t1 - t0:.1f} s)")  # fmt: skip
    sol = solve(problem, tol_k=args.tol, tol_source=10 * args.tol, verbose=True)
    print(f"  {sol.iterations} iterations in {time.perf_counter() - t1:.0f} s; "
          f"balance production/(absorption + leakage) = {sol.balance():.6f}")  # fmt: skip
    if model.name.startswith("C5G7"):
        print(c5g7_report(sol, model))
    else:
        print(f"  k = {sol.k:.6f}")


if __name__ == "__main__":
    main()
