"""A 17x17 lattice with MOC: k and pin powers against VERA problem 2 (step C4).

    uv run python -m fission.pincell.assembly 2A --rings --sph

The lattice calculation is the pin cell's (``heterogeneous.run``) on the
whole assembly: every fuel rod shares one shielded library per pellet ring
(one escape probability for all of them, so rods next to a guide tube or
at the edge are shielded as the average rod). Pin powers are fission rates
per rod normalized to a fuel average of 1, as in the KENO reference.
"""

from __future__ import annotations

import argparse
import time
import tomllib
from dataclasses import dataclass

import numpy as np

from fission.model import MODELS_DIR

from .heterogeneous import Result, run

PIN_POWERS = MODELS_DIR / "vera-p2-pin-powers.toml"


def reference_pin_powers(state: str, size: int = 17) -> np.ndarray:
    """KENO pin powers unfolded from the octant (NaN at the tubes)."""
    data = tomllib.loads(PIN_POWERS.read_text())[state]
    rows = [line.split() for line in data["octant"].strip().splitlines()]
    centre = size // 2
    out = np.full((size, size), np.nan)
    for i in range(size):
        for j in range(size):
            r, c = abs(i - centre), abs(j - centre)
            value = rows[max(r, c)][min(r, c)]
            if value != ".":
                out[i, j] = float(value)
    return out


def pin_powers(r: Result) -> np.ndarray:
    """Fission rate per cell / fuel-cell average, NaN in cells without fuel."""
    rates = r.solution.cell_rates("fission")
    out = np.where(rates > 0, rates, np.nan)
    return out / np.nanmean(out)


@dataclass
class PinComparison:
    rms: float  # % of the pin power
    max_abs: float  # %
    worst: tuple[int, int]  # (row, column) of the largest difference
    ours_max: float
    ref_max: float
    ours_min: float
    ref_min: float


def compare_pins(ours: np.ndarray, ref: np.ndarray) -> PinComparison:
    if not np.array_equal(np.isnan(ours), np.isnan(ref)):
        raise ValueError("fuel positions differ from the reference map")
    diff = 100 * (ours - ref) / ref
    worst = np.unravel_index(np.nanargmax(np.abs(diff)), diff.shape)
    return PinComparison(
        float(np.sqrt(np.nanmean(diff**2))), float(np.nanmax(np.abs(diff))),
        (int(worst[0]), int(worst[1])), float(np.nanmax(ours)), float(np.nanmax(ref)),
        float(np.nanmin(ours)), float(np.nanmin(ref)),
    )  # fmt: skip


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(prog="python -m fission.pincell.assembly", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)  # fmt: skip
    ap.add_argument("states", nargs="+")
    ap.add_argument("--model", default="vera-p2")
    ap.add_argument("--lattice", default="assembly")
    ap.add_argument("--shielding", default="sd")
    ap.add_argument("--rings", action="store_true")
    ap.add_argument("--sph", action="store_true")
    ap.add_argument(
        "--clad", action="store_true", help="clad slows down with the rings"
    )
    ap.add_argument("--moving", action="store_true",
                   help="U-238 scatters on moving targets (sd)")  # fmt: skip
    # 0.05 cm tracks are thicker than the 0.04 cm gap: -127 pcm on 2A
    ap.add_argument("--azim", type=int, default=64)
    ap.add_argument("--spacing", type=float, default=0.03)
    # the pin cell's converged defaults (heterogeneous.AZIM, FUEL_RINGS) cost
    # too much on 17x17; the polar set is Tabuchi-Yamamoto either way
    ap.add_argument("--fuel-rings", type=int, default=3)
    args = ap.parse_args(argv)
    for state in args.states:
        t0 = time.perf_counter()
        r = run(state, args.shielding, args.model, args.lattice, args.azim, args.spacing,
                fuel_rings=args.fuel_rings, ring_shielding=args.rings, sph=args.sph,
                coupled_clad=args.clad, moving_targets=args.moving)  # fmt: skip
        sol, ref = r.solution, r.state.reference
        versus = (
            f" (KENO {ref.keff:.6f}: {1e5 * (sol.k - ref.keff):+.0f} pcm)"
            if ref
            else ""
        )
        print(f"{args.model} {state}: k-infinity = {sol.k:.5f}{versus}, "
              f"{sol.problem.geometry.regions} regions, {sol.iterations} iterations")  # fmt: skip
        try:
            c = compare_pins(pin_powers(r), reference_pin_powers(state))
        except KeyError:
            c = None
        if c is not None:
            print(f"  pin powers vs KENO: RMS {c.rms:.2f}%, max {c.max_abs:.2f}% at "
                  f"row {c.worst[0]}, column {c.worst[1]}; max {c.ours_max:.4f} "
                  f"(KENO {c.ref_max:.4f}), min {c.ours_min:.4f} (KENO {c.ref_min:.4f})")  # fmt: skip
        print(f"  ({time.perf_counter() - t0:.0f} s)", flush=True)


if __name__ == "__main__":
    main()
