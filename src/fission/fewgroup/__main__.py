"""Two-group constants of a model pin cell for the nodal code.

    uv run python -m fission.fewgroup 1A --rings --sph

Runs the heterogeneous MOC pin cell (chapter 16-17 options), homogenizes it,
finds the critical B1 spectrum, condenses to two groups (0.625 eV) and writes
endf-data/fewgroup/<model>-<state>.json with the keyword arguments of the
nodal code's ``Material``.
"""

from __future__ import annotations

import argparse
import time
from pathlib import Path

import numpy as np

from fission.model import load_model
from fission.pincell.heterogeneous import AZIM, FUEL_RINGS, SPACING, run

from . import condense, critical, homogenize, mode
from .condense import THERMAL_CUT

OUT = Path("endf-data/fewgroup")


def table(fg) -> str:
    rows = [
        ("D (cm)", fg.D),
        ("Sigma_a (1/cm)", fg.sigma_a),
        ("nu Sigma_f (1/cm)", fg.nu_sigma_f),
        ("kappa Sigma_f (MeV/cm)", fg.kappa_sigma_f / 1e6),
        ("chi", fg.chi),
        ("discontinuity factor", fg.adf),
        ("flux (relative)", fg.flux / fg.flux.sum()),
    ]
    lines = [
        f"    {'':24s}" + "".join(f"{f'group {g + 1}':>13s}" for g in range(len(fg.D)))
    ]
    for name, v in rows:
        lines.append(f"    {name:24s}" + "".join(f"{x:13.6g}" for x in v))
    for g in range(len(fg.D)):
        lines.append(
            f"    {f'Sigma_s {g + 1} -> g':24s}"
            + "".join(f"{x:13.6g}" for x in fg.sigma_s[g])
        )
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(prog="python -m fission.fewgroup", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)  # fmt: skip
    ap.add_argument("states", nargs="+")
    ap.add_argument("--model", default="vera-p1")
    ap.add_argument("--lattice", default=None, help="default: the model's only lattice")
    ap.add_argument("--azim", type=int, default=AZIM)
    ap.add_argument("--spacing", type=float, default=SPACING)
    ap.add_argument("--fuel-rings", type=int, default=FUEL_RINGS)
    ap.add_argument("--shielding", default="sd")
    ap.add_argument("--rings", action="store_true")
    ap.add_argument("--sph", action="store_true")
    ap.add_argument(
        "--clad", action="store_true", help="clad slows down with the rings"
    )
    ap.add_argument("--moving", action="store_true",
                   help="U-238 scatters on moving targets (sd)")  # fmt: skip
    ap.add_argument("--spectrum", choices=["b1", "infinite"], default="b1")
    ap.add_argument("--cut", type=float, nargs="*", default=[THERMAL_CUT],
                    help="inner few-group boundaries (eV)")  # fmt: skip
    args = ap.parse_args(argv)
    model = load_model(args.model)
    for state in args.states:
        t0 = time.perf_counter()
        lattice = args.lattice or next(iter(model.lattices))
        r = run(state, args.shielding, model, lattice, args.azim, args.spacing,
                fuel_rings=args.fuel_rings,
                ring_shielding=args.rings, sph=args.sph,
                coupled_clad=args.clad, moving_targets=args.moving)  # fmt: skip
        c = homogenize(r, f"{model.name} {state}")
        infinite = mode(c, 0.0)
        weight = critical(c) if args.spectrum == "b1" else infinite
        fg = condense(c, weight, tuple(args.cut))
        path = OUT / f"{args.model}-{state}.json"
        fg.write(path, state=state, lattice_k_infinity=r.solution.k,
                 options=vars(args) | {"states": None})  # fmt: skip
        m2 = (infinite.k - 1.0) / weight.b2 if weight.b2 else float("nan")
        print(f"{model.name} {state}: lattice k-inf {r.solution.k:.5f}, "
              f"spectrum {args.spectrum}, B^2 = {weight.b2:.6g} /cm^2 "
              f"(M^2 = {m2:.1f} cm^2)")  # fmt: skip
        print(f"  few-group k-inf {fg.k(0.0):.5f}, k at B^2 {fg.k(weight.b2):.6f}")
        print(table(fg))
        print(f"  -> {path} ({time.perf_counter() - t0:.0f} s)", flush=True)
        assert np.isfinite(fg.D).all()


if __name__ == "__main__":
    main()
