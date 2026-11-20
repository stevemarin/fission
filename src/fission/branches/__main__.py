"""Run a state matrix of a model and write its few-group table.

    uv run python -m fission.branches vera
    uv run python -m fission.branches vera --solver moc --rings --sph --clad

The table goes to endf-data/branches/<model file>-<branches>-<solver>.json; an
interrupted run resumes from it (use --fresh to start over). --check N also
solves N random states inside the matrix (<...>-check.json), against which
``python -m fission.tables`` tests interpolation and fits.
"""

from __future__ import annotations

import argparse
import time
from pathlib import Path

from fission.fewgroup.condense import THERMAL_CUT
from fission.model import load_model

from . import OUT, SOLVERS, coefficients, run_matrix, run_states, sample_states

UNITS = {
    "fuel_temperature": "pcm/K",
    "moderator_temperature": "pcm/K",
    "moderator_density": "pcm per 0.01 g/cm3",
    "boron_ppm": "pcm/ppm",
}
SCALE = {"moderator_density": 0.01}


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(prog="python -m fission.branches", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)  # fmt: skip
    ap.add_argument("branches", help="name of a [branches.<name>] section")
    ap.add_argument("--model", default="vera-p1")
    ap.add_argument("--solver", choices=SOLVERS, default="infinite")
    ap.add_argument("--spectrum", choices=["b1", "infinite"], default="b1")
    ap.add_argument("--cut", type=float, nargs="*", default=[THERMAL_CUT],
                    help="inner few-group boundaries (eV)")  # fmt: skip
    ap.add_argument("--fresh", action="store_true", help="ignore an existing table")
    ap.add_argument("--keep-libraries", action="store_true",
                    help="keep each point's shaped libraries in the cache")  # fmt: skip
    ap.add_argument("--check", type=int, default=0, metavar="N",
                    help="also solve N random states off the grid (tests M3 tables)")  # fmt: skip
    moc = ap.add_argument_group("moc solver options")
    moc.add_argument("--rings", action="store_true")
    moc.add_argument("--sph", action="store_true")
    moc.add_argument("--clad", action="store_true")
    moc.add_argument("--moving", action="store_true")
    moc.add_argument("--azim", type=int, default=None)
    moc.add_argument("--spacing", type=float, default=None)
    moc.add_argument("--fuel-rings", type=int, default=None)
    args = ap.parse_args(argv)
    model = load_model(args.model)
    if args.branches not in model.branches:
        ap.error(f"no [branches.{args.branches}] in {model.name} "
                 f"(have {', '.join(model.branches) or 'none'})")  # fmt: skip
    options: dict = {}
    if args.solver == "moc":
        options = {"ring_shielding": args.rings, "sph": args.sph,
                   "coupled_clad": args.clad, "moving_targets": args.moving}  # fmt: skip
        for key, value in (("azim", args.azim), ("spacing", args.spacing),
                           ("fuel_rings", args.fuel_rings)):  # fmt: skip
            if value is not None:
                options[key] = value
    path = OUT / f"{Path(args.model).stem}-{args.branches}-{args.solver}.json"
    t0 = time.perf_counter()
    table = run_matrix(model, args.branches, args.solver, options, args.spectrum,
                       tuple(args.cut), path, args.fresh, args.keep_libraries,
                       progress=lambda line: print(line, flush=True))  # fmt: skip
    print(f"{len(table.points)} points -> {path} ({time.perf_counter() - t0:.0f} s)")
    if args.check:
        check = path.with_name(path.stem + "-check.json")
        states = sample_states(model.branches[args.branches], args.check)
        run_states(model, args.branches, states, args.solver, options, args.spectrum,
                   tuple(args.cut), check,
                   progress=lambda line: print(line, flush=True))  # fmt: skip
        print(f"{len(states)} check points -> {check}")
    base = table.header["base"]
    print(f"reactivity coefficients at {base['state']} (lattice k-infinity):")
    for c in coefficients(table):
        scale = SCALE.get(c.axis, 1.0)
        print(f"  {c.axis:22s} {c.value * scale:+9.3f} {UNITS[c.axis]}"
              f"   ({c.low:g} to {c.high:g})")  # fmt: skip


if __name__ == "__main__":
    main()
