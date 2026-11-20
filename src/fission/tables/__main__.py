"""Build the cross-section table of a state matrix and test it.

    uv run python -m fission.tables vera

Reads endf-data/branches/<model file>-<branches>-<solver>.json (from
``make branches``), writes endf-data/tables/<same name>.json for the nodal
code, and reports three tests:

1. leave one out: drop an interior node of an axis, interpolate there from
   the remaining nodes, compare with the computed constants (fuel axis in
   T and in sqrt(T));
2. the polynomial fit: terms chosen by leave-one-out error;
3. check points (if ``-check.json`` exists, see ``branches.run_states``):
   states off the grid, computed by the lattice code, against the
   interpolated and the fitted constants.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from fission.branches import OUT as BRANCHES_OUT
from fission.branches import Table
from fission.model import BRANCH_AXES

from . import KEYS, XSTable, describe, k_infinity, select

OUT = Path("endf-data/tables")


def errors(constants: dict[str, np.ndarray], record: dict) -> dict[str, float]:
    """Largest relative error of each constant, and k-infinity in pcm."""
    out = {}
    for k, v in constants.items():
        ref = np.asarray(record[k], dtype=float)
        scale = np.maximum(np.abs(ref), 1e-12 * np.abs(ref).max())
        mask = np.abs(ref) > 1e-12 * np.abs(ref).max()  # skip structural zeros
        out[k] = (
            float(np.max(np.abs(v - ref)[mask] / scale[mask])) if mask.any() else 0.0
        )
    out["k (pcm)"] = 1e5 * (k_infinity(constants) / record["few_group_k_infinity"] - 1)
    return out


def row(label: str, e: dict[str, float]) -> str:
    worst = max((v, k) for k, v in e.items() if k != "k (pcm)")
    return (f"  {label:34s} k {e['k (pcm)']:+7.1f} pcm   worst constant "
            f"{100 * worst[0]:6.3f}% ({worst[1]})")  # fmt: skip


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(prog="python -m fission.tables", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)  # fmt: skip
    ap.add_argument("branches")
    ap.add_argument("--model", default="vera-p1")
    ap.add_argument("--solver", default="infinite")
    args = ap.parse_args(argv)
    stem = f"{Path(args.model).stem}-{args.branches}-{args.solver}"
    table = Table.read(BRANCHES_OUT / f"{stem}.json")
    if not table.complete:
        ap.error(f"{stem}: the state matrix is not complete (make branches)")
    sqrt_t = XSTable.from_branches(table)
    linear_t = XSTable.from_branches(table, transforms={})
    path = OUT / f"{stem}.json"
    sqrt_t.write(path)
    print(f"{stem}: {int(np.prod(sqrt_t.shape))} nodes {sqrt_t.shape} -> {path}")

    print("1. leave one out (interior nodes, interpolated from the rest):")
    for n, axis in enumerate(BRANCH_AXES):
        values = sqrt_t.axes[axis]
        for i in range(1, len(values) - 1):
            for label, xs in (("sqrt(T)", sqrt_t), ("T", linear_t)):
                if axis != "fuel_temperature" and label == "T":
                    continue
                reduced = xs.without(axis, i)
                worst: dict[str, float] = {}
                for index, record in table.points.items():
                    if index[n] != i:
                        continue
                    e = errors(reduced.constants(record["state"]), record)
                    for k, v in e.items():
                        worst[k] = max(worst.get(k, 0.0), abs(v))
                name = f"{axis} = {values[i]:g}" + (
                    f", in {label}" if axis == "fuel_temperature" else ""
                )
                print(row(name, worst))

    print("2. polynomial fit (terms chosen by leave-one-out error over the grid):")
    fit, steps = select(sqrt_t)
    for term, score in steps:
        mark = " <" if len(fit.terms) == steps.index((term, score)) + 1 else ""
        print(f"     + {describe(term):9s} relative LOO error {score:.2e}{mark}")
    print(
        f"   {fit.size} terms per constant, against {int(np.prod(sqrt_t.shape))} nodes"
    )

    check = BRANCHES_OUT / f"{stem}-check.json"
    if not check.exists():
        print(f"3. no check points ({check} missing)")
        return
    records = json.loads(check.read_text())["points"]
    print(f"3. {len(records)} check points off the grid (max |error| over them):")
    for label, model in (("multilinear, fuel in sqrt(T)", sqrt_t),
                         ("multilinear, fuel in T", linear_t),
                         (f"fit, {fit.size} terms", fit)):  # fmt: skip
        worst = {}
        for record in records:
            e = errors(model.constants(record["state"]), record)
            for k, v in e.items():
                worst[k] = max(worst.get(k, 0.0), abs(v))
        print(row(label, worst))
    print("   per point, k error in pcm (multilinear sqrt(T) / fit):")
    for record in records:
        s = record["state"]
        a = errors(sqrt_t.constants(s), record)["k (pcm)"]
        b = errors(fit.constants(s), record)["k (pcm)"]
        print(f"     Tf {s['fuel_temperature']:6.0f} K, Tm {s['moderator_temperature']:4.0f} K, "
              f"rho {s['moderator_density']:.4f}, B {s['boron_ppm']:6.0f} ppm: "
              f"{a:+7.1f} / {b:+7.1f}")  # fmt: skip
    assert all(k in KEYS for k in sqrt_t.constants(records[0]["state"]))


if __name__ == "__main__":
    main()
