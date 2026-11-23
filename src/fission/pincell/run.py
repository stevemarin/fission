"""k-infinity of a homogenized lattice from a model file (VERA Problem 1).

    uv run python -m fission.pincell.run 1A --shielding nr
    uv run python -m fission.pincell.run 1A --model models/vera-p1.toml --lattice cell

The cases are the model's states. Shielding options:

* none: infinitely dilute libraries (chapter 10); U-238 resonance absorption
  is grossly overestimated, so k-infinity is far too low.
* nr: every library is re-averaged with the narrow-resonance flux of the
  mixture, phi(E) = w(E) / Sigma_t(E), where Sigma_t is the pointwise total
  macroscopic cross section of the homogeneous mixture. This is
  self-shielding in an infinite homogeneous medium.

The KENO-VI reference of the specification is for the *heterogeneous* pin
cell (fuel, clad and water separated), which has less resonance absorption
than the homogenized mixture, so it is printed for orientation only.
"""

from __future__ import annotations

import argparse
import time

import numpy as np

from fission.model import Component, Model, State, homogenize, load_model
from fission.processing.cache import cached_doppler
from fission.processing.groupr import GroupLibrary
from fission.processing.library import cached_groups

from .infinite import Mixture, Solution, solve

DEFAULT_MODEL = "vera-p1"
DEFAULT_LATTICE = "cell"

# Nuclides whose ENDF/B-VIII.1 resonance data this code cannot process yet:
# Fe-54 and Fe-57 (R-matrix limited format, LRF = 7), Fe-58 (an unresolved-
# range form, LRF = 1), Cr-53 (Reich-Moore with negative J selecting the
# channel spin). In the VERA cladding they are 1.9e-5 /b-cm (0.04% of its
# atoms); they are left out of the mixture.
UNSUPPORTED = frozenset({"Fe57", "Fe58", "Cr53", "Mn55"})


def components(model: Model, lattice: str, state: str) -> list[Component]:
    """The homogenized lattice at a state, without the unsupported nuclides."""
    return [
        c for c in homogenize(model, lattice, state) if c.nuclide not in UNSUPPORTED
    ]


def nr_shape(components: list[Component]) -> tuple[list[float], list[float]]:
    """(E, 1/Sigma_t(E)) of the mixture on the union of all pointwise grids."""
    tables = [(c, cached_doppler(c.nuclide, c.temperature)) for c in components]
    grid = np.unique(np.concatenate([np.asarray(t.energies) for _, t in tables]))
    total = np.zeros_like(grid)
    for c, t in tables:
        total += c.density * np.interp(
            grid, np.asarray(t.energies), np.asarray(t.xs[1])
        )
    return grid.tolist(), (1.0 / total).tolist()


def libraries(
    name: str, components: list[Component], shielding: str
) -> list[GroupLibrary]:
    if shielding == "none":
        return [
            cached_groups(c.nuclide, c.temperature, thermal=c.thermal)
            for c in components
        ]
    if shielding != "nr":
        raise ValueError(f"unknown shielding {shielding!r}")
    shape = nr_shape(components)
    return [
        cached_groups(c.nuclide, c.temperature, thermal=c.thermal, flux_shape=shape,
                      shape_name=f"NR-{name}")
        for c in components
    ]  # fmt: skip


def run(
    state: str,
    shielding: str,
    model: Model | str = DEFAULT_MODEL,
    lattice: str = DEFAULT_LATTICE,
) -> Solution:
    m = load_model(model) if isinstance(model, str) else model
    comps = components(m, lattice, state)
    return solve(Mixture(comps, libraries(state, comps, shielding)))


def report(model: Model, state: State, sol: Solution, shielding: str) -> str:
    m = sol.mixture
    ff = sol.four_factors()
    ref = state.reference
    lines = [
        (
            f"{model.name}, state {state.name}, homogenized, shielding = {shielding}: "
            f"fuel {state.fuel_temperature:g} K, moderator {state.moderator_temperature:g} K"
        ),
        f"  k-infinity = {sol.k:.5f}"
        + (f"   (reference {ref.keff:.5f}: {ref.source})" if ref else ""),
        "  four factors: " + ", ".join(f"{k} = {v:.4f}" for k, v in ff.items()),
        "  absorption per fission neutron (removal for (n,2n)):",
    ]
    phi = sol.flux
    rows = []
    for c, lib in zip(m.components, m.libraries, strict=True):
        a = c.density * (lib.sigma[1] - lib.scatter[0].sum(axis=1)) @ phi
        rows.append((a, f"{c.nuclide}({c.region})"))
    for a, name in sorted(rows, reverse=True)[:8]:
        lines.append(f"    {name:16s} {a:8.4f}")
    lines.append(
        f"  balance: production/removal = {(m.nu_fission @ phi) / (m.removal @ phi):.6f}"
    )
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(
        prog="python -m fission.pincell.run",
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("states", nargs="+", help="states of the model (e.g. 1A 1B)")
    ap.add_argument("--shielding", choices=["none", "nr"], default="nr")
    ap.add_argument(
        "--model", default=DEFAULT_MODEL, help="TOML file or name in models/"
    )
    ap.add_argument("--lattice", default=DEFAULT_LATTICE)
    args = ap.parse_args(argv)
    model = load_model(args.model)
    for name in args.states:
        if name not in model.states:
            ap.error(
                f"state {name!r} not in {model.name} (have {', '.join(model.states)})"
            )
    if args.lattice not in model.lattices:
        ap.error(f"lattice {args.lattice!r} not in {model.name}")
    for name in args.states:
        t0 = time.perf_counter()
        sol = run(name, args.shielding, model, args.lattice)
        print(report(model, model.states[name], sol, args.shielding))
        print(f"  ({time.perf_counter() - t0:.0f} s)\n")


if __name__ == "__main__":
    main()
