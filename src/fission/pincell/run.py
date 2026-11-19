"""k-infinity of the homogenized VERA Problem 1 pin cell.

    uv run python -m fission.pincell.run 1A --shielding nr

Shielding options:

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

from fission.processing.cache import cached_doppler
from fission.processing.groupr import GroupLibrary
from fission.processing.library import cached_groups

from .infinite import Mixture, Solution, solve
from .vera import CASES, Case, Component, homogenized


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
    case: Case, components: list[Component], shielding: str
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
                      shape_name=f"NR-{case.name}")
        for c in components
    ]  # fmt: skip


def run(case_name: str, shielding: str) -> Solution:
    case = CASES[case_name]
    components = homogenized(case)
    return solve(Mixture(components, libraries(case, components, shielding)))


def report(case: Case, sol: Solution, shielding: str) -> str:
    m = sol.mixture
    ff = sol.four_factors()
    lines = [
        (
            f"VERA problem {case.name}, homogenized, shielding = {shielding}: "
            f"fuel {case.fuel_temperature:g} K, moderator {case.moderator_temperature:g} K"
        ),
        f"  k-infinity = {sol.k:.5f}   (heterogeneous CE KENO-VI reference: {case.keff:.5f})",
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
    ap.add_argument("cases", nargs="+", choices=sorted(CASES))
    ap.add_argument("--shielding", choices=["none", "nr"], default="nr")
    args = ap.parse_args(argv)
    for name in args.cases:
        t0 = time.perf_counter()
        sol = run(name, args.shielding)
        print(report(CASES[name], sol, args.shielding))
        print(f"  ({time.perf_counter() - t0:.0f} s)\n")


if __name__ == "__main__":
    main()
