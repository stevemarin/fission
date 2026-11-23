"""Solve a model's core at power with thermal-hydraulic feedback.

    uv run python -m fission.feedback wbn1-core

Prints the iterations, k-effective, the radial assembly powers, the axial
power shape and the temperatures at the converged state.
"""

from __future__ import annotations

import argparse

import numpy as np

from fission.model import load_model
from fission.nodal.__main__ import assembly_powers

from .core import isothermal, solve


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(prog="python -m fission.feedback", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)  # fmt: skip
    ap.add_argument("model", nargs="?", default="wbn1-core")
    ap.add_argument("--refine", type=int, default=1)
    ap.add_argument("--axial", type=float, default=None, help="largest axial node (cm)")
    ap.add_argument("--power", type=float, default=1.0, help="fraction of rated power")
    ap.add_argument("--itc", action="store_true",
                    help="isothermal temperature coefficient at zero power, 560-570 K")  # fmt: skip
    args = ap.parse_args(argv)
    model = load_model(args.model)
    if args.itc:
        temps = (560.0, 565.0, 570.0)
        ks = isothermal(model, temps, args.refine, args.axial)
        for t, k in zip(temps, ks, strict=True):
            print(f"{model.name}: T = {t:.1f} K (fuel and water), k = {k:.6f}")
        itc = 1e5 * (1 / ks[0] - 1 / ks[-1]) / (temps[-1] - temps[0])
        print(
            f"  isothermal temperature coefficient: {itc:+.2f} pcm/K = {itc / 1.8:+.2f} pcm/F"
        )
        return
    r = solve(model, args.refine, args.axial, power_level=args.power, verbose=True)
    mesh = r.mesh
    fuel = mesh.field(lambda x: x.nu_fission).sum(axis=0) > 0
    print(f"{model.name}: k = {r.k:.5f} after {len(r.history)} feedback iterations "
          f"({r.seconds:.0f} s)")  # fmt: skip
    tf, tm = r.states["fuel_temperature"], r.states["moderator_temperature"]
    v = mesh.volumes()
    print(f"  fuel (Doppler) temperature: average {np.sum((tf * v)[fuel]) / np.sum(v[fuel]):.1f} K, "
          f"max {tf[fuel].max():.1f} K; peak centreline {np.nanmax(r.centre):.1f} K")  # fmt: skip
    if np.isfinite(r.film).any():
        print(f"  clad-to-water heat transfer: {np.nanmin(r.film) / 1e3:.1f}-"
              f"{np.nanmax(r.film) / 1e3:.1f} kW/(m2 K)")  # fmt: skip
    outlet = tm[-1][fuel.any(axis=0)]
    print(f"  water: inlet {model.core.thermal.inlet_temperature:.1f} K, outlet "  # type: ignore[union-attr]
          f"{outlet.min():.1f}-{outlet.max():.1f} K, density "
          f"{r.states['moderator_density'][fuel].min():.4f}-"
          f"{r.states['moderator_density'][fuel].max():.4f} g/cm3")  # fmt: skip
    p = assembly_powers(mesh, r.flux)
    print("  assembly powers (quarter core, centre top left):")
    for row in p:
        print("    " + " ".join("  .   " if np.isnan(x) else f"{x:6.3f}" for x in row))
    axial = r.power.sum(axis=(1, 2))
    dz = mesh.widths[2]
    rows = axial > 0
    # power per cm, normalized to 1 over the fuel height (nodes differ in height)
    shape = axial[rows] / dz[rows]
    shape /= axial[rows].sum() / dz[rows].sum()
    # offset about the midplane of the fuel; the node across it is split by height
    upper = np.cumsum(dz)[rows]
    lower = upper - dz[rows]
    half = 0.5 * (lower[0] + upper[-1])
    above = np.clip((upper - half) / dz[rows], 0.0, 1.0)
    top = np.sum(axial[rows] * above)
    bottom = np.sum(axial[rows] * (1.0 - above))
    print("  axial power (bottom to top): " + " ".join(f"{x:.2f}" for x in shape))
    print(
        f"  axial offset (top - bottom) / total: {100 * (top - bottom) / (top + bottom):+.1f}%"
    )


if __name__ == "__main__":
    main()
