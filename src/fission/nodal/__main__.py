"""Solve a model's core with nodal diffusion.

    uv run python -m fission.nodal iaea-3d                 # NEM, one node per assembly
    uv run python -m fission.nodal iaea-3d --method fd --refine 4

Prints k-effective and the radial assembly powers (average power density
of each assembly, normalized to 1 over the fuel assemblies of the full
core; power is taken proportional to nu Sigma_f phi).
"""

from __future__ import annotations

import argparse
import time

import numpy as np

from fission.model import load_model

from . import diffusion, nem
from .mesh import CoreMesh, core_mesh, parse_banks


def assembly_powers(mesh: CoreMesh, flux: np.ndarray) -> np.ndarray:
    """[rows, columns] of the map: average power density per position, NaN
    where there is no fuel; normalized to 1 over the full core."""
    power = (mesh.field(lambda x: x.nu_fission) * flux).sum(axis=0) * mesh.volumes()
    fuel = mesh.field(lambda x: x.nu_fission).sum(axis=0) > 0
    rows = int(mesh.position[..., 0].max()) + 1
    cols = int(mesh.position[..., 1].max()) + 1
    p = np.zeros((rows, cols))
    v = np.zeros((rows, cols))
    for j in range(mesh.shape[1]):
        for i in range(mesh.shape[2]):
            r, c = mesh.position[j, i]
            p[r, c] += power[:, j, i].sum()
            v[r, c] += (mesh.volumes()[:, j, i] * fuel[:, j, i]).sum()
    density = np.where(v > 0, p / np.where(v > 0, v, 1), np.nan)
    # multiplicity in the full core: quarter positions on the symmetry lines
    weight = np.ones((rows, cols))
    if mesh.symmetry == "quarter":
        weight[:] = 4.0
        weight[0, :] = weight[:, 0] = 2.0
        weight[0, 0] = 1.0
    mask = ~np.isnan(density)
    mean = np.sum((density * weight)[mask]) / np.sum(weight[mask])
    return density / mean


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(prog="python -m fission.nodal", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)  # fmt: skip
    ap.add_argument("model", nargs="?", default="iaea-3d")
    ap.add_argument("--method", choices=["nem", "fd"], default="nem")
    ap.add_argument("--refine", type=int, default=1, help="nodes per assembly width")
    ap.add_argument("--axial", type=float, default=None, help="largest axial node (cm)")
    ap.add_argument("--bank", action="append", metavar="NAME[=STEPS]",
                    help="insert a control rod bank (steps withdrawn, default 0 = fully in)")  # fmt: skip
    args = ap.parse_args(argv)
    model = load_model(args.model)
    mesh = core_mesh(model, args.refine, args.axial, parse_banks(args.bank))
    t0 = time.perf_counter()
    if args.method == "nem":
        r = nem.solve(mesh)
        k, flux, info = r.k, r.flux, f"{r.iterations} NEM iterations"
    else:
        s = diffusion.solve(mesh, tol_k=1e-8, tol_source=1e-7, omega=1.8, inner=6)
        k, flux, info = s.k, s.flux, f"{s.outer} outer iterations"
    nz, ny, nx = mesh.shape
    print(f"{model.name}: {args.method}, {nx} x {ny} x {nz} nodes, k = {k:.6f} "
          f"({info}, {time.perf_counter() - t0:.1f} s)")  # fmt: skip
    p = assembly_powers(mesh, flux)
    print("assembly powers (map orientation):")
    for row in p:
        print("  " + " ".join("  .   " if np.isnan(v) else f"{v:6.3f}" for v in row))
    print(f"  max {np.nanmax(p):.3f}")


if __name__ == "__main__":
    main()
