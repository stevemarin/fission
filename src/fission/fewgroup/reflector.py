"""Reflector constants from a fuel-reflector slab (M4b).

A reflector has no fission, so its spectrum is not its own: it is what
leaks in from the fuel. The water next to the core is full of fast neutrons
from the fuel slowing down, and its thermal flux peaks a few centimetres
out, where they have thermalized but are no longer absorbed by the fuel.
Constants for a reflector node therefore come from a calculation with fuel
next to the water: here one row of fuel cells and one of water cells
(``lattices.reflector`` of the model), solved by the heterogeneous MOC of
chapters 16-19. The water cells are homogenized with the slab's own flux
(no B1 spectrum: there is nothing to make critical), and the diffusion
coefficient is the transport one, D_g = 1 / (3 Sigma_tr,g) with
Sigma_tr = Sigma_t - sum_g' Sigma_s1,g->g', condensed so that the current
D_g phi_g of each fine group is kept.

    uv run python -m fission.fewgroup.reflector     # boron 0, 1300, 2000 ppm

The result is an M3 table (``tables.XSTable``) with boron as its only
axis, written to endf-data/tables/<model>-reflector.json.
"""

from __future__ import annotations

import argparse
import time
from dataclasses import replace
from pathlib import Path

import numpy as np

from fission.model import BRANCH_AXES, load_model

from .b1 import Mode
from .condense import THERMAL_CUT, FewGroup, condense
from .homogenize import FineConstants, homogenize

OUT = Path("endf-data/tables")


def water_cells(geometry, material: str = "fuel") -> np.ndarray:
    """Flat-source regions of the cells that contain no `material`."""
    cell = geometry.cell_of_regions()
    names = np.array(geometry.materials())
    has = np.zeros(geometry.nx * geometry.ny, dtype=bool)
    for c, n in zip(cell, names, strict=True):
        if n.split("/")[0] == material:
            has[c] = True
    return ~has[cell]


def transport_spectrum(c: FineConstants) -> Mode:
    """The slab's own flux as the weight; current D_g phi_g with the
    transport diffusion coefficient (B^2 = 0, k undefined)."""
    sigma_tr = c.total - c.scatter[1].sum(axis=1)
    d = 1.0 / (3.0 * sigma_tr)
    flux = c.flux / c.flux.sum()
    return Mode(0.0, float("nan"), flux, d * flux)


def reflector_constants(result, name: str, cuts=(THERMAL_CUT,)) -> FewGroup:
    regions = water_cells(result.solution.problem.geometry)
    fine = homogenize(result, name, regions)
    return condense(fine, transport_spectrum(fine), cuts)


def main(argv: list[str] | None = None) -> None:
    from fission.pincell.heterogeneous import run
    from fission.tables import XSTable

    ap = argparse.ArgumentParser(prog="python -m fission.fewgroup.reflector",
                                 description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)  # fmt: skip
    ap.add_argument("--model", default="vera-p1")
    ap.add_argument("--lattice", default="reflector")
    ap.add_argument("--state", default="1B", help="fuel and water conditions")
    ap.add_argument("--boron", type=float, nargs="+", default=[0.0, 1300.0, 2000.0])
    ap.add_argument("--azim", type=int, default=64)
    ap.add_argument("--spacing", type=float, default=0.03)
    ap.add_argument("--fuel-rings", type=int, default=3)
    args = ap.parse_args(argv)
    model = load_model(args.model)
    base = model.states[args.state]
    records = []
    for ppm in args.boron:
        t0 = time.perf_counter()
        state = replace(
            base, name=f"reflector-{ppm:g}ppm", boron_ppm=ppm, reference=None
        )
        m = replace(model, states={**model.states, state.name: state})
        r = run(state.name, "sd", m, args.lattice, args.azim, args.spacing,
                fuel_rings=args.fuel_rings, ring_shielding=True, sph=True,
                coupled_clad=True, cmfd=False)  # fmt: skip
        # (no CMFD: its coarse problem on a long fuel-water slab diverges)
        fg = reflector_constants(r, state.name)
        records.append(fg)
        print(f"boron {ppm:g} ppm: slab k {r.solution.k:.5f}, D {fg.D.round(4)}, "
              f"Sigma_a {fg.sigma_a.round(5)}, Sigma_1->2 {fg.sigma_s[0, 1]:.5f} "
              f"({time.perf_counter() - t0:.0f} s)", flush=True)  # fmt: skip
    axes: dict[str, tuple[float, ...]] = {a: (getattr(base, a),) for a in BRANCH_AXES}
    axes["boron_ppm"] = tuple(args.boron)
    n = len(records)
    data = {
        k: np.array([getattr(fg, k) for fg in records]).reshape((1, 1, 1, n) + getattr(records[0], k).shape)
        for k in ("D", "total", "sigma_a", "sigma_s", "nu_sigma_f", "sigma_f",
                  "kappa_sigma_f", "chi", "adf")
    }  # fmt: skip
    table = XSTable(axes, data, {}, {"model": model.name, "lattice": args.lattice,
                                     "solver": "moc", "region": "water cells",
                                     "boundaries_eV": list(records[0].boundaries)})  # fmt: skip
    path = OUT / f"{Path(args.model).stem}-reflector.json"
    table.write(path)
    print(f"-> {path}")


if __name__ == "__main__":
    main()
