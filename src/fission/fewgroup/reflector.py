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
from fission.transport import moc

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


def homogeneous_surface_flux(
    fg: FewGroup, average: np.ndarray, width: float, k: float = np.inf
) -> np.ndarray:
    """Flux at x = 0 of the homogeneous medium 0 < x < width with
    reflective x = width and the given average flux (per group).

    -D phi'' + (A - F / k) phi = 0 with the removal matrix
    A[g, g] = Sigma_a,g + out-scattering, A[g, h] = -Sigma_h->g, and the
    fission matrix F[g, h] = chi_g nu Sigma_f,h (zero in a reflector). With
    D^-1 (A - F / k) = V diag(kappa^2) V^-1 the solution is
    phi = V diag(cosh kappa_i (width - x)) a, whose average is
    V diag(sinh(kappa W) / (kappa W)) a, so

        phi(0) = V diag(kappa W coth kappa W) V^-1 phi_average.

    In fuel the fundamental mode has kappa^2 < 0: kappa is imaginary and
    the cosh is a cosine; complex arithmetic covers both. Fixing the
    average is the same as fixing the current at x = 0: by the medium's
    balance it is J = (A - F / k) phi_average W.
    """
    s = fg.sigma_s
    a = np.diag(fg.sigma_a + s.sum(axis=1) - np.diag(s)) - (s.T - np.diag(np.diag(s)))
    a = a - np.outer(fg.chi, fg.nu_sigma_f) / k
    kappa2, v = np.linalg.eig(np.diag(1.0 / fg.D) @ a)
    kw = np.sqrt(kappa2.astype(complex)) * width
    return (v @ np.diag(kw / np.tanh(kw)) @ np.linalg.solve(v, average)).real


def _condensed(fine: FineConstants, fg: FewGroup, x: np.ndarray) -> np.ndarray:
    """Sum of a fine-group quantity over each few group (fine group g spans
    [b[g+1], b[g]], boundaries descending)."""
    edges = np.asarray(fg.boundaries)
    b = np.asarray(fine.groups.boundaries)
    out = np.zeros(len(fg.D))
    for big in range(len(fg.D)):
        inside = (b[:-1] <= edges[big] * (1 + 1e-9)) & (
            b[1:] >= edges[big + 1] * (1 - 1e-9)
        )
        out[big] = x[inside].sum()
    return out


def reflector_constants(result, name: str, cuts=(THERMAL_CUT,)) -> FewGroup:
    """Constants of the non-fuel cells of a one-row slab (fuel cells first),
    with the discontinuity factor of the fuel-reflector face.

    The MOC flux phi_s on the face is shared by both sides. Each side, as a
    homogeneous diffusion medium with its own average flux (and the slab's
    k for the fuel), has its own flux there (``homogeneous_surface_flux``),
    so f_refl = phi_s / phi_refl(0) and f_fuel = phi_s / phi_fuel(W_fuel).
    The nodal interface condition sees only their ratio, and the core's
    fuel nodes carry the ADF of their infinite lattice (~1), so the
    reflector gets f_refl / f_fuel: the fuel side of the face is then right
    for fuel next to the reflector (K. S. Smith's equivalence theory)."""
    geometry = result.solution.problem.geometry
    regions = water_cells(geometry)
    fine = homogenize(result, name, regions)
    fg = condense(fine, transport_spectrum(fine), cuts)
    if geometry.ny != 1:
        return fg
    cell = geometry.cell_of_regions()
    first = int(cell[regions].min())  # the reflector's first cell; its left face
    if cell[~regions].max() != first - 1 or cell[regions].max() != geometry.nx - 1:
        raise ValueError("the slab needs fuel cells first, then the reflector")
    fuel_fine = homogenize(result, "fuel", ~regions)
    fuel = condense(fuel_fine, transport_spectrum(fuel_fine), cuts)
    surface = _condensed(fine, fg, moc.face_flux(result.solution, [first])[0])
    pitch = result.solution.problem.tracks.width / geometry.nx
    phi_refl = homogeneous_surface_flux(
        fg, _condensed(fine, fg, fine.flux), (geometry.nx - first) * pitch
    )
    phi_fuel = homogeneous_surface_flux(
        fuel, _condensed(fuel_fine, fuel, fuel_fine.flux), first * pitch,
        result.solution.k,
    )  # fmt: skip
    fg.adf = phi_fuel / phi_refl  # = f_refl / f_fuel
    print(f"  face flux / homogeneous: reflector side {surface / phi_refl}, "
          f"fuel side {surface / phi_fuel}")  # fmt: skip
    return fg


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
              f"Sigma_a {fg.sigma_a.round(5)}, Sigma_1->2 {fg.sigma_s[0, 1]:.5f}, "
              f"ADF {fg.adf.round(4)} "
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
