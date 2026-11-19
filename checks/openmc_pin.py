"""Continuous-energy OpenMC run of a pin cell described by a JSON spec.

Independent check, run in the tools environment:

    python checks/openmc_pin.py spec.json result.json

The spec (written by ``fission.pincell.reference``) gives the pitch, the
ring radii and materials, every material's temperature, number densities and
S(alpha, beta) tables, the energy groups for the tallies and the run
settings. The fuel is split into equal-area rings so that reaction rates can
be compared ring by ring. Results (per source neutron, with 1-sigma
uncertainties) go to result.json.
"""

from __future__ import annotations

import json
import math
import sys
import tempfile
from pathlib import Path

import numpy as np
import openmc

LIBRARY = Path(__file__).resolve().parents[1] / "endf-data" / "openmc"


def build(spec: dict) -> tuple[openmc.Model, list[tuple[str, int, openmc.Cell]]]:
    mats = {}
    for name, m in spec["materials"].items():
        mat = openmc.Material(name=name, temperature=m["temperature"])
        for nuclide, density in m["densities"].items():
            mat.add_nuclide(nuclide, density)  # atoms/b-cm
        mat.set_density("sum")
        for table in m.get("sab", []):
            mat.add_s_alpha_beta(table)
        mats[name] = mat
    pitch = spec["pitch"]
    half = pitch / 2
    box = (
        +openmc.XPlane(-half, boundary_type="reflective")
        & -openmc.XPlane(half, boundary_type="reflective")
        & +openmc.YPlane(-half, boundary_type="reflective")
        & -openmc.YPlane(half, boundary_type="reflective")
    )
    # rings: [(outer radius, material, pieces)], the last material fills the box
    cells, inner, region_cells = [], None, []
    radii = spec["radii"]
    for (r_out, name), pieces in zip(
        zip(radii, spec["ring_materials"], strict=False), spec["pieces"], strict=False
    ):
        r_in = 0.0 if inner is None else inner.r
        for k in range(pieces):  # equal-area sub-rings
            r = math.sqrt(r_in**2 + (k + 1) / pieces * (r_out**2 - r_in**2))
            cyl = openmc.ZCylinder(r=r)
            region = -cyl if inner is None else (+inner & -cyl)
            cell = openmc.Cell(name=f"{name}/{k}", fill=mats[name], region=region)
            cells.append(cell)
            region_cells.append((name, k, cell))
            inner = cyl
    last = spec["ring_materials"][-1]
    cell = openmc.Cell(name=f"{last}/0", fill=mats[last], region=box & +inner)
    cells.append(cell)
    region_cells.append((last, 0, cell))
    geometry = openmc.Geometry(openmc.Universe(cells=cells))

    run = spec["settings"]
    settings = openmc.Settings()
    settings.run_mode = "eigenvalue"
    settings.particles = run["particles"]
    settings.batches = run["batches"]
    settings.inactive = run["inactive"]
    settings.seed = run.get("seed", 1)
    settings.ptables = run.get("ptables", True)
    settings.temperature = {"method": "nearest", "tolerance": 1.0}
    settings.source = openmc.IndependentSource(
        space=openmc.stats.Point(), constraints={"fissionable": True}
    )
    settings.output = {"tallies": False}

    groups = openmc.EnergyFilter(np.sort(np.asarray(spec["groups"])))
    cell_filter = openmc.CellFilter([c for _, _, c in region_cells])
    tallies = openmc.Tallies()
    flux = openmc.Tally(name="flux")
    flux.filters = [cell_filter, groups]
    flux.scores = ["flux"]
    tallies.append(flux)
    for name, mat in mats.items():
        t = openmc.Tally(name=f"rates {name}")
        t.filters = [
            openmc.CellFilter([c for m, _, c in region_cells if m == name]),
            groups,
        ]
        t.nuclides = list(spec["materials"][name]["densities"])
        t.scores = ["absorption", "(n,gamma)", "fission", "nu-fission"]
        tallies.append(t)
    materials = openmc.Materials(mats.values())
    materials.cross_sections = str(LIBRARY / "cross_sections.xml")
    model = openmc.Model(geometry, materials, settings, tallies)
    return model, region_cells


def main(spec_path: str, out_path: str) -> None:
    spec = json.loads(Path(spec_path).read_text())
    model, region_cells = build(spec)
    with tempfile.TemporaryDirectory() as tmp:
        statepoint = model.run(cwd=tmp, threads=spec["settings"].get("threads"))
        with openmc.StatePoint(statepoint) as sp:
            k = sp.keff
            out = {
                "k": k.nominal_value,
                "k_std": k.std_dev,
                "regions": [f"{m}/{i}" for m, i, _ in region_cells],
                "flux": sp.get_tally(name="flux")
                .mean.reshape(len(region_cells), -1)[:, ::-1]
                .tolist(),  # group 0 = highest energy
                "rates": {},
            }
            for name in spec["materials"]:
                t = sp.get_tally(name=f"rates {name}")
                cells = [f"{m}/{i}" for m, i, _ in region_cells if m == name]
                nuclides = list(spec["materials"][name]["densities"])
                shape = (len(cells), -1, len(nuclides), 4)
                mean = t.mean.reshape(shape)[:, ::-1]
                std = t.std_dev.reshape(shape)[:, ::-1]
                for c, cell in enumerate(cells):
                    for j, n in enumerate(nuclides):
                        for s, score in enumerate(t.scores):
                            out["rates"].setdefault(cell, {}).setdefault(n, {})[
                                score
                            ] = {
                                "mean": mean[c, :, j, s].tolist(),
                                "std": std[c, :, j, s].tolist(),
                            }
    Path(out_path).write_text(json.dumps(out))
    print(f"k = {out['k']:.5f} +- {out['k_std']:.5f} -> {out_path}")


if __name__ == "__main__":
    main(*sys.argv[1:3])
