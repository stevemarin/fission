"""Continuous-energy OpenMC run of a pin lattice described by a JSON spec.

Independent check, run in the tools environment:

    python checks/openmc_pin.py spec.json result.json

The spec (written by ``fission.pincell.reference``) gives the pin pitch, the
map of pin types (rows from the top), each type's ring radii and materials,
the water gap around the lattice, every material's temperature, number
densities and S(alpha, beta) tables, the energy groups for the tallies and
the run settings. A pin cell is a 1 x 1 map. Rings marked with several
pieces (the fuel) are split into equal-area sub-rings so that reaction
rates can be compared ring by ring.

Results (per source neutron) go to result.json: k, the flux and the
reaction rates in the groups for every region label ("fuel/0" ... summed
over all pins; "<fill>/gap" for the gap strip around the lattice), and the
fission rate of every pin position.
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


def build(spec: dict) -> tuple[openmc.Model, list[tuple[str, openmc.Cell]]]:
    mats = {}
    for name, m in spec["materials"].items():
        mat = openmc.Material(name=name, temperature=m["temperature"])
        for nuclide, density in m["densities"].items():
            mat.add_nuclide(nuclide, density)  # atoms/b-cm
        mat.set_density("sum")
        for table in m.get("sab", []):
            mat.add_s_alpha_beta(table)
        mats[name] = mat
    labelled: list[tuple[str, openmc.Cell]] = []  # (label, cell) for the tallies
    universes = {}
    for pin_name, pin in spec["pins"].items():
        cells, inner = [], None
        for r_out, name, pieces in zip(pin["radii"], pin["materials"], pin["pieces"],
                                       strict=False):  # fmt: skip
            r_in = 0.0 if inner is None else inner.r
            for k in range(pieces):  # equal-area sub-rings
                r = math.sqrt(r_in**2 + (k + 1) / pieces * (r_out**2 - r_in**2))
                cyl = openmc.ZCylinder(r=r)
                region = -cyl if inner is None else (+inner & -cyl)
                cell = openmc.Cell(fill=mats[name], region=region)
                cells.append(cell)
                labelled.append((f"{name}/{k}", cell))
                inner = cyl
        fill = pin["materials"][-1]
        cell = openmc.Cell(fill=mats[fill], region=None if inner is None else +inner)
        cells.append(cell)
        labelled.append((f"{fill}/0", cell))
        universes[pin_name] = openmc.Universe(name=pin_name, cells=cells)
    pitch, n, gap = spec["pitch"], len(spec["map"]), spec.get("gap", 0.0)
    lattice = openmc.RectLattice()
    lattice.pitch = (pitch, pitch)
    lattice.lower_left = (-n * pitch / 2, -n * pitch / 2)
    lattice.universes = [[universes[p] for p in row] for row in spec["map"]]
    half_lattice, half = n * pitch / 2, n * pitch / 2 + gap

    def box(h, boundary="transmission"):
        return (
            +openmc.XPlane(-h, boundary_type=boundary)
            & -openmc.XPlane(h, boundary_type=boundary)
            & +openmc.YPlane(-h, boundary_type=boundary)
            & -openmc.YPlane(h, boundary_type=boundary)
        )

    outer = box(half, "reflective")
    root_cells = []
    if gap > 0:
        inside = box(half_lattice)
        root_cells.append(openmc.Cell(fill=lattice, region=inside))
        fill = spec["gap_material"]
        gap_cell = openmc.Cell(fill=mats[fill], region=outer & ~inside)
        root_cells.append(gap_cell)
        labelled.append((f"{fill}/gap", gap_cell))
    else:
        root_cells.append(openmc.Cell(fill=lattice, region=outer))
    geometry = openmc.Geometry(openmc.Universe(cells=root_cells))

    run = spec["settings"]
    settings = openmc.Settings()
    settings.run_mode = "eigenvalue"
    settings.particles = run["particles"]
    settings.batches = run["batches"]
    settings.inactive = run["inactive"]
    settings.seed = run.get("seed", 1)
    settings.ptables = run.get("ptables", True)
    if "free_gas_threshold" in run:
        settings.free_gas_threshold = run["free_gas_threshold"]
    if "resonance_scattering" in run:
        settings.resonance_scattering = run["resonance_scattering"]
    settings.temperature = {"method": "nearest", "tolerance": 1.0}
    settings.source = openmc.IndependentSource(
        space=openmc.stats.Box(
            (-half_lattice,) * 2 + (-1,), (half_lattice,) * 2 + (1,)
        ),
        constraints={"fissionable": True},
    )
    settings.output = {"tallies": False}

    groups = openmc.EnergyFilter(np.sort(np.asarray(spec["groups"])))
    every = openmc.CellFilter([c for _, c in labelled])
    tallies = openmc.Tallies()
    flux = openmc.Tally(name="flux")
    flux.filters = [every, groups]
    flux.scores = ["flux"]
    tallies.append(flux)
    for name in mats:
        cells = [c for label, c in labelled if label.split("/")[0] == name]
        if not cells:
            continue
        t = openmc.Tally(name=f"rates {name}")
        t.filters = [openmc.CellFilter(cells), groups]
        t.nuclides = list(spec["materials"][name]["densities"])
        t.scores = ["absorption", "(n,gamma)", "fission", "nu-fission"]
        tallies.append(t)
    mesh = openmc.RegularMesh()
    mesh.dimension = (n, n)
    mesh.lower_left = (-half_lattice, -half_lattice)
    mesh.upper_right = (half_lattice, half_lattice)
    pins = openmc.Tally(name="pins")
    pins.filters = [openmc.MeshFilter(mesh)]
    pins.scores = ["fission"]
    tallies.append(pins)
    materials = openmc.Materials(mats.values())
    materials.cross_sections = str(LIBRARY / "cross_sections.xml")
    return openmc.Model(geometry, materials, settings, tallies), labelled


def main(spec_path: str, out_path: str) -> None:
    spec = json.loads(Path(spec_path).read_text())
    model, labelled = build(spec)
    labels = list(dict.fromkeys(label for label, _ in labelled))
    n = len(spec["map"])
    with tempfile.TemporaryDirectory() as tmp:
        statepoint = model.run(cwd=tmp, threads=spec["settings"].get("threads"))
        with openmc.StatePoint(statepoint) as sp:
            k = sp.keff
            per_cell = sp.get_tally(name="flux").mean.reshape(len(labelled), -1)[
                :, ::-1
            ]
            flux = {label: np.zeros(per_cell.shape[1]) for label in labels}
            for (label, _), row in zip(labelled, per_cell, strict=True):
                flux[label] += row  # group 0 = highest energy
            rates: dict = {}
            for name in spec["materials"]:
                cells = [label for label, _ in labelled if label.split("/")[0] == name]
                if not cells:
                    continue
                t = sp.get_tally(name=f"rates {name}")
                nuclides = list(spec["materials"][name]["densities"])
                shape = (len(cells), -1, len(nuclides), 4)
                mean = t.mean.reshape(shape)[:, ::-1]
                std = t.std_dev.reshape(shape)[:, ::-1]
                for c, label in enumerate(cells):
                    for j, nuc in enumerate(nuclides):
                        for s, score in enumerate(t.scores):
                            entry = rates.setdefault(label, {}).setdefault(nuc, {})
                            old = entry.get(score)
                            m, v = mean[c, :, j, s], std[c, :, j, s] ** 2
                            if old is not None:  # the same label in several pins
                                m = m + np.asarray(old["mean"])
                                v = v + np.asarray(old["std"]) ** 2
                            entry[score] = {
                                "mean": m.tolist(),
                                "std": np.sqrt(v).tolist(),
                            }
            pins = sp.get_tally(name="pins")
            # mesh bins: x fastest, y from the bottom -> rows from the top
            pin_fission = pins.mean.reshape(n, n)[::-1]
            pin_std = pins.std_dev.reshape(n, n)[::-1]
            out = {
                "name": spec.get("name", ""),
                "k": k.nominal_value,
                "k_std": k.std_dev,
                "regions": labels,
                "flux": [flux[label].tolist() for label in labels],
                "rates": rates,
                "pin_fission": pin_fission.tolist(),
                "pin_fission_std": pin_std.tolist(),
            }
    Path(out_path).write_text(json.dumps(out))
    print(f"k = {out['k']:.5f} +- {out['k_std']:.5f} -> {out_path}")


if __name__ == "__main__":
    main(*sys.argv[1:3])
