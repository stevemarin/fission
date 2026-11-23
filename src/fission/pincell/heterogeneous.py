"""k-infinity of the heterogeneous VERA pin cell: MOC on ENDF/B-VIII.1 libraries.

    uv run python -m fission.pincell.heterogeneous 1A --shielding equivalence

Each region keeps its own material (fuel, gap, clad, water), with group
constants from chapters 10-11 at the region's temperature. Shielding:

* none: infinitely dilute libraries;
* equivalence: fuel and clad nuclides averaged with the NR flux of their
  equivalent homogeneous mixture (``shielding.py``: Dancoff correction from
  MOC, escape cross section, background sigma_0 per nuclide);
* cp: the same slowing-down balance without the rational approximation:
  the escape probability P(Sigma) of each resonant region is computed by
  MOC and the region's flux shape follows from its pointwise total cross
  section (one shape for all nuclides of the region), narrow resonances;
* sd: as cp, but the fuel's flux between 1 eV and 25 keV comes from solving
  the slowing-down equation on the ultra-fine pointwise grid, with the exact
  elastic kernel of every fuel nuclide (no NR/WR assumption for the fuel).

MOC scatters isotropically, so anisotropic scattering is folded into a
**transport correction** (outflow approximation): Sigma_tr,g = Sigma_t,g -
sum_g' Sigma_s1,g->g', and the same amount is removed from the within-group
P0 scattering, which leaves absorption unchanged.
"""

from __future__ import annotations

import argparse
import time
from dataclasses import dataclass, field, replace

import numpy as np

from fission.model import MacroXS, Model, State, load_model, number_densities
from fission.processing.groupr import GroupLibrary
from fission.processing.groups import XMAS_172
from fission.processing.library import cached_groups
from fission.transport import Geometry, cyclic_tracks
from fission.transport.moc import Solution, solve, tracking, with_materials

from .run import UNSUPPORTED
from .shielding import (
    Escape,
    collision_matrix,
    coupled_collisions,
    coupled_flux,
    dancoff,
    escape_probability,
    potential_scattering,
    region_flux,
    ring_escape,
    shielded,
    sigma0,
    sph_factors,
    with_sph,
)

RESONANT = ("fuel", "clad")
MOVING = ("U238",)  # nuclides that scatter on moving targets with moving_targets
# Default discretization, converged on VERA 1D to about 20 pcm (docs ch. 19):
# 256 azimuthal directions over 2 pi, 0.02 cm tracks, 3 Tabuchi-Yamamoto
# polar angles (the default of `tracking`), 10 equal-area pellet rings.
# The old 64 / 0.02 / 3 Gauss-Legendre / 3 rings was about 270 pcm low.
AZIM = 256
SPACING = 0.02
POLAR = 3
FUEL_RINGS = 10


def macroscopic(name: str, items: list[tuple[float, GroupLibrary]]) -> MacroXS:
    """Transport-corrected macroscopic constants of a material (1/cm)."""
    groups = items[0][1].groups.size
    total = np.zeros(groups)
    scatter = np.zeros((2, groups, groups))
    fission = np.zeros(groups)
    nu_fission = np.zeros(groups)
    chi = np.zeros(groups)
    for density, lib in items:
        total += density * lib.sigma[1]
        s = lib.scatter
        scatter[0] += density * s[0]
        if s.shape[0] > 1:
            scatter[1] += density * s[1]
        fission += density * lib.sigma.get(18, np.zeros(groups))
        if lib.nu_fission is not None and lib.chi is not None:
            nu_fission += density * lib.nu_fission
            # chi weighted by this nuclide's fission neutrons in the library's
            # weighting spectrum
            chi += density * float(lib.nu_fission @ lib.flux) * lib.chi
    correction = scatter[1].sum(axis=1)  # sum_g' Sigma_s1,g->g'
    transport = total - correction
    s0 = scatter[0] - np.diag(correction)
    if chi.sum() > 0:
        chi /= chi.sum()
    return MacroXS(
        name, transport, transport - s0.sum(axis=1), s0, fission, nu_fission, chi
    )


@dataclass
class Result:
    state: State
    shielding: str
    solution: Solution
    escapes: dict[str, Escape]
    backgrounds: dict[str, dict[str, float]]  # material -> nuclide -> sigma_0
    densities: dict[str, dict[str, float]]
    libraries: dict[str, dict[str, GroupLibrary]]
    sph: dict[str, np.ndarray] = field(default_factory=dict)  # region -> mu_g

    def absorption(self) -> dict[tuple[str, str], float]:
        """Absorption (sigma_t - out-scattering, P0) by (material, nuclide),
        per fission neutron produced."""
        sol = self.solution
        g = sol.problem.geometry
        materials = np.array(g.materials())
        production = sol.rates("nu_fission").sum()
        out = {}
        for m, libs in self.libraries.items():
            mask = materials == m
            flux = (sol.problem.areas[mask, None] * sol.flux[mask]).sum(axis=0)
            for n, lib in libs.items():
                removal = lib.sigma[1] - lib.scatter[0].sum(axis=1)
                out[(m, n)] = self.densities[m][n] * float(removal @ flux) / production
        return out


def run(
    state: str,
    shielding: str = "sd",
    model: Model | str = "vera-p1",
    lattice: str = "cell",
    azim: int = AZIM,
    spacing: float = SPACING,
    polar: int = POLAR,
    fuel_rings: int = FUEL_RINGS,
    sectors: int = 8,
    outer_ring: float | None = 0.6,
    resonant: tuple[str, ...] = RESONANT,
    ring_shielding: bool = False,
    sph: bool = False,
    coupled_clad: bool = False,
    moving_targets: bool = False,
    verbose: bool = False,
    cmfd: bool = True,
) -> Result:
    """k-infinity of a model pin cell with MOC.

    ring_shielding (cp and sd only): every fuel ring gets its own flux shape
    (``ring_escape``) and its own libraries; MOC sees the rings as the
    materials ``fuel/0`` (centre) ... ``fuel/{fuel_rings - 1}`` (rim).
    sph (cp and sd only): superhomogenization factors (``sph_factors``) for
    every shielded region, so that the multigroup flux of the self-shielding
    model reproduces its group-averaged pointwise flux.
    coupled_clad (sd with ring_shielding): the pellet rings and the clad
    slow down together (``coupled_flux``), instead of the clad being part of
    the narrow-resonance outside of the fuel and getting its own NR shape.
    moving_targets (sd): U-238 scatters on nuclei in thermal motion at the
    fuel temperature below 1 keV (``moving.free_gas_kernel``), instead of
    on nuclei at rest.
    """
    m = load_model(model) if isinstance(model, str) else model
    st = m.states[state]
    densities = {
        name: {
            n: d for n, d in number_densities(mat, st).items() if n not in UNSUPPORTED
        }
        for name, mat in m.materials.items()
    }
    geometry = Geometry.from_lattice(
        m, lattice, {"fuel": fuel_rings}, sectors, outer_ring
    )
    tracks = cyclic_tracks(geometry.width, geometry.height, azim, spacing)
    base = tracking(geometry, tracks, polar)
    escapes: dict[str, Escape] = {}
    backgrounds: dict[str, dict[str, float]] = {}
    shapes: dict[str, tuple[list[float], list[float]]] = {}  # by MOC material
    # for SPH, one entry per self-shielding model: (pointwise energy grid,
    # psi_M, label -> (pointwise flux, collision source))
    sph_inputs: list[tuple[np.ndarray, np.ndarray, dict]] = []
    if shielding not in ("none", "equivalence", "cp", "sd"):
        raise ValueError(f"unknown shielding {shielding!r}")
    split = ring_shielding and shielding in ("cp", "sd")
    if ring_shielding and not split:
        raise ValueError("ring shielding needs shielding cp or sd")
    if sph and shielding not in ("cp", "sd"):
        raise ValueError("SPH needs shielding cp or sd")
    if coupled_clad and not (split and shielding == "sd"):
        raise ValueError("coupled_clad needs shielding sd with ring_shielding")
    if moving_targets and shielding != "sd":
        raise ValueError("moving targets need shielding sd")
    moving = MOVING if moving_targets else ()
    # materials the model marks shielded (absorber rods) join the fuel and clad
    present = set(geometry.materials())
    resonant = tuple(resonant) + tuple(
        k for k, mat in m.materials.items()
        if mat.shielded and k in present and k not in resonant
    )  # fmt: skip
    if shielding != "none":
        sigma_p = {
            name: sum(d * potential_scattering(n) for n, d in dens.items())
            for name, dens in densities.items()
        }
        areas = dict.fromkeys(m.materials, 0.0)
        for region_material, area in zip(geometry.materials(), base.areas, strict=True):
            areas[region_material] += area
        cell = [
            (areas[k], densities[k], st.temperature(m.materials[k].temperature))
            for k in m.materials
        ]
        if coupled_clad:
            coupled = coupled_flux(
                (densities["fuel"], st.temperature(m.materials["fuel"].temperature)),
                (densities["clad"], st.temperature(m.materials["clad"].temperature)),
                coupled_collisions(base, sigma_p), cell,
                (list(m.materials).index("fuel"), list(m.materials).index("clad")),
                moving=moving,
            )  # fmt: skip
            for label, phi in coupled.flux.items():
                shapes[label] = (coupled.energy.tolist(), phi.tolist())
            sph_inputs.append((coupled.energy, coupled.level, {
                x: (coupled.flux[x], coupled.source[x]) for x in coupled.flux
            }))  # fmt: skip
        for name in resonant:
            escapes[name] = dancoff(base, sigma_p, name)
            backgrounds[name] = sigma0(densities[name], escapes[name].sigma_e)
            if shielding in ("cp", "sd") and not (
                coupled_clad and name in ("fuel", "clad")
            ):
                table = escape_probability(base, sigma_p, name)
                temperature = st.temperature(m.materials[name].temperature)
                marched = shielding == "sd" and name == "fuel"
                ring_table, collisions = None, None
                if split and name == "fuel" and marched:
                    _, matrix, volumes = collision_matrix(base, sigma_p, name)
                    collisions = (matrix, volumes)
                    ring_table = (1.0 - matrix.sum(axis=2)).T  # E_j
                elif split and name == "fuel":
                    ring_table = ring_escape(base, sigma_p, name)[1]
                details: dict[str, np.ndarray] = {}
                shape = region_flux(
                    densities[name], temperature, table, slowing_down=marched,
                    cell=cell, region=list(m.materials).index(name),
                    rings=ring_table, details=details, collisions=collisions,
                    moving=moving if name == "fuel" else (),
                )  # fmt: skip
                references = {}
                if ring_table is None:
                    shapes[name] = shape  # type: ignore[assignment]
                    references[name] = (
                        details["flux"],
                        details["source"] * details["total"],
                    )
                else:
                    for k, ring_shape in enumerate(shape):
                        shapes[f"{name}/{k}"] = ring_shape  # type: ignore[assignment]
                        own = details.get(f"source/{k}", details["source"])
                        references[f"{name}/{k}"] = (details[f"flux/{k}"],
                                                     own * details["total"])  # fmt: skip
                sph_inputs.append((details["energy"], details["level"], references))
    full = base  # the tracking with model materials (SPH works on it)
    if split:
        base = replace(base, geometry=geometry.split_rings("fuel"))
    libraries: dict[str, dict[str, GroupLibrary]] = {}
    for label in sorted(set(base.geometry.materials())):
        name = label.split("/")[0]
        mat = m.materials[name]
        temperature = st.temperature(mat.temperature)
        libraries[label] = {}
        for n in densities[name]:
            law = mat.thermal_law(n)
            if label in shapes:
                tag = (
                    f"{shielding}{'c' if coupled_clad else ''}{'m' if moving else ''}"
                    f"-{label.replace('/', '-')}"
                )
                lib = cached_groups(n, temperature, thermal=law, flux_shape=shapes[label],
                                    shape_name=f"{tag}-{st.name}")  # fmt: skip
            elif name in backgrounds:
                lib = shielded(n, temperature, backgrounds[name][n], law)
            else:
                lib = cached_groups(n, temperature, thermal=law)
            libraries[label][n] = lib
    factors: dict[str, np.ndarray] = {}
    names = np.array(base.geometry.materials())
    for energy, level, references in sph_inputs if sph else ():
        sigma_t = {}
        for x in references:
            dens = densities[x.split("/")[0]]
            sigma_t[x] = np.zeros(XMAS_172.size)
            for n, lib in libraries[x].items():
                sigma_t[x] += dens[n] * lib.sigma[1]
        mu = sph_factors(
            full, sigma_p, names, energy, level, references, sigma_t, XMAS_172
        )
        factors.update(mu)
        for x in references:
            libraries[x] = {n: with_sph(lib, mu[x]) for n, lib in libraries[x].items()}
    region_densities = {label: densities[label.split("/")[0]] for label in libraries}
    xs = {
        label: macroscopic(
            label, [(region_densities[label][n], lib) for n, lib in libs.items()]
        )
        for label, libs in libraries.items()
    }
    sol = solve(with_materials(base, xs), tol_k=1e-6, tol_source=1e-5, verbose=verbose,
                cmfd=cmfd)  # fmt: skip
    return Result(st, shielding, sol, escapes, backgrounds, region_densities, libraries,
                  factors)  # fmt: skip


def report(r: Result) -> str:
    st, sol = r.state, r.solution
    ref = st.reference
    versus = (
        f"   (reference {ref.keff:.5f}: {1e5 * (sol.k - ref.keff):+.0f} pcm)"
        if ref
        else ""
    )
    lines = [
        (
            f"VERA {st.name}, heterogeneous pin cell (MOC), shielding = {r.shielding}: "
            f"fuel {st.fuel_temperature:g} K, moderator {st.moderator_temperature:g} K"
        ),
        f"  k-infinity = {sol.k:.5f}{versus}",
        (
            f"  {sol.problem.geometry.regions} regions, {sol.iterations} iterations, "
            f"balance {sol.balance():.5f}"
        ),
    ]
    for name, e in r.escapes.items():
        lines.append(f"  {name}: Dancoff C = {e.dancoff:.4f}, mean chord {e.chord:.4f} cm, "
                     f"Sigma_e = {e.sigma_e:.4f} /cm")  # fmt: skip
    if "fuel" in r.backgrounds and "U238" in r.backgrounds["fuel"]:
        lines.append(f"  sigma_0: U-238 {r.backgrounds['fuel']['U238']:.1f} b, "
                     f"U-235 {r.backgrounds['fuel'].get('U235', 0):.0f} b")  # fmt: skip
    lines.append("  absorption per fission neutron (removal for (n,2n)):")
    rates = sorted(r.absorption().items(), key=lambda kv: -kv[1])[:8]
    for (m, n), a in rates:
        lines.append(f"    {n + '(' + m + ')':18s} {a:8.4f}")
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(prog="python -m fission.pincell.heterogeneous",
                                 description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)  # fmt: skip
    ap.add_argument("states", nargs="+")
    ap.add_argument(
        "--shielding", choices=["none", "equivalence", "cp", "sd"], default="sd"
    )
    ap.add_argument("--model", default="vera-p1")
    ap.add_argument("--azim", type=int, default=AZIM)
    ap.add_argument("--spacing", type=float, default=SPACING)
    ap.add_argument("--fuel-rings", type=int, default=FUEL_RINGS,
                    help="equal-area pellet rings")  # fmt: skip
    ap.add_argument(
        "--rings", action="store_true", help="shield each fuel ring separately"
    )
    ap.add_argument("--sph", action="store_true", help="SPH factors (cp, sd)")
    ap.add_argument("--cp", action="store_true",
                    help="also solve the cell by collision probabilities (cross-check)")  # fmt: skip
    ap.add_argument("--clad", action="store_true",
                    help="slow down the clad with the pellet rings (sd --rings)")  # fmt: skip
    ap.add_argument("--moving", action="store_true",
                   help="U-238 scatters on moving targets (sd)")  # fmt: skip
    args = ap.parse_args(argv)
    model = load_model(args.model)
    for name in args.states:
        t0 = time.perf_counter()
        r = run(name, args.shielding, model, azim=args.azim, spacing=args.spacing,
                fuel_rings=args.fuel_rings,
                ring_shielding=args.rings, sph=args.sph,
                coupled_clad=args.clad, moving_targets=args.moving)  # fmt: skip
        print(report(r))
        if args.cp:
            from fission.transport.cp import pin_cell

            k_cp, _, _ = pin_cell(r.solution.problem)
            print(f"  CP (Wigner-Seitz cell, white boundary), same cross sections: "
                  f"k = {k_cp:.5f} ({1e5 * (k_cp - r.solution.k):+.0f} pcm vs MOC)")  # fmt: skip
        print(f"  ({time.perf_counter() - t0:.0f} s)\n", flush=True)


if __name__ == "__main__":
    main()
