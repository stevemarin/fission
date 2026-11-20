"""Continuous-energy reference for the VERA pin cell (step C2), via OpenMC.

OpenMC is an independent check, never part of the pipeline: this module only
writes a JSON description of the pin (``spec``), and reads back the JSON that
``checks/openmc_pin.py`` writes from the tools environment (``compare``):

    make openmc-pin STATES=1A       # spec -> OpenMC -> comparison

Both codes start from the same ENDF/B-VIII.1 files (OpenMC's through NJOY,
``checks/openmc_library.py``), the same number densities (by default without
the nuclides this code cannot read yet, ``UNSUPPORTED``) and the same
temperatures and S(alpha, beta). What remains between them is method:
multigroup, self-shielding, flat sources, transport correction. KENO's
reference used ENDF/B-VII.0, so OpenMC minus KENO is mostly data.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from fission.model import Model, load_model, number_densities
from fission.processing.doppler import BOLTZMANN_EV_PER_K
from fission.processing.groups import XMAS_172

from .run import UNSUPPORTED

RUNS = Path("endf-data/openmc/runs")
BANDS = (  # eV: thermal, the slowing-down range of C1, unresolved range, fast
    (1e-5, 0.625),
    (0.625, 0.996),
    (0.996, 24787.5),
    (24787.5, 1.5e5),
    (1.5e5, 2e7),
)


def spec(
    state: str,
    model: Model | str = "vera-p1",
    lattice: str = "cell",
    fuel_rings: int = 3,
    all_nuclides: bool = False,
    particles: int = 100_000,
    batches: int = 250,
    inactive: int = 50,
    ptables: bool = True,
    free_gas_cutoff: float | None = None,
    resonance_scattering: bool = False,
) -> dict:
    """JSON-ready description of a pin lattice (or a single pin cell) at one state."""
    m = load_model(model) if isinstance(model, str) else model
    st = m.states[state]
    lat = m.lattices[lattice]
    if lat.pins is None or lat.pitch is None:
        raise ValueError("the OpenMC check needs a pin lattice")
    used = {p for row in lat.pins for p in row}
    pins = {
        name: {
            "radii": list(m.pins[name].radii),
            "materials": list(m.pins[name].materials),
            "pieces": [
                fuel_rings if mm == "fuel" else 1 for mm in m.pins[name].materials[:-1]
            ],
        }
        for name in sorted(used)
    }
    materials = {}
    for name, mat in m.materials.items():
        dens = number_densities(mat, st)
        if not all_nuclides:
            dens = {n: d for n, d in dens.items() if n not in UNSUPPORTED}
        sab = [
            "c_H_in_H2O"
            for n in dens
            if mat.thermal_law(n) == "HinH2O"  # the only law processed so far
        ]
        materials[name] = {
            "temperature": st.temperature(mat.temperature),
            "densities": dens,
            "sab": sab,
        }
    return {
        "name": f"{m.name} {state}",
        "state": state,
        "pitch": lat.pitch,
        "gap": lat.gap,
        "gap_material": m.pins[lat.pins[0][0]].materials[-1],
        "map": [list(row) for row in lat.pins],
        "pins": pins,
        "materials": materials,
        "groups": list(XMAS_172.boundaries),
        "settings": {
            "particles": particles,
            "batches": batches,
            "inactive": inactive,
            "ptables": ptables,
            # OpenMC samples the target's thermal motion for elastic
            # scattering below free_gas_threshold * kT (default 400); a
            # cutoff in eV (e.g. 4, this code's thermal cutoff) makes it
            # scatter on targets at rest above, as the slowing-down march
            **(
                {}
                if free_gas_cutoff is None
                else {
                    "free_gas_threshold": free_gas_cutoff
                    / (
                        BOLTZMANN_EV_PER_K
                        * st.temperature(m.materials["fuel"].temperature)
                    )
                }
            ),
            # exact target motion in U-238's resonances: relative velocity
            # sampling against the 0 K elastic cross section (OpenMC 'rvs')
            **(
                {
                    "resonance_scattering": {
                        "enable": True,
                        "method": "rvs",
                        "energy_min": 1.0,
                        "energy_max": 1000.0,
                        "nuclides": ["U238"],
                    }
                }
                if resonance_scattering
                else {}
            ),
        },
    }


# ---------------------------------------------------------------- comparison


def reference_rates(result: dict) -> dict[tuple[str, str, str], np.ndarray]:
    """(region, nuclide, score) -> group rates per fission neutron produced."""
    production = sum(
        np.sum(scores["nu-fission"]["mean"])
        for nucs in result["rates"].values()
        for scores in nucs.values()
    )
    return {
        (region, n, score): np.asarray(v["mean"]) / production
        for region, nucs in result["rates"].items()
        for n, scores in nucs.items()
        for score, v in scores.items()
    }


def our_rates(r) -> dict[tuple[str, str, str], np.ndarray]:
    """The same table from a ``heterogeneous.Result``: regions are its MOC
    materials (``fuel/0``... when the pellet is split into rings)."""
    sol = r.solution
    materials = np.array(sol.problem.geometry.materials())
    production = sol.rates("nu_fission").sum()
    out = {}
    for region, libs in r.libraries.items():
        mask = materials == region
        flux = (sol.problem.areas[mask, None] * sol.flux[mask]).sum(axis=0)
        for n, lib in libs.items():
            d = r.densities[region][n] * flux / production
            out[(region, n, "(n,gamma)")] = d * lib.sigma.get(102, 0 * flux)
            out[(region, n, "fission")] = d * lib.sigma.get(18, 0 * flux)
    return out


def _merge(rates: dict, keep_rings: bool) -> dict:
    """Sum region labels ``fuel/0``, ``fuel/1`` into ``fuel`` unless kept."""
    out: dict = {}
    for (region, n, s), v in rates.items():
        base = region.split("/")[0]
        key = (region if keep_rings and base == "fuel" else base, n, s)
        out[key] = out.get(key, 0) + v
    return out


def compare(r, result: dict) -> str:
    """Report: k, then capture and fission by nuclide, ring and energy band."""
    ours_raw = our_rates(r)
    ref_raw = reference_rates(result)

    def fuel_labels(rates: dict) -> set[str]:
        return {region for region, _, _ in rates if region.startswith("fuel")}

    # ring by ring only if both have the same rings (OpenMC tallies 3)
    split = any("/" in region for region, _, _ in ours_raw)
    rings = split and fuel_labels(ours_raw) == fuel_labels(ref_raw)
    ref = _merge(ref_raw, rings)
    ours = _merge(ours_raw, rings)
    k, k_std = result["k"], result["k_std"]
    lines = [
        (
            f"{result.get('name', '')}: this code k = {r.solution.k:.5f}, "
            f"OpenMC k = {k:.5f} +- {k_std:.5f}: "
            f"{1e5 * (r.solution.k - k):+.0f} +- {1e5 * k_std:.0f} pcm"
        ),
        "  per fission neutron      this code    OpenMC   difference",
    ]
    big = sorted(
        (key for key in ours if key in ref),
        key=lambda key: -float(np.sum(ref[key])),
    )
    for key in big[:12]:
        a, b = float(np.sum(ours[key])), float(np.sum(ref[key]))
        region, n, s = key
        label = f"{n} {s} ({region})"
        lines.append(f"  {label:26s} {a:9.5f} {b:9.5f} {a - b:+9.5f}")
    e = np.asarray(XMAS_172.boundaries)
    mid = np.sqrt(e[:-1] * e[1:])
    lines.append("  by energy band (sum over fuel regions):")
    for n, s in (("U238", "(n,gamma)"), ("U235", "fission"), ("U235", "(n,gamma)")):
        cols = []
        for lo, hi in BANDS:
            band = (mid > lo) & (mid < hi)
            a = sum(float(np.sum(v[band])) for (rg, nn, ss), v in ours.items()
                    if rg.startswith("fuel") and nn == n and ss == s)  # fmt: skip
            b = sum(float(np.sum(v[band])) for (rg, nn, ss), v in ref.items()
                    if rg.startswith("fuel") and nn == n and ss == s)  # fmt: skip
            cols.append(f"{lo:g}-{hi:g} eV: {a - b:+.5f}")
        lines.append(f"    {n} {s}: " + ", ".join(cols))
    if len(result.get("pin_fission", [[0]])) > 1:
        from .assembly import compare_pins, pin_powers

        omc = np.asarray(result["pin_fission"])
        ours_pins = pin_powers(r)
        omc = np.where(np.isnan(ours_pins), np.nan, omc)
        c = compare_pins(ours_pins, omc / np.nanmean(omc))
        lines.append(f"  pin powers vs OpenMC: RMS {c.rms:.2f}%, max {c.max_abs:.2f}% "
                     f"at {c.worst}; max {c.ours_max:.4f} ({c.ref_max:.4f}), "
                     f"min {c.ours_min:.4f} ({c.ref_min:.4f})")  # fmt: skip
    return "\n".join(lines)


def few_group_reference(result: dict, cut: float = 0.625) -> dict[str, np.ndarray]:
    """Cell-homogenized two-group constants from OpenMC's tallies, condensed
    with OpenMC's own (infinite-lattice) spectrum: flux ratio and
    Sigma_x,G = sum of reaction rates / sum of track-length flux."""
    e = np.asarray(XMAS_172.boundaries)
    mid = np.sqrt(e[:-1] * e[1:])
    owner = (mid < cut).astype(int)  # 0 fast, 1 thermal
    flux = np.asarray(result["flux"]).sum(axis=0)  # all regions
    out = {"flux": np.bincount(owner, flux, 2)}
    for score, key in (("absorption", "sigma_a"), ("nu-fission", "nu_sigma_f"),
                       ("fission", "sigma_f")):  # fmt: skip
        rate = sum(
            np.asarray(v[score]["mean"])
            for nucs in result["rates"].values()
            for v in nucs.values()
        )
        out[key] = np.bincount(owner, rate, 2) / out["flux"]
    return out


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(prog="python -m fission.pincell.reference")
    sub = ap.add_subparsers(dest="command", required=True)
    w = sub.add_parser("spec", help="write the JSON spec for checks/openmc_pin.py")
    w.add_argument("state")
    w.add_argument("--model", default="vera-p1")
    w.add_argument("--out", type=Path)
    w.add_argument("--all-nuclides", action="store_true")
    w.add_argument("--no-ptables", action="store_true")
    w.add_argument("--particles", type=int, default=100_000)
    w.add_argument("--batches", type=int, default=250)
    w.add_argument("--tag", default="", help="suffix of the file names")
    w.add_argument("--lattice", default=None, help="default: the model's only lattice")
    w.add_argument("--free-gas-cutoff", type=float, default=None,
                   help="eV: target at rest above (OpenMC default: 400 kT)")  # fmt: skip
    w.add_argument("--resonance-scattering", action="store_true",
                   help="exact U-238 resonance scattering kernel (OpenMC rvs)")  # fmt: skip
    c = sub.add_parser("compare", help="run this code and compare with OpenMC")
    c.add_argument("state")
    c.add_argument("--model", default="vera-p1")
    c.add_argument("--result", type=Path)
    c.add_argument("--shielding", default="sd")
    c.add_argument("--rings", action="store_true", help="shield each pellet ring")
    c.add_argument("--sph", action="store_true", help="SPH factors")
    c.add_argument("--clad", action="store_true", help="clad slows down with the rings")
    c.add_argument("--moving", action="store_true",
                  help="U-238 scatters on moving targets (sd)")  # fmt: skip
    c.add_argument("--lattice", default=None)
    c.add_argument("--azim", type=int, default=None, help="default: heterogeneous.AZIM")
    c.add_argument("--spacing", type=float, default=None)
    c.add_argument("--fuel-rings", type=int, default=None)
    args = ap.parse_args(argv)
    if args.command == "spec":
        model = load_model(args.model)
        s = spec(args.state, model, args.lattice or next(iter(model.lattices)),
                 all_nuclides=args.all_nuclides,
                 particles=args.particles, batches=args.batches,
                 ptables=not args.no_ptables,
                 free_gas_cutoff=args.free_gas_cutoff,
                 resonance_scattering=args.resonance_scattering)  # fmt: skip
        out = args.out or RUNS / f"{args.state}{args.tag}.spec.json"
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(json.dumps(s, indent=1))
        print(out)
    else:
        from .heterogeneous import AZIM, FUEL_RINGS, SPACING, run

        result = json.loads(
            (args.result or RUNS / f"{args.state}.result.json").read_text()
        )
        model = load_model(args.model)
        r = run(args.state, args.shielding, model, args.lattice or next(iter(model.lattices)),
                args.azim or AZIM, args.spacing or SPACING,
                fuel_rings=args.fuel_rings or FUEL_RINGS, ring_shielding=args.rings,
                sph=args.sph, coupled_clad=args.clad, moving_targets=args.moving)  # fmt: skip
        print(compare(r, result))


if __name__ == "__main__":
    main()
