"""Multigroup libraries on disk, and a command-line summary.

A library (`groupr.GroupLibrary`) for one nuclide at one temperature takes
15-25 s to build, so it is cached in <data dir>/groups/ as .npz, keyed on the
Doppler table it starts from, the group structure, the weight function, the
Legendre order and the code (fission.processing + endf.distributions).

    uv run python -m fission.processing.library U235 --temps 293.6 900
"""

from __future__ import annotations

import argparse
import hashlib
import json
import time
from pathlib import Path
from typing import Any

import numpy as np

from ..data import EndfLibrary
from ..endf import Evaluation
from ..endf import cache as endf_cache
from ..endf.distributions.sab import scattering_law
from .cache import cached_doppler, code_fingerprint
from .groupr import GroupLibrary, group_library, shape_table
from .groups import XMAS_172, GroupStructure
from .thermal import FreeGas, Tabulated, with_thermal
from .weights import ThermalSlowingDownFission

FORMAT_VERSION = 1
_OPTIONAL = (
    "nu_fission",
    "nu_prompt",
    "nu_delayed",
    "chi_prompt",
    "chi_delayed",
    "chi",
    "decay_constants",
)


def save_library(path: Path, lib: GroupLibrary, meta: dict[str, Any]) -> None:
    payload: dict[str, Any] = {
        "boundaries": np.asarray(lib.groups.boundaries),
        "flux": lib.flux,
        "scatter_elastic": lib.scatter_elastic,
        "scatter_inelastic": lib.scatter_inelastic,
        **{f"sigma{mt}": v for mt, v in lib.sigma.items()},
        **{k: getattr(lib, k) for k in _OPTIONAL if getattr(lib, k) is not None},
        "meta": np.array(json.dumps({
            **meta, "nuclide": lib.nuclide, "temperature": lib.temperature,
            "groups": lib.groups.name, "weight": lib.weight, "lmax": lib.lmax, "notes": lib.notes,
            "thermal": lib.thermal, "thermal_cutoff": lib.thermal_cutoff,
        })),
    }  # fmt: skip
    tmp = path.with_name(path.name + ".tmp.npz")
    np.savez(tmp, **payload)
    tmp.replace(path)


def load_library(path: Path, key: str | None = None) -> GroupLibrary | None:
    try:
        with np.load(path) as data:
            meta = json.loads(str(data["meta"]))
            if key is not None and meta.get("key") != key:
                return None
            groups = GroupStructure(meta["groups"], tuple(data["boundaries"].tolist()))
            sigma = {int(n[5:]): data[n] for n in data.files if n.startswith("sigma")}
            optional = {k: data[k] for k in _OPTIONAL if k in data.files}
            return GroupLibrary(
                meta["nuclide"], meta["temperature"], groups, meta["weight"], meta["lmax"],
                data["flux"], sigma, data["scatter_elastic"], data["scatter_inelastic"],
                notes=list(meta["notes"]), thermal=meta.get("thermal", "target at rest"),
                thermal_cutoff=meta.get("thermal_cutoff", 0.0), **optional,
            )  # fmt: skip
    except OSError, ValueError, KeyError, EOFError:
        return None


def cached_groups(
    nuclide: str,
    temperature: float,
    *,
    groups: GroupStructure = XMAS_172,
    lmax: int = 3,
    thermal: str | None = "free-gas",
    flux_shape: tuple[list[float], list[float]] | None = None,
    shape_name: str = "",
    lib: EndfLibrary | None = None,
    refresh: bool = False,
) -> GroupLibrary:
    """Infinitely dilute multigroup data, weighted with a Maxwellian at the
    library temperature, 1/E and a fission spectrum (`weights`).

    thermal: elastic scattering below 4 eV from a free gas ("free-gas", the
    default), a thermal scattering law of the ENDF thermal sublibrary (e.g.
    "HinH2O" for H1), or None for a target at rest.

    flux_shape = (energies, F): weight w(E) F(E) instead of w(E) inside each
    group (self-shielding, e.g. a narrow-resonance flux); shape_name labels
    the cache file, the key includes the shape itself.
    """
    lib = lib or EndfLibrary()
    weight = ThermalSlowingDownFission.at(temperature)
    meta = {
        "format": FORMAT_VERSION,
        "nuclide": nuclide,
        "temperature": temperature,
        "evaluation_sha256": endf_cache._sha256_file(lib.neutron(nuclide)),
        "boundaries": hashlib.sha256(
            json.dumps(groups.boundaries).encode()
        ).hexdigest(),
        "weight": weight.describe(),
        "lmax": lmax,
        "thermal": thermal,
        "shape": None
        if flux_shape is None
        else hashlib.sha256(
            np.asarray(flux_shape[0]).tobytes() + np.asarray(flux_shape[1]).tobytes()
        ).hexdigest(),
        "code": code_fingerprint(),
    }
    key = hashlib.sha256(json.dumps(meta, sort_keys=True).encode()).hexdigest()
    directory = lib.root / "groups"
    stem = f"{nuclide}-{temperature:g}K-{groups.size}g-{thermal or 'static'}"
    if flux_shape is not None:
        stem += f"-{shape_name or 'shaped'}"
    path = directory / f"{stem}-{key[:16]}.npz"
    if not refresh and path.exists():
        cached = load_library(path, key)
        if cached is not None:
            return cached
    ev = Evaluation.from_file(lib.neutron(nuclide))
    table = cached_doppler(nuclide, temperature, lib=lib)
    if flux_shape is not None:
        table = shape_table(table, *flux_shape)
    gl = group_library(
        ev, table, groups, weight, nuclide=nuclide, temperature=temperature,
        lmax=lmax, flux_shape=flux_shape,
    )  # fmt: skip
    if thermal == "free-gas":
        model: FreeGas | Tabulated = FreeGas(ev.header.awr, temperature)
        gl = with_thermal(gl, model, weight, (table.energies, table.xs[2]))
    elif thermal:
        law = scattering_law(Evaluation.from_file(lib.thermal(thermal)), temperature)
        gl = with_thermal(gl, Tabulated(law), weight, flux_shape=flux_shape)
    directory.mkdir(parents=True, exist_ok=True)
    save_library(path, gl, {**meta, "key": key})
    # every thermal variant of this nuclide/temperature/structure: only entries
    # built by other code are deleted
    endf_cache.prune_stale(
        directory,
        f"{stem}-*.npz"
        if flux_shape is not None
        else f"{nuclide}-{temperature:g}K-{groups.size}g-*.npz",
        keep=path,
        code=code_fingerprint(),
    )
    return gl


def summary(gl: GroupLibrary) -> str:
    g = gl.groups
    lines = [
        f"{gl.nuclide} at {gl.temperature:g} K, {g.name}, weight {gl.weight}, P{gl.lmax}",
        f"thermal scattering: {gl.thermal}"
        + (f" below {gl.thermal_cutoff:g} eV" if gl.thermal_cutoff else ""),
    ]
    head = f"{'g':>4} {'E_low (eV)':>11} {'E_high (eV)':>11} {'total':>9} {'elastic':>9} {'capture':>9}"
    if gl.nu_fission is not None:
        head += f" {'fission':>9} {'nu-fis':>9} {'chi':>8}"
    lines.append(head)
    picks = sorted(
        {0, 10, 20, 30, 45, 60, 75, 87, 100, 120, 140, 164, g.size - 1}
        & set(range(g.size))
    )
    for k in picks:
        lo, hi = g.bounds(k)
        row = f"{k:4d} {lo:11.4g} {hi:11.4g} {gl.sigma[1][k]:9.4g} {gl.sigma[2][k]:9.4g} {gl.sigma.get(102, np.zeros(g.size))[k]:9.4g}"
        if gl.nu_fission is not None and gl.chi is not None:
            row += f" {gl.sigma[18][k]:9.4g} {gl.nu_fission[k]:9.4g} {gl.chi[k]:8.2e}"
        lines.append(row)
    el = gl.scatter_elastic[0].sum(axis=1)
    worst = np.nanmax(np.abs(el[:-3] / gl.sigma[2][:-3] - 1))
    lines.append(
        f"elastic P0 rows reproduce sigma_el to {worst:.1e} (all but the 3 lowest groups)"
    )
    th = g.find(0.0253)
    if th is not None:  # group 0 is the fastest: g' < g is upscattering
        up = np.tril(gl.scatter_elastic[0], -1)[th].sum()
        lines.append(
            f"thermal group {th}: upscatter {up:.3g} b of {el[th]:.3g} b scattering"
        )
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(
        prog="python -m fission.processing.library",
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("nuclides", nargs="+")
    ap.add_argument("--temps", type=float, nargs="+", default=[293.6])
    ap.add_argument("--refresh", action="store_true")
    ap.add_argument(
        "--thermal", default="free-gas", help="free-gas, none, or e.g. HinH2O"
    )
    ap.add_argument("--quiet", action="store_true", help="only report build times")
    args = ap.parse_args(argv)
    lib = EndfLibrary()
    for nuclide in args.nuclides:
        for temp in args.temps:
            t0 = time.perf_counter()
            thermal = None if args.thermal == "none" else args.thermal
            gl = cached_groups(
                nuclide, temp, thermal=thermal, lib=lib, refresh=args.refresh
            )
            print(f"{nuclide} {temp:g} K: {time.perf_counter() - t0:.1f} s")
            if not args.quiet:
                print(summary(gl), "\n")


if __name__ == "__main__":
    main()
