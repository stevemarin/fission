"""Disk cache for Doppler-broadened tables.

Broadening a full U-238 table takes about a minute, so results are kept in
<data dir>/doppler/<nuclide>-<T>K-<key>.npz, like the 0 K tables. The key
hashes everything the result depends on:

* the 0 K table's own key (evaluation contents, its parameters, and the
  fission.endf code that built it),
* the temperature and the broadening parameters (tol, atol, emax),
* a code fingerprint of fission.processing plus the Rust extension.

    uv run python -m fission.processing.cache U238 H1 --temps 293.6 600 900
    uv run python -m fission.processing.cache --list
"""

from __future__ import annotations

import argparse
import functools
import hashlib
import json
import time
from pathlib import Path
from typing import Any

from ..data import EndfLibrary
from ..endf import Evaluation
from ..endf import cache as endf_cache
from ..endf.pointwise import PointwiseTable
from .doppler import broaden_table

FORMAT_VERSION = 1
COLD: dict[str, Any] = {
    "tol": 1e-3,
    "atol": 1e-5,
    "lo": None,
    "hi": None,
    "mts": None,
    "use_rust": True,
}


@functools.cache
def code_fingerprint() -> str:
    """Hash of fission.processing sources and of the code behind the 0 K tables."""
    digest = hashlib.sha256(endf_cache.code_fingerprint().encode())
    for source in sorted(Path(__file__).parent.glob("*.py")):
        digest.update(source.name.encode())
        digest.update(source.read_bytes())
    return digest.hexdigest()


def cached_doppler(
    nuclide: str,
    temperature: float,
    *,
    lib: EndfLibrary | None = None,
    cache_dir: Path | None = None,
    cold_cache_dir: Path | None = None,
    tol: float = 1e-3,
    atol: float = 1e-5,
    emax: float | None = None,
    refresh: bool = False,
) -> PointwiseTable:
    """`broaden_table` of the cached 0 K table, loaded from the cache when possible."""
    lib = lib or EndfLibrary()
    directory = cache_dir or lib.root / "doppler"
    endf_file = lib.neutron(nuclide)
    cold_key, _ = endf_cache.cache_key(endf_file, endf_cache.pointwise_params(**COLD))
    meta = {
        "format": FORMAT_VERSION,
        "nuclide": nuclide,
        "temperature": temperature,
        "cold": cold_key,
        "params": {"tol": tol, "atol": atol, "emax": emax},
        "code": code_fingerprint(),
    }
    key = hashlib.sha256(json.dumps(meta, sort_keys=True).encode()).hexdigest()
    stem = f"{nuclide}-{temperature:g}K"
    path = directory / f"{stem}-{key[:16]}.npz"

    if not refresh and path.exists():
        table = endf_cache.load_table(path, key)
        if table is not None:
            return table

    cold = endf_cache.cached_pointwise(
        nuclide, lib=lib, cache_dir=cold_cache_dir, **COLD
    )
    awr = Evaluation.from_file(endf_file).header.awr
    hot = broaden_table(cold, awr, temperature, tol=tol, atol=atol, emax=emax)
    directory.mkdir(parents=True, exist_ok=True)
    endf_cache.save_table(path, hot, {**meta, "key": key, "points": len(hot)})
    endf_cache.prune_stale(
        directory, f"{stem}-*.npz", keep=path, code=code_fingerprint()
    )
    return hot


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(
        prog="python -m fission.processing.cache",
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("nuclides", nargs="*", help="e.g. U238 H1")
    ap.add_argument(
        "--temps", type=float, nargs="+", default=[293.6, 600.0, 900.0, 1200.0]
    )
    ap.add_argument("--refresh", action="store_true", help="rebuild even if cached")
    ap.add_argument("--list", action="store_true", help="list cached tables")
    args = ap.parse_args(argv)

    lib = EndfLibrary()
    for nuclide in args.nuclides:
        cold = endf_cache.cached_pointwise(nuclide, lib=lib, **COLD)
        for temp in args.temps:
            t0 = time.perf_counter()
            hot = cached_doppler(nuclide, temp, lib=lib, refresh=args.refresh)
            print(
                f"{nuclide:6s} {temp:7.1f} K {len(hot):9d} points"
                f" (0 K: {len(cold):9d})  {time.perf_counter() - t0:7.2f} s"
            )
    if args.list:
        current = code_fingerprint()
        for path in sorted((lib.root / "doppler").glob("*.npz")):
            meta = endf_cache.read_meta(path) or {}
            status = "current" if meta.get("code") == current else "stale"
            size = path.stat().st_size / 1e6
            print(
                f"{path.name:40s} {meta.get('points', '?'):>9} points {size:8.1f} MB  {status}"
            )


if __name__ == "__main__":
    main()
