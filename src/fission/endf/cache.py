"""Disk cache for pointwise tables (`pointwise.pointwise`).

Building a full table takes from milliseconds (H-1) to ~90 s (U-238), and
every later step (Doppler broadening, group averaging) starts from these
tables, so they are cached in <data dir>/pointwise/ as .npz files.

The danger with any cache is a stale result. The key is therefore a hash of
everything that determines the table:

* the evaluation file's *contents* (not its name or date),
* the parameters (tol, atol, lo, hi, requested MTs, Rust or Python),
* a code fingerprint: every Python source file of fission.endf plus the
  compiled Rust extension. This is deliberately conservative: any edit, even
  a comment, triggers a rebuild, but a cached table can never silently
  outlive a change to the physics.

When a table built by newer code is written, entries for the same nuclide
built by older code are deleted.

    uv run python -m fission.endf.cache U238 Zr90      # build or load
    uv run python -m fission.endf.cache --list
"""

from __future__ import annotations

import argparse
import functools
import hashlib
import json
import sys
import time
import zipfile
from pathlib import Path
from typing import Any

import numpy as np

from .. import _rust
from ..data import EndfLibrary
from .evaluation import Evaluation
from .pointwise import PointwiseTable, pointwise

FORMAT_VERSION = 1

# what a damaged, truncated or foreign file can raise when loaded
_UNREADABLE = (OSError, ValueError, KeyError, EOFError, zipfile.BadZipFile)


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


@functools.cache
def code_fingerprint() -> str:
    """Hash of the code that produces tables: fission.endf sources + Rust extension."""
    digest = hashlib.sha256()
    package = Path(__file__).parent
    for source in sorted(package.glob("*.py")):
        digest.update(source.name.encode())
        digest.update(source.read_bytes())
    rust_file = getattr(_rust, "__file__", None)
    if rust_file:
        digest.update(Path(rust_file).read_bytes())
    return digest.hexdigest()


def cache_key(endf_file: Path, params: dict) -> tuple[str, dict]:
    """(hex key, metadata) for one evaluation file and parameter set."""
    meta = {
        "format": FORMAT_VERSION,
        "evaluation": endf_file.name,
        "evaluation_sha256": _sha256_file(endf_file),
        "params": params,
        "code": code_fingerprint(),
    }
    key = hashlib.sha256(json.dumps(meta, sort_keys=True).encode()).hexdigest()
    return key, meta


def _save(path: Path, table: PointwiseTable, meta: dict) -> None:
    # entries: "energies", "meta" (JSON string) and one "mt<N>" array per reaction
    payload: dict[str, Any] = {
        "energies": np.asarray(table.energies),
        "meta": np.array(json.dumps(meta)),
        **{f"mt{mt}": np.asarray(values) for mt, values in table.xs.items()},
    }
    tmp = path.with_name(path.name + ".tmp.npz")
    np.savez(tmp, **payload)
    tmp.replace(path)  # atomic: a reader never sees a half-written file


def _read_meta(path: Path) -> dict | None:
    try:
        with np.load(path) as data:
            return json.loads(str(data["meta"]))
    except _UNREADABLE:  # treat as absent
        return None


def _load(path: Path, key: str) -> PointwiseTable | None:
    try:
        with np.load(path) as data:
            meta = json.loads(str(data["meta"]))
            if meta.get("key") != key:
                return None
            xs = {
                int(name[2:]): data[name].tolist()
                for name in data.files
                if name.startswith("mt")
            }
            return PointwiseTable(energies=data["energies"].tolist(), xs=xs)
    except _UNREADABLE:
        return None


def _prune_stale(directory: Path, nuclide: str, keep: Path) -> None:
    """Delete this nuclide's entries that were built by different code."""
    for other in directory.glob(f"{nuclide}-*.npz"):
        if other == keep:
            continue
        meta = _read_meta(other)
        if meta is None or meta.get("code") != code_fingerprint():
            other.unlink(missing_ok=True)


def cached_pointwise(
    nuclide: str,
    *,
    lib: EndfLibrary | None = None,
    cache_dir: Path | None = None,
    tol: float = 1e-3,
    atol: float = 1e-5,
    lo: float | None = None,
    hi: float | None = None,
    mts: list[int] | None = None,
    use_rust: bool = True,
    refresh: bool = False,
) -> PointwiseTable:
    """`pointwise(...)` for a nuclide, loaded from the cache when possible."""
    lib = lib or EndfLibrary()
    endf_file = lib.neutron(nuclide)
    directory = cache_dir or lib.root / "pointwise"
    params = {
        "tol": tol,
        "atol": atol,
        "lo": lo,
        "hi": hi,
        "mts": sorted(mts) if mts else None,
        "use_rust": use_rust,
    }
    key, meta = cache_key(endf_file, params)
    path = directory / f"{nuclide}-{key[:16]}.npz"

    if not refresh and path.exists():
        table = _load(path, key)
        if table is not None:
            return table

    ev = Evaluation.from_file(endf_file)
    table = pointwise(ev, mts=mts, lo=lo, hi=hi, tol=tol, atol=atol, use_rust=use_rust)
    directory.mkdir(parents=True, exist_ok=True)
    _save(path, table, {**meta, "key": key, "points": len(table)})
    _prune_stale(directory, nuclide, keep=path)
    return table


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(
        prog="python -m fission.endf.cache",
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("nuclides", nargs="*", help="e.g. U238 Zr90 H1")
    ap.add_argument("--refresh", action="store_true", help="rebuild even if cached")
    ap.add_argument("--list", action="store_true", help="list cached tables")
    args = ap.parse_args(argv)

    lib = EndfLibrary()
    for nuclide in args.nuclides:
        t0 = time.perf_counter()
        table = cached_pointwise(nuclide, lib=lib, refresh=args.refresh)
        print(
            f"{nuclide:6s} {len(table):9d} points  {len(table.mts):3d} reactions"
            f"  {time.perf_counter() - t0:7.2f} s"
        )
    if args.list:
        directory = lib.root / "pointwise"
        current = code_fingerprint()
        for path in sorted(directory.glob("*.npz")):
            meta = _read_meta(path) or {}
            status = "current" if meta.get("code") == current else "stale"
            size = path.stat().st_size / 1e6
            print(
                f"{path.name:40s} {meta.get('points', '?'):>9} points {size:8.1f} MB  {status}"
                f"  {meta.get('params')}",
                file=sys.stdout,
            )


if __name__ == "__main__":
    main()
