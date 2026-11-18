"""Fetch ENDF/B-VIII.1 evaluations from the NNDC on demand.

The NNDC distributes each ENDF/B-VIII.1 sublibrary as one tarball and lists its
size and checksums in a release ``metadata.json``. Individual evaluations are not
published separately, so "on demand" works in two levels:

1. the first time a sublibrary is needed, its tarball is downloaded once and
   verified against the published SHA-1;
2. each evaluation (one ENDF-6 text file per material) is extracted from the
   tarball only when it is asked for.

Everything lives under one data directory (``ENDF_DATA_DIR``, default
``<project>/endf-data``)::

    endf-data/
        metadata.json                               release metadata (cached)
        archives/neutrons-version.VIII.1.tar.gz     sublibrary tarballs
        neutrons/index.txt                          member list of the tarball
        neutrons/n-092_U_238.endf                   extracted evaluations
        thermal_scatt/tsl-HinH2O.endf

Typical use::

    lib = EndfLibrary()
    u238 = lib.neutron("U238")          # -> Path to n-092_U_238.endf
    paths = lib.neutrons(["U235", "U238", "O16", "H1"])   # one pass over the tarball
    h2o = lib.thermal("HinH2O")         # S(alpha, beta) for hydrogen bound in water
    lib.isotopes("Zr")                  # ['Zr90', 'Zr91', ...]; there are no
                                        # natural-element evaluations in VIII.1

or from the shell::

    uv run python -m nodal.data U235 U238 O16 H1 --thermal HinH2O
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import shutil
import sys
import tarfile
import urllib.request
from pathlib import Path

RELEASE = "B-VIII.1"
BASE_URL = f"https://www.nndc.bnl.gov/endf-releases/releases/{RELEASE}"
USER_AGENT = "nodal-endf-fetch/0.1"
CHUNK = 1 << 20  # 1 MiB


_NUCLIDE = re.compile(r"^([A-Za-z]{1,2})-?(\d{1,3})(m\d?)?$")
_NEUTRON_FILE = re.compile(r"^n-(\d{3})_([A-Za-z]{1,2})_(\d{3})(m\d)?\.endf$")


def _neutron_file_stem(name: str) -> tuple[str, int, str]:
    """'U238' / 'U-238' / 'am242m' -> ('U', 238, 'm1') style key (symbol, A, state)."""
    m = _NUCLIDE.match(name.strip())
    if not m:
        raise ValueError(
            f"cannot parse nuclide name {name!r} (try e.g. 'U238', 'H1', 'Am242m')"
        )
    symbol, mass, state = m.groups()
    if state == "m":
        state = "m1"
    return symbol.capitalize(), int(mass), state or ""


def default_data_dir() -> Path:
    env = os.environ.get("ENDF_DATA_DIR")
    if env:
        return Path(env).expanduser()
    return Path(__file__).resolve().parents[2] / "endf-data"


def _download(url: str, dest: Path, sha1: str | None = None) -> None:
    """Stream url to dest via a .part file; verify SHA-1 before the final rename."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    part = dest.with_name(dest.name + ".part")
    digest = hashlib.sha1()
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(request) as response, open(part, "wb") as out:
        total = int(response.headers.get("Content-Length") or 0)
        done, last_pct = 0, -1
        while chunk := response.read(CHUNK):
            out.write(chunk)
            digest.update(chunk)
            done += len(chunk)
            pct = 100 * done // total if total else -1
            if total > 10 * CHUNK and pct != last_pct and pct % 5 == 0:
                print(
                    f"\r  {dest.name}: {done / 1e6:7.1f} / {total / 1e6:.1f} MB ({pct}%)",
                    end="",
                    file=sys.stderr,
                    flush=True,
                )
                last_pct = pct
    if total > 10 * CHUNK:
        print(file=sys.stderr)
    if sha1 is not None and digest.hexdigest() != sha1.lower():
        part.unlink()
        raise OSError(
            f"SHA-1 mismatch for {url}: got {digest.hexdigest()}, expected {sha1}"
        )
    part.replace(dest)


class EndfLibrary:
    """Local, lazily populated copy of the ENDF/B-VIII.1 library (ENDF-6 format)."""

    def __init__(self, root: str | Path | None = None):
        self.root = Path(root) if root is not None else default_data_dir()

    # -------------------------------------------------------------- release info
    def metadata(self) -> dict:
        """Release metadata: sublibrary tarball names, sizes and checksums."""
        path = self.root / "metadata.json"
        if not path.exists():
            _download(f"{BASE_URL}/metadata.json", path)
        return json.loads(path.read_text())

    def sublibraries(self) -> list[str]:
        """Sublibrary keys, e.g. 'neutrons', 'thermal_scatt', 'decay', 'nfy'."""
        return sorted(self.metadata()["sublibraries"])

    # ------------------------------------------------------------------ archives
    def archive(self, sublibrary: str) -> Path:
        """Path to the sublibrary tarball, downloading and verifying it if needed."""
        subs = self.metadata()["sublibraries"]
        if sublibrary not in subs:
            raise KeyError(
                f"unknown sublibrary {sublibrary!r}; choose from {sorted(subs)}"
            )
        info = subs[sublibrary]["files"]["ENDF-6"]
        path = self.root / "archives" / info["name"]
        if not path.exists():
            print(
                f"Downloading {sublibrary} sublibrary ({info['size']}) from NNDC ...",
                file=sys.stderr,
            )
            _download(
                f"{BASE_URL}/{sublibrary}/{info['name']}",
                path,
                sha1=info["checksums"]["SHA1"],
            )
        return path

    def members(self, sublibrary: str) -> list[str]:
        """File names inside the sublibrary tarball (cached in <sublibrary>/index.txt)."""
        index = self.root / sublibrary / "index.txt"
        if not index.exists():
            with tarfile.open(self.archive(sublibrary), "r:gz") as tar:
                names = sorted(m.name for m in tar.getmembers() if m.isfile())
            index.parent.mkdir(parents=True, exist_ok=True)
            index.write_text("\n".join(names) + "\n")
        return index.read_text().split()

    def extract(self, sublibrary: str, members: list[str]) -> list[Path]:
        """Extract the given tarball members (one pass) and return their local paths."""
        dest_dir = self.root / sublibrary
        paths = {m: dest_dir / Path(m).name for m in members}
        missing = {m for m, p in paths.items() if not p.exists()}
        if missing:
            with tarfile.open(self.archive(sublibrary), "r:gz") as tar:
                for member in tar:
                    if member.name in missing:
                        src = tar.extractfile(member)
                        assert src is not None
                        tmp = paths[member.name].with_suffix(".part")
                        with src, open(tmp, "wb") as out:
                            shutil.copyfileobj(src, out)
                        tmp.replace(paths[member.name])
                        missing.discard(member.name)
                        if not missing:
                            break
        if missing:
            raise FileNotFoundError(f"not in {sublibrary} archive: {sorted(missing)}")
        return [paths[m] for m in members]

    # ------------------------------------------------------- friendly lookups
    def _neutron_table(self) -> dict[tuple[str, int, str], str]:
        table = {}
        for member in self.members("neutrons"):
            m = _NEUTRON_FILE.match(Path(member).name)
            if m:
                _, symbol, mass, state = m.groups()
                table[(symbol.capitalize(), int(mass), state or "")] = member
        return table

    def isotopes(self, element: str) -> list[str]:
        """Nuclides of an element with neutron evaluations, e.g. 'Zr' -> ['Zr90', ...]."""
        symbol = element.strip().capitalize()
        return [
            f"{s}{a}{st[:1]}"
            for (s, a, st) in sorted(self._neutron_table())
            if s == symbol
        ]

    def neutrons(self, names: list[str]) -> list[Path]:
        """Incident-neutron evaluations for several nuclides (single tarball pass)."""
        table = self._neutron_table()
        members = []
        for name in names:
            key = _neutron_file_stem(name)
            if key not in table:
                hint = self.isotopes(key[0])
                raise KeyError(
                    f"no ENDF/B-VIII.1 neutron evaluation for {name!r}"
                    + (f"; available: {hint}" if hint else "")
                )
            members.append(table[key])
        return self.extract("neutrons", members)

    def neutron(self, name: str) -> Path:
        """Incident-neutron evaluation (ENDF-6 file) for one nuclide, e.g. 'U238'."""
        return self.neutrons([name])[0]

    def thermal_materials(self) -> list[str]:
        """Thermal scattering law materials, e.g. 'HinH2O', 'OinUO2', 'crystalline-graphite'."""
        return sorted(
            Path(m).name[4:-5]
            for m in self.members("thermal_scatt")
            if Path(m).name.startswith("tsl-") and m.endswith(".endf")
        )

    def thermal(self, name: str) -> Path:
        """Thermal scattering law (S(alpha, beta), ENDF MF7) evaluation, e.g. 'HinH2O'."""
        for member in self.members("thermal_scatt"):
            if Path(member).name.lower() == f"tsl-{name}.endf".lower():
                return self.extract("thermal_scatt", [member])[0]
        raise KeyError(
            f"no thermal scattering evaluation {name!r}; "
            f"see EndfLibrary().thermal_materials()"
        )


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(
        prog="python -m nodal.data",
        description="Fetch ENDF/B-VIII.1 evaluations from the NNDC.",
    )
    ap.add_argument(
        "nuclides", nargs="*", help="incident-neutron evaluations, e.g. U235 H1"
    )
    ap.add_argument(
        "--thermal",
        nargs="*",
        default=[],
        metavar="MAT",
        help="thermal scattering laws, e.g. HinH2O",
    )
    ap.add_argument("--isotopes", metavar="ELEMENT", help="list nuclides of an element")
    ap.add_argument(
        "--list-thermal", action="store_true", help="list thermal materials"
    )
    ap.add_argument(
        "--root", help="data directory (default: $ENDF_DATA_DIR or ./endf-data)"
    )
    args = ap.parse_args(argv)

    lib = EndfLibrary(args.root)
    if args.isotopes:
        print(" ".join(lib.isotopes(args.isotopes)))
    if args.list_thermal:
        print("\n".join(lib.thermal_materials()))
    for path in lib.neutrons(args.nuclides) if args.nuclides else []:
        print(path)
    for name in args.thermal:
        print(lib.thermal(name))


if __name__ == "__main__":
    main()
