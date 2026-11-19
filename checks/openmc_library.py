"""Build an OpenMC HDF5 library from our ENDF/B-VIII.1 files (independent check).

Runs in the tools environment, never in the pipeline:

    make openmc-library TEMPS="565 600 900 1200"

NJOY (reconr, broadr, heatr, gaspr, purr, acer) processes each evaluation in
endf-data/neutrons at the requested temperatures; ``openmc.data`` converts the
ACE files to HDF5 in endf-data/openmc/ and writes cross_sections.xml. The
thermal scattering law of H in H2O is processed at the temperatures it
tabulates among those requested.
"""

from __future__ import annotations

import argparse
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import openmc.data

ROOT = Path(__file__).resolve().parents[1] / "endf-data"
OUT = ROOT / "openmc"


def neutron(path: Path, temperatures: list[float]) -> tuple[str, float]:
    t0 = time.perf_counter()
    data = openmc.data.IncidentNeutron.from_njoy(
        str(path), temperatures=temperatures, error=0.001
    )
    data.export_to_hdf5(OUT / f"{data.name}.h5", "w")
    return data.name, time.perf_counter() - t0


def thermal(temperatures: list[float]) -> tuple[str, float]:
    t0 = time.perf_counter()
    data = openmc.data.ThermalScattering.from_njoy(
        str(ROOT / "neutrons" / "n-001_H_001.endf"),
        str(ROOT / "thermal_scatt" / "tsl-HinH2O.endf"),
        temperatures=temperatures,
        iwt=0,
    )
    data.export_to_hdf5(OUT / f"{data.name}.h5", "w")
    return data.name, time.perf_counter() - t0


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--temps", type=float, nargs="+", default=[565, 600, 900, 1200])
    ap.add_argument("--sab-temps", type=float, nargs="+", default=[565, 600])
    ap.add_argument("--only", nargs="*", help="file stems, e.g. n-092_U_238")
    ap.add_argument("--jobs", type=int, default=8)
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    files = sorted((ROOT / "neutrons").glob("n-*.endf"))
    if args.only:
        files = [f for f in files if f.stem in args.only]
    with ProcessPoolExecutor(args.jobs) as pool:
        jobs = {pool.submit(neutron, f, args.temps): f.name for f in files}
        if not args.only:
            jobs[pool.submit(thermal, args.sab_temps)] = "tsl-HinH2O"
        failed = []
        for job in as_completed(jobs):
            try:
                name, seconds = job.result()
                print(f"{name:10s} {seconds:6.0f} s", flush=True)
            except Exception as exc:  # noqa: BLE001  (report, continue with the rest)
                failed.append(jobs[job])
                print(f"FAILED {jobs[job]}: {exc}", flush=True)
    lib = openmc.data.DataLibrary()
    for h5 in sorted(OUT.glob("*.h5")):
        lib.register_file(h5)
    lib.export_to_xml(OUT / "cross_sections.xml")
    print(f"{len(lib.libraries)} files in {OUT / 'cross_sections.xml'}")
    if failed:
        sys.exit(f"failed: {failed}")


if __name__ == "__main__":
    main()
