"""Summarize an ENDF/B-VIII.1 evaluation: header and MF3 cross sections.

    uv run python -m nodal.endf U238
    uv run python -m nodal.endf H1 --energy 0.0253 1e6
"""

import argparse

from ..data import EndfLibrary
from .evaluation import LRP_MEANING, Evaluation

THERMAL = 0.0253  # eV, neutron energy at 2200 m/s


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(prog="python -m nodal.endf", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("nuclide", help="e.g. U238, H1, O16")
    ap.add_argument("--energy", type=float, nargs="+", default=[THERMAL, 1e3, 1e6],
                    help="energies (eV) at which to print each cross section")
    args = ap.parse_args(argv)

    ev = Evaluation.from_file(EndfLibrary().neutron(args.nuclide))
    h = ev.header
    print(f"{h.zsymam}   MAT {ev.mat}   ZA {h.za}   AWR {h.awr:.6g}   ({h.sublibrary})")
    print(f"  lab {h.lab}, {h.edate}, {h.authors}")
    print(f"  fissionable: {'yes' if h.lfi else 'no'}   E max: {h.emax:.3g} eV   "
          f"ENDF-{h.nfor}, library version {h.nver}.{h.lrel}")
    print(f"  LRP = {h.lrp}: {LRP_MEANING.get(h.lrp, '?')}")
    print(f"  {len(h.description)} description lines, {len(ev.sections)} sections\n")

    cols = "".join(f"{f'{e:.4g} eV':>13}" for e in args.energy)
    print(f"  {'MT':>4}  {'reaction':<24}{'points':>7}{'threshold':>12}{'Q (MeV)':>10}{cols}")
    for xs in ev.cross_sections():
        values = "".join(f"{xs(e):13.5g}" for e in args.energy)
        print(f"  {xs.mt:4d}  {xs.name:<24}{len(xs.table):7d}{xs.threshold:12.4g}"
              f"{xs.qi / 1e6:10.3f}{values}")
    print("\n  cross sections in barns")


if __name__ == "__main__":
    main()
