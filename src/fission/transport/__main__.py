"""Track a model lattice and compare track-estimated areas with exact ones.

uv run python -m fission.transport vera-p1 --azim 64 --spacing 0.01
"""

from __future__ import annotations

import argparse
import time

import numpy as np

from fission.model import load_model

from .geometry import Geometry
from .trace import trace, track_areas
from .tracks import cyclic_tracks


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(prog="python -m fission.transport", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)  # fmt: skip
    ap.add_argument("model", help="TOML file or name in models/")
    ap.add_argument("--lattice", default="cell")
    ap.add_argument("--azim", type=int, default=64, help="directions in [0, 2 pi)")
    ap.add_argument("--spacing", type=float, default=0.01, help="track spacing (cm)")
    ap.add_argument("--fuel-rings", type=int, default=3)
    ap.add_argument("--sectors", type=int, default=8)
    args = ap.parse_args(argv)
    model = load_model(args.model)
    g = Geometry.from_lattice(
        model, args.lattice, {"fuel": args.fuel_rings}, args.sectors
    )
    t0 = time.perf_counter()
    tracks = cyclic_tracks(g.width, g.height, args.azim, args.spacing)
    t1 = time.perf_counter()
    seg = trace(g, tracks)
    t2 = time.perf_counter()
    exact = g.areas()
    rel = np.abs(tracks.weights @ track_areas(g, tracks, seg) / exact - 1)
    print(
        f"{model.name}, lattice {args.lattice}: {g.nx}x{g.ny} cells, {g.regions} regions"
    )
    print(f"  {len(tracks)} tracks in {args.azim} directions ({t1 - t0:.2f} s), "
          f"{len(seg.length)} segments ({t2 - t1:.2f} s)")  # fmt: skip
    a = tracks.angles[0]
    print(f"  first angle {np.degrees(a.phi):.3f} deg, spacing {a.spacing:.5f} cm, "
          f"nx {a.nx}, ny {a.ny}")  # fmt: skip
    print("  largest relative area error per material:")
    materials = np.array(g.materials())
    for m in dict.fromkeys(materials):
        print(f"    {m:10s} {rel[materials == m].max():.1e}")


if __name__ == "__main__":
    main()
