"""Feed the C3 two-group constants to the old nodal code (formerly branch claude-nodal).

    make nodal-check STATE=1A

Runs in the tools environment (the nodal code needs scipy, which OpenMC's
environment already has) on the copy in endf-data/nodal-src (the branch
claude-nodal is deleted; that gitignored extract is the only copy left).
Three checks: an all-reflective core reproduces the few-group k-infinity; the
same with D B^2 added to the absorption gives k = 1 at the critical buckling;
a bare square core whose geometric buckling is B^2 shows the effect of the
extrapolation distance.
"""

import json
import math
import sys
from pathlib import Path

import numpy as np
from nodal import Core, Material, NodalSolver, refine

data = json.loads(Path(sys.argv[1]).read_text())
m = data["material"]
b2 = data["buckling_cm-2"]
reflective = {s: "reflective" for s in ("x-", "x+", "y-", "y+", "z-", "z+")}


def solve(mat, widths, bc):
    shape = tuple(len(w) for w in widths)
    core = Core(widths, np.zeros(shape, int), [mat], bc)
    return NodalSolver(core).solve(tol_k=1e-10, tol_phi=1e-9).k


base = Material(**m)
k_inf = solve(base, [[1.26] * 3, [1.26] * 3, [1.0]], reflective)
print(
    f"reflective 3x3 pin nodes: k = {k_inf:.6f} (few-group k-inf {data['k_infinity']:.6f})"
)
leaky = Material(
    **{**m, "sigma_a": (np.array(m["sigma_a"]) + np.array(m["D"]) * b2).tolist()}
)
print(
    f"with D B^2 (B^2 = {b2:.6g}): k = {solve(leaky, [[1.26] * 3, [1.26] * 3, [1.0]], reflective):.6f}"
)
side = math.pi * math.sqrt(2.0 / b2)  # bare square whose geometric buckling is B^2
dx, _ = refine([side], 5.0)
bare = {**reflective, "x-": "vacuum", "x+": "vacuum", "y-": "vacuum", "y+": "vacuum"}
print(
    f"bare square {side:.1f} cm, vacuum sides: k = {solve(base, [dx, dx, [1.0]], bare):.5f}"
    " (above 1: the flux vanishes only at the extrapolated edge, so the core acts larger)"
)
