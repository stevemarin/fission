"""IAEA PWR benchmarks (ANL-7416 Suppl. 2, problem 11), 2-group, quarter core.

    uv run python benchmarks/iaea.py --dim 2 --h 20
    uv run python benchmarks/iaea.py --dim 3 --h 10 --hz 20
    uv run python benchmarks/iaea.py --dim 2 --h 20 --method fdm

Quarter core: 9 x 9 assemblies, 20 cm pitch, the first row/column being half
assemblies cut by the symmetry planes (10 cm). Reflective on x-/y-; on the
outer surfaces d(phi)/dn = -0.4692 phi / D, i.e. an albedo of 0.4692.
"""

import argparse

import numpy as np

from nodal import Core, Material, NodalSolver, refine

K_REF = {2: 1.029585, 3: 1.02903}
ALBEDO = 0.4692

# 1 = fuel 1, 2 = fuel 2, 3 = fuel 2 + rod, 4 = reflector, 0 = outside.
# Row index = y assembly, column index = x assembly, (0, 0) at the core centre.
ASSEMBLY_MAP = np.array([
    [3, 2, 2, 2, 3, 2, 2, 1, 4],
    [2, 2, 2, 2, 2, 2, 2, 1, 4],
    [2, 2, 2, 2, 2, 2, 1, 1, 4],
    [2, 2, 2, 2, 2, 2, 1, 4, 4],
    [3, 2, 2, 2, 3, 1, 1, 4, 0],
    [2, 2, 2, 2, 1, 1, 4, 4, 0],
    [2, 2, 1, 1, 1, 4, 4, 0, 0],
    [1, 1, 1, 4, 4, 4, 0, 0, 0],
    [4, 4, 4, 4, 0, 0, 0, 0, 0],
])
PITCH = np.array([10.0] + [20.0] * 8)
# 3D: four partially inserted rods (80 cm from the top of the core).
PARTIAL_RODS = [(2, 2)]


def materials(buckling=0.0):
    """Regions 1-5 (index 0 unused); buckling adds D*B^2 to absorption."""
    def mat(D1, D2, sa1, sa2, s12, nsf2, name):
        return Material.two_group(D1, D2, sa1 + D1 * buckling, sa2 + D2 * buckling,
                                  s12, 0.0, nsf2, name=name)
    return [
        None,
        mat(1.5, 0.4, 0.01, 0.080, 0.02, 0.135, "fuel 1"),
        mat(1.5, 0.4, 0.01, 0.085, 0.02, 0.135, "fuel 2"),
        mat(1.5, 0.4, 0.01, 0.130, 0.02, 0.135, "fuel 2 + rod"),
        mat(2.0, 0.3, 0.00, 0.010, 0.04, 0.000, "reflector"),
        mat(2.0, 0.3, 0.00, 0.055, 0.04, 0.000, "reflector + rod"),
    ]


def build(dim, h, hz=None):
    """Return (core, ix_owner, iy_owner, iz_owner)."""
    dx, ox = refine(PITCH, h)
    radial = ASSEMBLY_MAP.T  # -> [ix, iy]
    if dim == 2:
        mats = materials(buckling=0.8e-4)
        matmap = radial[np.ix_(ox, ox)][:, :, None]
        dz, oz = np.array([1.0]), np.array([0])
        bc = {"x+": ALBEDO, "y+": ALBEDO}
    else:
        mats = materials()
        rodded = radial == 3
        partial = np.zeros_like(rodded)
        for i, j in PARTIAL_RODS:
            partial[i, j] = partial[j, i] = True
        fuelcol = (radial >= 1) & (radial <= 3)
        upper = np.where(partial, 3, radial)
        bottom = np.where(fuelcol, 4, radial)
        top = np.where(fuelcol, np.where(rodded | partial, 5, 4), radial)
        # axial: 20 cm reflector, 260 cm core, 80 cm core (partial rods), 20 cm reflector
        layers = [bottom, radial, upper, top]
        dz, oz = refine([20.0, 260.0, 80.0, 20.0], hz or h)
        stack = np.stack(layers, axis=2)  # [ix, iy, layer]
        matmap = stack[np.ix_(ox, ox, oz)]
        bc = {"x+": ALBEDO, "y+": ALBEDO, "z-": ALBEDO, "z+": ALBEDO}
    matmap = matmap - 1  # region r -> material index r-1, outside -> -1
    core = Core([dx, dx, dz], matmap, mats[1:], bc)
    return core, ox, ox, oz


def assembly_power(core, result, ox, oy):
    """Radial assembly powers (axially integrated), fuel-area mean = 1. NaN = no fuel."""
    nx = ox.max() + 1
    num = np.zeros((nx, nx))
    area = np.zeros((nx, nx))
    i, j, kz = core.ijk.T
    np.add.at(num, (ox[i], oy[j]), result.power * core.vol)
    fuel = core.nu_sigma_f.sum(axis=1) > 0
    a_xy = core.h[:, 0] * core.h[:, 1]
    first_plane = kz == kz[fuel].min()  # count radial area once per (x, y) node
    np.add.at(area, (ox[i[fuel & first_plane]], oy[j[fuel & first_plane]]),
              a_xy[fuel & first_plane])
    with np.errstate(invalid="ignore", divide="ignore"):
        p = num / area
    p[area == 0] = np.nan
    w = PITCH[:, None] * PITCH[None, :]
    ok = ~np.isnan(p)
    return p / ((p * w)[ok].sum() / w[ok].sum())


def print_map(p):
    print("  assembly power (row = y from centre, col = x from centre)")
    for row in p.T:
        print("   " + " ".join("  ---" if np.isnan(v) else f"{v:5.3f}" for v in row))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("--dim", type=int, choices=(2, 3), default=2)
    ap.add_argument("--h", type=float, default=20.0, help="radial node size (cm)")
    ap.add_argument("--hz", type=float, default=None, help="axial node size (cm), 3D only")
    ap.add_argument("--method", choices=("nem", "fdm"), default="nem")
    ap.add_argument("-v", "--verbose", action="store_true")
    args = ap.parse_args()

    core, ox, oy, _ = build(args.dim, args.h, args.hz)
    print(f"IAEA-{args.dim}D  {args.method.upper()}  mesh {core.shape}  "
          f"{core.N} nodes x {core.G} groups")
    r = NodalSolver(core, args.method).solve(verbose=args.verbose)
    err = (r.k - K_REF[args.dim]) * 1e5
    print(f"  k_eff = {r.k:.6f}   reference {K_REF[args.dim]:.6f}   error {err:+.1f} pcm")
    print(f"  {r.nodal_iterations} nodal updates, {r.outer_iterations} CMFD outers, "
          f"{r.seconds:.2f} s, converged={r.converged}")
    print_map(assembly_power(core, r, ox, oy))


if __name__ == "__main__":
    main()
