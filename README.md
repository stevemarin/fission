# nodal

A compact multigroup, 3D Cartesian nodal diffusion code in Python (numpy + scipy),
built the way production LWR core simulators are: a nodal expansion method (NEM)
for accuracy, with coarse-mesh finite difference (CMFD) for the global solve.

```
uv run python benchmarks/iaea.py --dim 2 --h 20         # IAEA-2D, one node per assembly
uv run python benchmarks/iaea.py --dim 3 --h 10 --hz 20 # IAEA-3D
uv run python benchmarks/iaea.py --dim 2 --method fdm   # same mesh, plain finite difference
uv run pytest
```

## Method

| Piece | What it does |
|---|---|
| `nem.py` | NEM one-node kernel: 4th-order polynomial per direction, weighted-residual 1st/2nd moments, quadratic transverse leakage (one-sided at vacuum edges), partial-current (response-matrix) coupling, ADFs |
| `cmfd.py` | Partial-current CMFD (pCMFD) global eigenproblem, Wielandt-shifted power iteration, ILU-preconditioned BiCGSTAB with a reused preconditioner |
| `solver.py` | Nonlinear iteration: CMFD gives NEM node fluxes, k and incoming partial currents; NEM sweeps return partial currents that set the pCMFD corrections |
| `geometry.py` | Structured mesh, material map with `-1` for outside nodes (staircase core edges), albedo boundary conditions (`reflective`, `vacuum`, `zero_flux`, or a number) |

Why pCMFD rather than classic net-current CMFD: the classic form goes unstable
once nodes are several diffusion lengths thick. For the thermal group in a 20 cm
PWR node that is already the case, and it diverges on IAEA-2D. pCMFD stays
stable at any node size.

## Results

| Case | Mesh | k_eff | Error vs ref | FDM error, same mesh | Time |
|---|---|---|---|---|---|
| IAEA-2D (ref 1.029585) | 20 cm | 1.029486 | −9.9 pcm | +249.5 pcm | 0.1 s |
| IAEA-2D | 10 cm | 1.029606 | +2.1 pcm | | 0.1 s |
| IAEA-3D (ref 1.02903) | 20 cm | 1.028966 | −6.4 pcm | +263.7 pcm | 0.9 s |
| IAEA-3D | 10 cm radial × 20 cm axial | 1.029087 | +5.7 pcm | | 3.3 s |
| Bare 2-group slab, analytic | 50 cm nodes | | < 0.1 pcm | +75 pcm | |
| Bare 2-group square, analytic | 20 cm nodes | | −0.35 pcm | +56.5 pcm | |

## Using it

```python
import numpy as np
from nodal import Core, Material, NodalSolver

fuel = Material.two_group(D1=1.5, D2=0.4, sa1=0.01, sa2=0.08, s12=0.02, nsf1=0.0, nsf2=0.135)
refl = Material.two_group(2.0, 0.3, 0.0, 0.01, 0.04, 0.0, 0.0)
matmap = np.zeros((10, 10, 1), int); matmap[8:, :] = 1; matmap[:, 8:] = 1
core = Core([[20.0] * 10, [20.0] * 10, [1.0]], matmap, [fuel, refl],
            bc={"x+": "vacuum", "y+": "vacuum"})   # unspecified sides are reflective
r = NodalSolver(core).solve(verbose=True)
r.k, core.to_grid(r.power)
```

`Material` takes general G-group data (`sigma_s[g_from, g_to]`, `chi`,
`kappa_sigma_f`, `adf`). For a 2D problem use `dz=[1.0]`, leave z reflective,
and put any axial buckling into the absorption (`D·B²`).

## Not included (natural next steps)

Cross-section feedback (T/H, boron, burnup), pin power reconstruction,
time-dependent kinetics, hexagonal geometry, and a per-face ADF (currently one
per node per group).
