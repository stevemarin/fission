# Nodal diffusion of the core

A PWR core has about 200 assemblies, each 4 m tall. Transport in every pin
(chapter 19 did one assembly in two dimensions) is far too expensive for
the many core calculations of a cycle. The lattice code therefore hands the
core code few-group constants per assembly (chapters 18 and 20), and the
core code solves **two-group diffusion** on coarse **nodes**: one per
assembly radially and about 20 cm axially. This chapter builds that solver
and checks it on the classic benchmark designed to break coarse-mesh
methods.

## Two-group diffusion

In each group the neutron balance with Fick's law, $J = -D\nabla\phi$, reads

$$
-\nabla\cdot D_g\nabla\phi_g + \Sigma_{r,g}\,\phi_g
= \frac{\chi_g}{k}\sum_{g'}\nu\Sigma_{f,g'}\phi_{g'} + \sum_{g'\neq g}\Sigma_{s,g'\to g}\,\phi_{g'},
$$

with the removal $\Sigma_r$ = absorption + out-scattering. At a vacuum
boundary no neutrons come back in. In diffusion theory the incoming partial
current is $J^- = \phi/4 - (D/2)\,\partial\phi/\partial n$ (Marshak), so
$J^- = 0$ becomes the Robin condition

$$
\frac{\partial\phi}{\partial n} = -\frac{c}{D}\,\phi, \qquad c = \tfrac12 .
$$

The IAEA benchmark prescribes $c = 0.4692$ instead, the transport
extrapolation distance $0.7104\,\lambda_{tr}$ written for diffusion codes
(the model file's `vacuum_coefficient`).

## The core as nodes

`nodal.core_mesh` turns the model's core (chapter 13) into a 3D grid.
Every map position becomes `refine` × `refine` nodes. The axial node
boundaries are those of all assembly stacks (rods inserted part way, for
example), with nodes of at most one radial node width in between. With
quarter symmetry the first row and column lie on the symmetry planes,
through the middle of the centre assembly: those nodes are half as wide,
and their outer faces are reflective. Positions outside the map are void:
their faces get the vacuum condition.

## Finite difference: the reference

Integrating the balance over node $i$ gives
$\sum_f A_f J_f + \Sigma_r V\phi_i = V\,Q_i$. If the flux is linear inside
each node and the current is continuous at a face, the current between
nodes lo and hi (hi further along the axis) is

$$
J = -\tilde D\,(\phi_\text{hi} - \phi_\text{lo}), \qquad
\tilde D = \frac{2 D_\text{lo} D_\text{hi}}{D_\text{lo} h_\text{hi} + D_\text{hi} h_\text{lo}},
$$

and on a vacuum face, eliminating the surface flux with the Robin condition,
$J_\text{out} = \tilde D_b\,\phi$ with $\tilde D_b = 2Dc/(2D + ch)$. The
eigenvalue comes from power iteration, and each group's equation from
red-black successive over-relaxation. The nodes are coloured like a 3D
chess board, so that each colour depends only on the other.

A bare homogeneous one-group cube has an exact answer. With a reflective
plane at 0 and a vacuum face at $a$, the flux along each axis is
$\cos(Bx)$ with $B\tan(Ba) = c/D$, and $k = \nu\Sigma_f/(\Sigma_a + D\sum B^2)$.
For a 60 × 60 × 60 cm quarter of such a cube (tests), finite difference is
off by 849, 224, 57 and 14 pcm with 30, 15, 7.5 and 3.75 cm nodes, falling
by a factor of four per halving: second order, as it should be.

**The IAEA 3D benchmark** (ANL-7416 Supplement 2, problem 11, reproduced
in [[FeenoX](#ref-feenox)]) is a quarter core of 20 cm assemblies, 340 cm
of core between 20 cm reflectors, two fuel types, five rods (one of them
only in the top 80 cm) and a reflector, all in two groups. Its strong
thermal-flux gradients at the reflector and the rods were chosen to break
coarse-mesh methods, and finite difference shows why:

| node size | 20 cm | 10 cm | 5 cm | 2.5 cm | 1.25 cm |
| --- | --- | --- | --- | --- | --- |
| finite difference, $k$ | 1.031666 | 1.029059 | 1.028662 | 1.028895 | 1.029022 |

The error is not even monotonic: different regions err in different
directions, and they cancel differently on each mesh. A 1.25 cm mesh is 5.6
million nodes in a quarter core and takes 7.5 minutes (with the sweeps in
Rust, below); it ends 5 pcm from the converged NEM value.

## The nodal expansion method

NEM keeps 20 cm nodes and describes the flux inside each one better.
Integrating the diffusion equation over the two transverse directions of a
node gives a one-dimensional equation along $x$ (local coordinate
$\xi \in [-\tfrac12, \tfrac12]$, node width $h$):

$$
-\frac{D}{h^2}\,\phi''(\xi) + \Sigma_r\,\phi(\xi) = Q(\xi) - L(\xi),
$$

where $Q$ is the fission and in-scattering source and $L$ the
**transverse leakage**, what flows out through the $y$ and $z$ faces per
unit volume. The flux is expanded in polynomials of zero mean,

$$
\phi(\xi) = \bar\phi + a_1P_1 + a_2P_2 + a_3P_3 + a_4P_4,\quad
P_1 = \xi,\; P_2 = 3\xi^2 - \tfrac14,\; P_3 = \xi(\xi^2 - \tfrac14),\;
P_4 = (\xi^2 - \tfrac1{20})(\xi^2 - \tfrac14).
$$

$P_3$ and $P_4$ vanish at the faces, so they bend the flux without changing
its face values. The source is a combination of the group fluxes, so it is
expanded in the same polynomials. The transverse leakage is taken as a
quadratic, $L = \bar L + l_1P_1 + l_2P_2$, the one that keeps the average
leakage of the node and of its two neighbours (a missing neighbour repeats
the node's own value, which is exact at a symmetry plane). Weighting the
equation with $P_1$ and $P_2$ (the first two moments) gives two equations
per group for the four coefficients.

**Two-node problems.** Across an interface the unknowns are $a_1..a_4$ of
both nodes in every group. There are eight equations per group: the four
moments; continuity of the flux at the interface (with discontinuity
factors $f$: $f_\text{lo}\phi_\text{lo} = f_\text{hi}\phi_\text{hi}$); continuity
of the current, $J = -(D/h)\,\phi'$; and the currents on the two far faces,
taken from the coarse problem. The groups are coupled through the source,
so each interface is one small linear system ($8G$ unknowns), and all
interfaces are solved together as a batch. On an external face the same is
done with one node and the boundary condition in place of the neighbour.

**Coarse mesh finite difference.** The current a two-node problem finds is
put back into the finite-difference problem through a correction $\hat D$:

$$
J = -\tilde D\,(\phi_\text{hi} - \phi_\text{lo}) - \hat D\,(\phi_\text{hi} + \phi_\text{lo}).
$$

The coarse problem with $\hat D$ reproduces the NEM currents exactly and
gives the new eigenvalue and node averages, which feed the next two-node
problems (Smith's non-linear iteration). This is the same idea as the CMFD
acceleration of MOC in chapter 15: the expensive method computes currents,
and the cheap one turns them into an eigenvalue.

**Two things that had to be fixed.** The first version took, on the far
faces, the face *fluxes* found by the previous two-node problems. That
lagged one problem behind the other, and the iteration oscillated and
diverged, even in one dimension; under-relaxation only delayed it. The
coarse problem's *currents* satisfy every node's balance by construction,
and with them as far-face conditions the iteration converges in about 15
steps without relaxation. The second fix: where $|\hat D| > \tilde D$, a
coefficient of the coarse operator turns negative (the fuel-reflector and
rod interfaces of the IAEA core). There the face current is written
one-sidedly instead: for $J > 0$, $\tilde D' = -\hat D' = J/(2\phi_\text{lo})$.
The current is the same, so the converged answer is too, and every
coefficient stays non-negative.

A third fix came later, with the Pyrex core of chapter 22. Its axial mesh
has a 1.27 cm layer between 21 cm nodes, under assemblies 21.5 cm wide.
There the iteration oscillated with period two: $k$ alternated between
0.99119 and 0.99130, the flux change grew from 4.4% to 6% over ten
steps, and then it overflowed. The current through a thin node's faces is
very sensitive to its neighbours' fluxes, and each update overshot. Taking
only half of each new $\hat D$, $\hat D \leftarrow \tfrac12 \hat D_\text{old}
+ \tfrac12 \hat D_\text{new}$, converges it in 16 steps. The cores that
converged before take about as many steps as without it (IAEA 3D 18 and 18,
the zoned core 16 and 17), so half steps are now the default
(`nem.RELAX`).

## Results

The exact cube with NEM:

| | 30 cm | 15 cm | 7.5 cm |
| --- | --- | --- | --- |
| 1D (only $x$ leaks) | 0.005 pcm | | |
| 2D | +54 pcm | +5.6 pcm | +0.7 pcm |
| 3D | +510 pcm | +75 pcm | +9.4 pcm |

In one dimension the fourth-order polynomial is all but exact. With
transverse leakage the error falls about eightfold per halving (third
order). This small cube is a severe case: two nodes per direction, with the
flux falling to 6% of its maximum at the vacuum faces.

The IAEA 3D benchmark with NEM:

| node size | 20 cm | 10 cm | 5 cm |
| --- | --- | --- | --- |
| NEM, $k$ | 1.028988 | 1.029095 | 1.029074 |
| time | 0.6 s | 2.8 s | 28 s |

NEM with one node per assembly is 9 pcm from its own converged value; finite
difference on the same mesh is 260 pcm off. FeenoX's second-order finite
elements give 1.02918 (quarter core, 70,779 nodes) and 1.02912 (eighth
core, 47,798 nodes) [[FeenoX](#ref-feenox)]. Their difference shows that
mesh error is still of that size, and our converged value lies 5–11 pcm
below them.

The radial assembly powers (average power density, normalized to 1 over
the full core; $\nu\Sigma_f\phi$ for the power, since both fuels have the
same $\nu\Sigma_f$), quarter core, centre at top left:

```
0.734  1.287  1.427  1.199  0.611  0.955  0.955  0.771
1.287  1.403  1.438  1.295  1.075  1.054  0.970  0.750
1.427  1.438  1.374  1.317  1.184  1.088  0.993  0.703
1.199  1.295  1.317  1.182  0.975  0.919  0.858
0.611  1.075  1.184  0.975  0.472  0.697  0.604
0.955  1.054  1.088  0.919  0.697  0.587
0.955  0.970  0.993  0.858  0.604
0.771  0.750  0.703
```

The rods show as dips: the centre (0.734), the assemblies on the axes
(0.611), and the diagonal one (0.472). The maximum, 1.438, lies next to the
partially rodded assembly.

## Rust

Profiling showed one hot spot: the red-black SOR sweep took 85–95% of the
time of both methods (NEM spends its time in the coarse finite-difference
problem, not in the two-node problems, which are batched linear solves).
`rust/src/diffusion.rs` sweeps one group in place. It evaluates the same
expression as the Python reference, in the same order, so the result is
identical bit for bit. A half-sweep only updates nodes of one colour, from
neighbours of the other, so on large meshes it is split over threads by
planes of $z$ without changing a single bit (tests: `test_diffusion_sweeps`,
on one thread and on several). With the default `use_rust=True`:

| | Python | Rust |
| --- | --- | --- |
| NEM, 20 cm nodes | 2.5 s | 0.6 s |
| NEM, 10 cm nodes | 10 s | 2.8 s |
| finite difference, 5 cm | 23 s | 4.8 s |
| finite difference, 2.5 cm | 269 s | 34 s |

The eigenvalues are unchanged to every digit. The rest of the time is the
Python outer loop and NumPy.

## What is still missing

* **Our own constants.** The IAEA core gives its constants directly.
  Chapter 22 builds a core from the M3 tables, with thermal-hydraulic
  feedback and reflector constants of its own.

## Further reading

* [[FeenoX](#ref-feenox)] (free): the IAEA 2D and 3D benchmarks with the
  original problem description (figures 3–5), and finite-element solutions.
* [[Hébert](#ref-hebert)] (not free): nodal methods and the transverse
  integration.
* Code: [`nodal/mesh.py`](../../src/fission/nodal/mesh.py),
  [`nodal/diffusion.py`](../../src/fission/nodal/diffusion.py) (finite
  difference, coarse operator), [`diffusion.rs`](../../rust/src/diffusion.rs), [`nodal/nem.py`](../../src/fission/nodal/nem.py)
  (two-node problems, non-linear iteration),
  [`nodal/__main__.py`](../../src/fission/nodal/__main__.py) (`make nodal`),
  `models/iaea-3d.toml`.
