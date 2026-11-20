# The state matrix and its tables

Chapter 18 produced two-group constants for one cell at one state. A core
is not at one state. Every node of a nodal calculation has its own fuel
temperature (it follows the local power), its own coolant temperature and
density (the coolant heats up along the channel), and the same soluble boron
everywhere (the operators change it through the cycle). The nodal code
needs constants at each node's conditions. The lattice code cannot run
inside every node at every iteration, so it computes the constants once on
a grid of conditions, and the nodal code interpolates. That grid is the
**state matrix**, and each of its points is a **branch**.

## The four variables and their feedbacks

The state variables are also where a PWR's reactivity feedbacks come from:

* **Fuel temperature.** The U-238 resonances broaden (chapter 8). With the
  pellet's flux shielding them (chapters 12 and 16), the broadened
  resonances absorb more, and $k$ drops. This is the Doppler feedback: it
  acts as soon as the power rises, before any heat reaches the water.
* **Moderator density.** Less water moderates fewer neutrons to thermal
  energies, but it also absorbs fewer. A PWR lattice is *under-moderated*
  (it has less water than would maximize $k$), so losing water lowers $k$.
* **Moderator temperature at fixed density.** A hotter Maxwellian moves
  the thermal spectrum up (chapter 11), where U-235 fission and the $1/v$
  absorbers are weaker.
* **Soluble boron.** B-10 is a thermal absorber in the water. Its worth
  depends on the water density, because the boron is carried by the water.

The last point couples two axes. With much boron, removing water also
removes absorber, and the moderator feedback can change sign.

## Branches in the model file

A state matrix is declared next to the states it starts from (chapter 13):

```toml
[branches.vera]
base = "1B"
fuel_temperature = [565.0, 600.0, 900.0, 1200.0]
moderator_temperature = [565.0, 600.0]
moderator_density = [0.661, 0.743]
boron_ppm = [0.0, 1300.0, 2000.0]
```

The points are the **tensor product** of the axes: 4 × 2 × 2 × 3 = 48 here.
A full grid is what multilinear interpolation (M3) needs. It is also what
makes the cost grow quickly: every new value of an axis multiplies the
number of points. An axis the file leaves out keeps the base state's value.
The loader checks that every list increases, because the values are
interpolation nodes. It also checks that every list contains the base value,
because the coefficients below are taken around it. Moderator temperatures
must be temperatures of the H in H2O scattering law, which is tabulated every
5 K from 273 to 650 K. All four VERA states are grid points of this matrix
(1A is fuel 565 K, water 565 K at 0.743 g/cm³; 1B–1D are fuel 600, 900 and
1200 K with water at 600 K and 0.661 g/cm³; all at 1300 ppm).

## The driver

Each point is the lattice calculation of chapter 18:

1. build the state, and from it the number densities and temperatures of
   every material;
2. solve the lattice, either the homogenized cell with narrow-resonance
   shielding of chapter 12 (`infinite`, the default) or the heterogeneous pin
   of chapters 16–19 (`moc`);
3. homogenize to fine-group constants;
4. find the B1 critical spectrum and condense to two groups.

For the infinite medium, step 3 has nothing to weight: the mixture is
already homogeneous, and its constants are its own cross sections. Only
$\chi$ needs a weight. Each fissile nuclide has its own fission spectrum
$\chi_i$, and the medium's spectrum is their average weighted by fission
production, $\chi = \sum_i \chi_i P_i / \sum_i P_i$ with
$P_i = \nu\Sigma_{f,i}\cdot\phi$. With the converged flux this reproduces
the fission source exactly, $\sum_i \chi_i P_i = \chi \sum_i P_i$, so the
constants have the solver's $k_\infty$ (to 10⁻¹⁵).

The table is rewritten after every point, so an interrupted run resumes
where it stopped. A change of solver, options, spectrum or axes starts a
new table. The self-shielded libraries of a point are deleted after it:
every point has its own flux shape, so they would only fill the cache with
about 2 MB per nuclide and point.

## Results for VERA problem 1

The infinite-medium solver takes about 90 s per point. Nearly all of it is
building the 40 self-shielded libraries, because the narrow-resonance flux
shape changes with every composition and temperature. The whole matrix took
76 minutes. The points that are VERA states reproduce the standalone runs
(1A: $k_\infty$ = 1.13960).

Reactivity coefficients at the base state 1B, from the neighbouring grid
points ($\rho = (k-1)/k$, central differences where the base has two
neighbours):

| variable | coefficient | between |
| --- | --- | --- |
| fuel temperature | −3.96 pcm/K | 565 and 900 K |
| moderator temperature (fixed density) | −1.73 pcm/K | 565 and 600 K |
| moderator density | +22.5 pcm per 0.01 g/cm³ | 0.661 and 0.743 g/cm³ |
| boron | −5.94 pcm/ppm | 0 and 2000 ppm |

All signs are those of a safe PWR lattice at 1300 ppm. The density
coefficient is positive: more water, more $k$, so a lattice that loses
water (heats up, boils) loses reactivity. Boron worth falls as boron is
added (−7784 pcm for the first 1300 ppm, −4104 for the next 700), because
each added atom competes with those already there for the thermal neutrons.

The coupling of density and boron shows in the table. Raising the water
density from 0.661 to 0.743 g/cm³ (fuel and water at 600 K):

| boron | 0 ppm | 1300 ppm | 2000 ppm |
| --- | --- | --- | --- |
| $\Delta\rho$ | +1092 pcm | +185 pcm | −293 pcm |

At 2000 ppm the extra water brings enough boron to outweigh its moderation:
the moderator feedback is positive. This is why a PWR's boron
concentration at the start of a cycle is limited.

**The homogeneous Doppler is too strong.** From 1B to 1D the homogenized
cell loses 2201 pcm. KENO's heterogeneous pin loses 1422 pcm (the VERA
references, chapter 13). In the homogenized mixture the U-238 resonances
are shielded only by the mixture's own cross section. The pellet's spatial
self-shielding (chapter 16) is gone, so there is more resonance absorption
for Doppler broadening to increase. The infinite-medium matrix shows how
the driver works and gives the right shape of the feedbacks. For constants
that a core calculation can trust, the same matrix should run with the
heterogeneous solver:

```sh
make branches BRANCH_ARGS="--solver moc --rings --sph --clad"
```

At the converged MOC defaults (chapter 19) that is about 3.5 minutes per
point once the libraries are cached, three hours for this matrix.

## What the table holds

For every point (`endf-data/branches/<model>-<branches>-<solver>.json`):
the state, the lattice $k_\infty$, the critical buckling, and the two-group
$D$, $\Sigma_t$, $\Sigma_a$, $\Sigma_{s,g\to g'}$, $\nu\Sigma_f$,
$\Sigma_f$, $\kappa\Sigma_f$, $\chi$ and the discontinuity factors (1 for
the infinite medium). `Table.array(name)` returns any of them shaped like
the grid, (4, 2, 2, 3) plus the quantity's own shape, which is the form M3
interpolates.

## Constants between the nodes

The nodal code needs constants at any state inside the grid. Two ways to
get them are compared here (`fission.tables`, `make xs-tables`).

**Multilinear interpolation.** In one variable, between nodes $x_i$ and
$x_{i+1}$,

$$
\Sigma(x) = (1-t)\,\Sigma_i + t\,\Sigma_{i+1}, \qquad t = \frac{x - x_i}{x_{i+1}-x_i}.
$$

In four variables the same is done axis after axis. Each of the 16 corners
of the grid cell around the state gets a weight equal to the product of its
four factors $(1-t)$ or $t$. Outside the grid the end segment is continued.
The result is linear in the tabulated data. So any linear relation that
holds at every node, such as $\Sigma_t = \Sigma_a + \sum_{g'}\Sigma_{s,g\to g'}$,
also holds at every interpolated state, exactly.

**Which variable?** Interpolation is only as good as the data is linear in
the variable used. Doppler broadening widens a resonance by
$\Gamma_D \propto \sqrt{T}$ (chapter 8), and the shielded resonance
integral grows roughly like $\sqrt{T}$. So the fuel axis is interpolated in
$\sqrt{T_f}$. Leaving out one interior node and interpolating it from its
neighbours shows the difference in the lattice $k$:

| left out | in $T_f$ | in $\sqrt{T_f}$ |
| --- | --- | --- |
| $T_f$ = 600 K (from 565 and 900 K) | +15.9 pcm | +2.3 pcm |
| $T_f$ = 900 K (from 600 and 1200 K) | +93.8 pcm | +13.0 pcm |
| boron 1300 ppm (from 0 and 2000) | +63.5 pcm | |

No single variable is right for every constant. The thermal-to-fast
upscatter (0.0019 cm⁻¹ at 565 K, 0.0024 at 1200 K) comes from the thermal
motion of the fuel nuclei, whose energy is $kT$. It is linear in $T$: its
leave-one-out error at 900 K is 0.02% in $T$ and 0.9% in $\sqrt{T}$. It
hardly affects $k$; the resonance absorption does, and it wants $\sqrt{T}$.

**A fitted polynomial.** The other approach fits each constant with a few
terms,

$$
\Sigma(u) = \sum_j c_j\, m_j(u),
$$

where the $m_j$ are monomials of the variables mapped to $[-1, 1]$ (fuel in
$\sqrt{T_f}$). The candidates are the linear terms, squares and cubes where
the grid has enough nodes, and products of two variables. Terms are added
greedily. At each step the candidate that most lowers the leave-one-out
error over the 48 nodes wins. That error needs no refitting. For least
squares with the hat matrix $H = X(X^TX)^{-1}X^T$, the residual at node $i$
when $i$ is left out is $e_i/(1-H_{ii})$ (the PRESS residual). One set of
terms serves all constants, so the linear relations between constants hold
for the fit too. The order in which terms are chosen is itself a summary of
the physics:

| step | term added | relative leave-one-out error |
| --- | --- | --- |
| 1 | constant | 4.5% |
| 2 | $\rho$ | 3.3% |
| 3 | $T_f$ | 2.3% |
| 4 | $B$ | 1.1% |
| 5 | $T_m$ | 0.25% |
| 6 | $\rho\,B$ | 0.20% |
| 7 | $T_f\,B$ | 0.17% |
| 8 | $T_f^2$ | 0.13% |
| 9–13 | $T_f\rho$, $T_mB$, $T_m\rho$, $B^2$, $T_fT_m$ | 0.010% |

Water density moves the constants most, since they scale with the water's
number densities. Next come the Doppler effect and boron, then the moderator
temperature. The first cross term is density × boron, the coupling seen
above. The selection stops at the fewest terms within 1% of the best score,
here 13 coefficients per constant instead of 48 nodes.

**Check points.** Ten random states inside the box, computed by the lattice
code (`make branches BRANCH_ARGS="--check 10"`), test both methods off the
grid:

| method | largest error in $k$ | largest error in a constant |
| --- | --- | --- |
| multilinear, fuel in $T_f$ | +56.7 pcm | 0.21% ($D$) |
| multilinear, fuel in $\sqrt{T_f}$ | +28.7 pcm | 0.35% ($\Sigma_s$) |
| fit, 13 terms | +4.2 pcm | 0.21% ($D$) |

The multilinear errors all have the same sign, and they are largest where
the boron is between 0 and 1300 ppm. Boron worth saturates (−7784 pcm for
the first 1300 ppm, −4104 for the next 700). The straight line between
those two nodes lies below the real curve of absorption, so it puts too
little boron worth between them. The fit's $B^2$ term captures that
curvature. A third boron node (650 ppm) would bring multilinear
interpolation closer, at the price of 16 more lattice calculations.

The table for the nodal code is written to
`endf-data/tables/<model>-<branches>-<solver>.json`: axes, transforms, and
every constant on the grid. Discontinuity factors are interpolated like the
cross sections; they are 1 for the infinite medium, and the MOC solver
provides real ones (chapter 18).

## Further reading

* [[Hébert](#ref-hebert)] (not free): multiparameter reconstruction of
  cross sections in the lattice-to-core chain.
* Code: [`tables/__init__.py`](../../src/fission/tables/__init__.py)
  (multilinear interpolation, PRESS selection, fit),
  [`tables/__main__.py`](../../src/fission/tables/__main__.py) (tests,
  `make xs-tables`).
* Code: [`branches/__init__.py`](../../src/fission/branches/__init__.py)
  (driver, table, coefficients),
  [`branches/__main__.py`](../../src/fission/branches/__main__.py) (CLI,
  `make branches`), [`model/spec.py`](../../src/fission/model/spec.py)
  (`Branches`), `from_mixture` in
  [`fewgroup/homogenize.py`](../../src/fission/fewgroup/homogenize.py).
