# Describing a reactor: the model file

Chapters 7–12 turned evaluated data into a multigroup library and used it on
one homogenized pin cell, whose geometry and number densities were written
into the code. A real calculation needs many more inputs: for a core, a
few-group nodal model is built from lattice calculations of every kind of
fuel assembly at many operating conditions. Keeping that description
consistent by hand, one input deck per code, is where most errors come from.
Modern workflows describe the reactor **once, as data**, and generate
every calculation from that description. VERA's common input format
(VERAIn) drives both its transport and thermal-hydraulics codes from one
deck [[VERA](#ref-vera)]. TerraPower's ARMI framework keeps one hierarchical
reactor model that every solver plugin reads [[ARMI](#ref-armi)]. This
chapter describes this code's model file: what it contains, how number
densities follow from it, and what it leaves out.

## Layers

The file has five layers. Each refers only to the ones before it:

| Layer | What it is | Example (VERA) |
| --- | --- | --- |
| materials | compositions only | UO$_2$, Zircaloy-4, borated water |
| pins | concentric rings in a square cell | fuel rod: pellet, gap, clad, water |
| lattices | one **2D** xy slice: a square map of pins, or one homogeneous material | the 17×17 fuel lattice |
| assemblies | lattices **stacked in z**, with heights | reflector / fuel / reflector |
| core | an xy map of assemblies, full or quarter | 193 assemblies |

Beside the layers are the **states**: fuel temperature, moderator
temperature and density, and soluble boron. A material stores only *what*
it is made of; its temperature, and for the coolant its density and boron
content, come from the state. One material can therefore be evaluated at any
state point. This is what the branch calculations of a nodal model need:
the same lattice is computed at many fuel temperatures, water densities and
boron concentrations.

## Slices and stacks

The core is three-dimensional, but the file builds it from 2D slices
stacked in z, not from general 3D shapes. The reason is the calculation
chain:

* A **lattice calculation** (the MOC of the B steps) solves one xy slice
  that is assumed infinitely tall. Each distinct lattice, at each state,
  gives one set of homogenized group constants.
* A **nodal calculation** divides the core into boxes, one assembly wide
  and one axial zone tall. Each box takes the constants of the lattice in
  that zone.
* The **axial stack** of an assembly links the two: it says which lattice
  sits at which height.

Anything that is not uniform in z becomes its own axial zone: end
reflectors, plenum, spacer grids, axial blankets, and lattices with and
without inserted control rods. Assemblies may have different zone
boundaries; the nodal mesh is cut at the union of all of them. A control
rod whose tip lies inside a node is the classic difficulty. Simply
averaging rodded and unrodded constants by volume ignores the flux
depression in the rodded part and overestimates the rod's worth; this
error is called **rod cusping**. Nodal codes correct it by weighting with
the neighbouring fluxes or by refining the mesh axially.

## Number densities from what a specification gives

Specifications give mass densities and compositions; transport needs atoms
per barn-cm. For a compound of mass density $\rho$ whose molecule contains
$n_e$ atoms of element $e$,

$$
N_e = \frac{\rho\, N_A}{M}\, n_e \times 10^{-24},
\qquad
M = \sum_e n_e\, \bar m_e,
$$

with $N_A$ Avogadro's number, $M$ the molar mass of the molecule and
$10^{-24}$ cm²/barn. The mean atomic mass $\bar m_e$ of an element follows
from its isotopics. By atom fractions $a_i$: $\bar m_e = \sum_i a_i m_i$.
By weight fractions $w_i$: $1/\bar m_e = \sum_i w_i/m_i$, and the atom
fractions are $a_i = \bar m_e w_i / m_i$. Enrichments are usually stated
by weight, natural boron by atom, so the file accepts both.

**Soluble boron** is given in ppm by weight of the solution: a fraction
$w_B = \text{ppm}\times 10^{-6}$ of the mass is boron, $1-w_B$ is water.
This is the convention of the VERA specification, which inputs 1300 ppm as
weight fractions 0.0013 and 0.9987 [[VERA](#ref-vera), Problem 1,
"Materials"].

**Atomic masses** come from the data already in use: every ENDF evaluation
stores the mass ratio AWR of the nuclide to the neutron (MF1/MT451), so
$m_i = \text{AWR}\cdot m_n$ with the CODATA 2018 neutron mass
$m_n = 1.00866491595$ u. For H-1 this gives 1.007825 u and for U-238
238.0508 u.

**Check against VERA.** The VERA specification gives the fuel as UO$_2$ of
10.257 g/cm³ with the uranium isotopics of its Table P1-3 by weight. It
gives the moderator as water of 0.743 g/cm³ (1A) or 0.661 g/cm³ (1B–1D),
with 1300 ppm boron at 19.9 at% B-10. It also lists the resulting number
densities (Table P1-4). Computed from the first description, they agree
with the second:

| Material | Nuclide | computed (/b-cm) | Table P1-4 | difference |
| --- | --- | --- | --- | --- |
| fuel | U-235 | 7.18131e-4 | 7.18132e-4 | −1.3e-6 |
| fuel | U-238 | 2.21546e-2 | 2.21546e-2 | −1.5e-6 |
| fuel | O-16 | 4.57642e-2 | 4.57642e-2 | +0.7e-6 |
| water 1A | H-1 | 4.96224e-2 | 4.96224e-2 | −0.7e-6 |
| water 1A | B-10 | 1.07071e-5 | 1.07070e-5 | +5.4e-6 |
| water 1B | H-1 | 4.41459e-2 | 4.41459e-2 | −0.7e-6 |
| water 1B | B-10 | 9.52539e-6 | 9.52537e-6 | +2.1e-6 |

The table prints six digits. The differences are at that level,
consistent with rounding and with the specification's atomic masses (from
NIST) differing slightly from the ENDF ones. The tests check every nuclide
to $10^{-5}$. The gap helium and the Zircaloy are given only as number
densities, so the file states them that way.

## Homogenizing a slice

A lattice is infinitely tall, so volume fractions are area fractions. A
ring between radii $r_{i-1}$ and $r_i$ has area $\pi(r_i^2 - r_{i-1}^2)$,
and the last material fills the cell up to its edge, $p^2 - \pi r_n^2$ for
pitch $p$. The homogenized density of a nuclide is $\sum_r N_r A_r / A$,
which conserves the number of atoms in the slice. For VERA's cell the fuel
covers 33.2% of the area, the gap 1.4%, the clad 10.1% and the water 55.4%.
These are the components that chapter 12 used.

## The file

The model is written in TOML, a configuration format that Python reads with
its standard library (`tomllib`). Here is an abridged version of
`models/vera-p1.toml`:

```toml
[model]
format = 1
name = "VERA problem 1"

[materials.fuel]                    # UO2 (Tables P1-2, P1-3)
temperature = "fuel"
density = 10.257
formula = { U = 1, O = 2 }
weight_percent.U = { U234 = 0.0263, U235 = 3.1, U236 = 0.0143, U238 = 96.8594 }
atom_percent.O = { O16 = 100.0 }

[materials.gap]                     # explicit number densities (Table P1-4)
temperature = "moderator"
atom_densities = { He4 = 2.68714e-05 }

[materials.moderator]               # density and boron from the state
temperature = "moderator"
density = "moderator"
formula = { H = 2, O = 1 }
atom_percent.H = { H1 = 100.0 }
atom_percent.O = { O16 = 100.0 }
soluble_boron.atom_percent.B = { B10 = 19.9, B11 = 80.1 }
thermal = { H1 = "HinH2O" }         # S(alpha, beta), chapter 11

[pins.fuel_rod]
radii = [0.4096, 0.418, 0.475]
materials = ["fuel", "gap", "clad", "moderator"]

[lattices.cell]
pitch = 1.26
legend = { F = "fuel_rod" }
map = "F"

[states.1A]
fuel_temperature = 565.0
moderator_temperature = 565.0
moderator_density = 0.743
boron_ppm = 1300
reference = { keff = 1.187038, sigma = 0.000054, source = "Table P1-5" }
```

Larger maps are multi-line strings with one legend letter per pin
(`F F G F F ...`). Assemblies list their axial segments from bottom to
top, and the core map uses `.` for positions outside the core.

**Checks.** Reading the file checks everything, and every error names its
place in the file:

* unknown keys, which catch typos;
* references to undefined materials, pins, lattices or assemblies;
* radii that do not increase, or rings that do not fit in the pitch;
* ragged or non-square maps;
* isotopic percentages that do not add up to 100, or an isotope listed
  under the wrong element;
* assemblies of different heights in one core;
* lattices wider than the assembly pitch.

The choices that are part of the format but not implemented yet are
rejected explicitly:

* `geometry = "hex"`;
* `symmetry = "eighth"`.

Adding them later extends the code without changing existing files.

## What the format leaves out

* **Hexagonal lattices** (VVER, fast reactors, prismatic HTGRs) use the same
  layers, with hex-z nodal methods in the core. The format reserves
  `geometry = "hex"` for them.
* **Strong heterogeneity**: MOX next to UO$_2$, small or leaky cores, strong
  absorbers. Here the problem is how group constants are made, not how the
  core is described. The usual remedies are discontinuity factors,
  calculations on a group of neighbouring assemblies, and constants from 3D
  Monte Carlo.
* **Geometry that is not extruded** (pebble beds, conical or spherical
  regions) does not fit slices and stacks. It needs regions described by a
  3D mesh or by constructive solid geometry, finite-element diffusion, and
  constants computed region by region. That would be a different kind of
  model, outside this project's scope.
* **Branch matrices and burnup histories** come with the next steps: M2
  runs a lattice calculation for each state of a branch matrix, and M3
  tabulates the resulting constants for the nodal code.

## Further reading

* [[VERA](#ref-vera)] (free): the benchmark specification and the VERAIn
  common input; its §2 gives the material definitions used here.
* [[ARMI](#ref-armi)] (free): an open-source reactor model framework built
  around one shared description.
* Code: [`spec.py`](../../src/fission/model/spec.py) (the layers),
  [`load.py`](../../src/fission/model/load.py) (reading and checks),
  [`materials.py`](../../src/fission/model/materials.py) (number densities),
  [`homogenize.py`](../../src/fission/model/homogenize.py) (area weighting),
  [`vera-p1.toml`](../../models/vera-p1.toml); `make model MODEL=vera-p1
  STATE=1A`.
