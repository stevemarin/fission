# Preface {.unnumbered}

These notes explain the theory behind the first part of a lattice physics
code written from scratch: turning *evaluated nuclear data* (the ENDF/B-VIII.1
library) into *pointwise cross sections* $\sigma(E)$ that a transport solver can
use. Most of the space goes to the physics and mathematics of neutron
resonances. The software appears only where it helps to see a formula at work;
each chapter ends with pointers into the code and to further reading.

## Where this fits {.unnumbered}

A reactor core calculation is usually done in two steps:

1. **Lattice physics.** A single fuel pin or assembly is solved with neutron
   transport, using detailed energy-dependent cross sections, to produce a few
   *homogenized, energy-averaged* cross sections.
2. **Core analysis.** Those few-group constants feed a fast diffusion solver
   for the whole three-dimensional core.

Step 1 needs cross sections first. They come from *evaluations*: files in which
specialists (evaluators) have condensed experiments and nuclear models into
parameters. These notes cover the path from such a file to a table of
$\sigma(E)$:

$$
\text{ENDF file} \;\xrightarrow{\text{read}}\; \text{parameters}
\;\xrightarrow{\text{resonance theory}}\; \sigma(E)
\;\xrightarrow{\text{adaptive grid}}\; \text{linear table}
$$

The tables are first built at **0 K** (no thermal motion of the target
nuclei), then Doppler-broadened to the temperatures of a reactor. The
unresolved resonance range, self-shielding and thermal scattering come next;
chapter 16 treats self-shielding in a pin cell, chapter 17 checks it against
continuous-energy Monte Carlo and corrects it, and chapter 18 previews the rest.

## How to read these notes {.unnumbered}

* Chapters 1 and 2 are background: what a cross section is, and how ENDF
  stores one.
* Chapter 3 introduces the physics shared by all resonance formulas.
* Chapters 4–6 build the three resonance formalisms used by the nuclides in a
  PWR fuel pin, from the simplest (SLBW) to the most accurate (Reich–Moore).
  Each one fixes a defect of the previous one.
* Chapter 7 turns formulas into tables; chapter 8 heats them up.
* Chapter 9 turns to what reactions emit: neutrons per fission, spectra and
  angular distributions.
* Chapter 10 condenses all of it into multigroup constants for transport,
  and chapter 11 adds thermal scattering. Chapter 12 puts them to use: the
  spectrum and k-infinity of a homogenized PWR pin cell.
* Chapter 13 describes the model file: how a reactor is written down once,
  as data, from materials to the core map. Chapter 14 builds the geometry of
  the transport solver: flat-source regions, cyclic tracks and ray tracing;
  chapter 15 solves the transport equation on them (MOC) and checks the
  result on the C5G7 benchmark. Chapter 16 joins the two tracks: resonance
  self-shielding of a heterogeneous pin cell, solved with MOC on the
  libraries of chapters 10 and 11. Chapter 17 compares that pin cell with a
  continuous-energy Monte Carlo calculation on the same data, separates data
  from method, and fixes the two largest method errors it finds
  (superhomogenization and spatially dependent shielding).
* Sections marked **(derivation)** can be skipped on a first reading; the
  results they derive are stated before them.

Numbers quoted as "computed" were produced by this project's code from
ENDF/B-VIII.1; they are compared with independent published values wherever
possible.

## Notation and units {.unnumbered}

The notation follows the ENDF-6 formats manual [[ENDF-102](#ref-endf102)], so
the manual can be read alongside. Units are the ones used in ENDF files:

| Quantity | Symbol | Unit |
| --- | --- | --- |
| neutron energy (laboratory frame) | $E$ | eV |
| cross section | $\sigma$ | barn, $1\,\text{b} = 10^{-24}\,\text{cm}^2$ |
| radii (channel, scattering) | $a$, $AP$ | $10^{-12}$ cm $= 10$ fm |
| wave number | $k$ | $(10^{-12}\,\text{cm})^{-1}$ |
| resonance energy | $E_r$ | eV |
| partial widths (neutron, gamma, fission) | $\Gamma_n, \Gamma_\gamma, \Gamma_f$ | eV |
| target spin, resonance spin | $I$, $J$ | $\hbar$ |
| orbital angular momentum | $l$ | $\hbar$ |
| ratio of target mass to neutron mass | AWR, AWRI | – |

With radii in $10^{-12}$ cm and $k$ in $(10^{-12}\,\text{cm})^{-1}$, the
combination $\pi/k^2$ comes out directly in barns, because
$(10^{-12}\,\text{cm})^2 = 10^{-24}\,\text{cm}^2$.

## Building these notes {.unnumbered}

The chapters are plain Markdown with LaTeX mathematics. To build a single
self-contained HTML file (no network needed to view it):

```sh
make docs          # writes docs/build/theory.html
make figures       # regenerates the figures from ENDF data (needs make data)
```
