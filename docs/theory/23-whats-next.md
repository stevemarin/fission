# What comes next

Chapters 7 and 8 give pointwise tables at any temperature, with averaged
("infinitely dilute") cross sections in the unresolved range; chapter 9 adds
the data on the neutrons that reactions emit, chapter 10 condenses
everything into infinitely dilute multigroup constants, chapter 11 adds
thermal scattering, chapter 12 shields them in a homogeneous mixture, and
chapter 13 describes the reactor model they are applied to; chapter 14
builds the geometry of the transport solver, chapter 15 the solver itself,
chapter 16 shields the cross sections of a real pin cell, and chapter 17
checks that pin against continuous-energy Monte Carlo on the same data and
corrects it (SPH factors, rings in the pellet), chapter 18 condenses it
into two-group constants for a nodal code, and chapter 19 adds a second
transport method (collision probabilities), the clad to the slowing-down
calculation, and the 17×17 assembly, chapter 20 repeats the lattice
calculation at every point of a state matrix and builds the tables a nodal
code reads, chapter 21 solves the core with nodal diffusion, and
chapter 22 brings it to power with thermal-hydraulic feedback. This chapter previews what
remains.

## Beyond chapter 22

* **The assembly at converged tracking.** The pin cell's defaults (256
  directions, 10 rings) are too expensive for 17×17. Rod-by-rod tracking
  (modular rays) or a coarser azimuthal set with more rings would bring the
  assembly to the same accuracy.
* **The slowing-down march with moving hydrogen.** The march still treats
  hydrogen at rest above 1 eV. U-238 capture at 4–10 eV is 0.5% low, and at
  25–100 eV 0.5–0.7% high.
* **Rod-by-rod shielding.** Every rod of the assembly shares one shielded
  library. Rods next to guide tubes and at the edge see more water. Their
  own Dancoff factors would refine the pin powers.
* **The state matrix with MOC.** Chapter 20 ran the matrix with the
  homogenized cell, whose Doppler feedback is 55% too strong. The same
  matrix with the heterogeneous solver takes about three hours.
* **Burning the absorbers.** Chapter 22's Pyrex rods are as loaded. Their
  B-10 burns during the cycle, which needs depletion. Their worth is
  1.1–1.2% weaker than KENO's, and the rod gets no shielding of its own.
* **Control rods, further.** Chapter 22's bank worths agree with KENO
  within −7..+5% (total −1.7%). What is left out: a boron axis in the
  tables, so that critical boron and the 1170 ppm of the rod-worth
  measurements can be computed directly; a full-core (or rotational)
  quarter, for banks SC and SD; and face-dependent discontinuity factors.
  The reflector's factor belongs to its fuel face, but NEM applies it on
  every face of the node, the outer vacuum face too. The axial mesh is
  cut at every rod tip, so a partly inserted rod causes no cusping, at
  the price of thin nodes.
* **The remaining resonance formats.** LRF = 7 with an inelastic channel
  (Fe-57), Reich–Moore with negative channel spins (Mn-55, Cr-53) and the
  LRU = 2, LRF = 1 unresolved format (Fe-58) are still missing from the
  steel and the clad. The LRF = 7 kernel is Python only. A Rust port pays
  off once a full-range table needs rebuilding often.
* NEM runs in seconds, the fine-mesh finite
  difference in minutes; its red-black sweep is the hot path to port.
* **Depletion.** Ring-wise shielding matters most once plutonium builds up
  at the pellet rim.

## The unresolved range

Between ~20 keV and ~150 keV (U-238), only average resonance parameters are
known: mean level spacing, average widths, and the number of degrees of
freedom of their Porter–Thomas-type distributions [[ENDF-102](#ref-endf102)
§D.2]. Averaged cross sections follow from statistical (Hauser–Feshbach type)
formulas [[Fröhner](#ref-frohner) §4.2]. Self-shielding there needs either
Bondarenko factors (NJOY's UNRESR, [[NJOY](#ref-njoy) ch. 5]) or
**probability tables**, distributions of the cross section sampled from
statistical resonance ladders (NJOY's PURR, ch. 23). The unresolved range
contributes little to U-238's dilute resonance integral (0.9 b of 275 b);
in the VERA pin its shielding is worth +20 ± 20 pcm (chapter 17). It
matters more for fast spectra.

## Further reading

* [[NJOY](#ref-njoy)] ch. 5 UNRESR, ch. 23 PURR (free).
* [[Fröhner](#ref-frohner)] §4 (unresolved range, Bondarenko factors)
  (free).
* [[Bell & Glasstone](#ref-bell)] ch. 8 (resonance absorption: unresolved parameters §8.2, the narrow-resonance
  approximation §8.3, heterogeneous systems §8.4) (free).
* [[OpenMC](#ref-openmc)] methods documentation: how a modern Monte Carlo
  code represents cross sections (free).
* [[Hébert](#ref-hebert)], resonance self-shielding in lattice codes (not
  free).
