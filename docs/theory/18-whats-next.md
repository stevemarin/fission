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
corrects it (SPH factors, rings in the pellet). This chapter previews what
remains.

## Beyond chapter 17

The pin cell is now within 160–350 pcm of OpenMC on the same data, and
the error sits in U-238 capture in resonances with large neutron widths.
In rough order of importance:

* **The clad in the slowing-down calculation.** Fuel, gap and clad should
  be marched together, each with its own composition, coupled by
  first-collision probabilities computed at every energy instead of
  tabulated against one cross section. This is the pointwise
  collision-probability method, and it is also the cross-check of MOC
  planned for step C4.
* **Few groups for the core.** The pin and then the 17×17 assembly are
  homogenized and condensed to two groups for nodal diffusion (C3). SPH
  factors appear there again, in their original role: making the
  homogenized assembly reproduce the heterogeneous reaction rates.
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
