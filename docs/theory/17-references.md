# References {.unnumbered}

Links were checked when these notes were written (November 2026). **(free)**
means the full text is freely available at the link; **(not free)** means a
paywalled paper or a commercial book, cited because it is the original
source or a standard text.

## Data formats, libraries and processing {.unnumbered}

- []{#ref-endf102}**[ENDF-102]** D. A. Brown (ed.), *ENDF-6 Formats Manual:
  Data Formats and Procedures for the Evaluated Nuclear Data Files ENDF/B-VI,
  ENDF/B-VII and ENDF/B-VIII*, CSEWG Document ENDF-102, Report
  BNL-224854-2023-INRE, Brookhaven National Laboratory (2023). The definition
  of the format; Appendix D contains all resonance formulas.
  <https://www.nndc.bnl.gov/endfdocs/ENDF-102-2023.pdf> **(free)**
- []{#ref-endfb81}**[ENDF/B-VIII.1]** G. P. A. Nobre et al., "ENDF/B-VIII.1:
  Updated Nuclear Reaction Data Library for Science and Applications",
  *Nuclear Data Sheets* **210**, 1–224 (2026). Preprint:
  <https://arxiv.org/abs/2511.03564> **(free)**. Library download:
  <https://www.nndc.bnl.gov/endf-releases/?version=B-VIII.1> **(free)**
- []{#ref-endfb80}**[ENDF/B-VIII.0]** D. A. Brown et al., "ENDF/B-VIII.0: The
  8th Major Release of the Nuclear Reaction Data Library", *Nuclear Data
  Sheets* **148**, 1–142 (2018). <https://doi.org/10.1016/j.nds.2018.02.001>
- []{#ref-njoy}**[NJOY]** R. E. MacFarlane (original author), with
  D. W. Muir, R. M. Boicourt, A. C. Kahler, J. L. Conlin, W. Haeck, *The NJOY
  Nuclear Data Processing System, Version 2016*, LA-UR-17-20093, Los Alamos
  National Laboratory (2016; updated 2019 for NJOY2016.53).
  Ch. 3 RECONR (reconstruction), ch. 4 BROADR (Doppler), ch. 5 UNRESR,
  ch. 7 THERMR, ch. 23 PURR. <https://github.com/njoy/NJOY2016-manual>
  **(free)**
- []{#ref-sammy}**[SAMMY]** N. M. Larson, *Updated Users' Guide for SAMMY:
  Multilevel R-Matrix Fits to Neutron Data Using Bayes' Equations*,
  ORNL/TM-9179/R8, ENDF-364/R2, Oak Ridge National Laboratory (2008).
  §II: R-matrix theory and Reich–Moore in practice.
  <https://www.osti.gov/biblio/941054> **(free)**
- []{#ref-openmc}**[OpenMC]** OpenMC documentation, Methods, ch. 3 "Cross
  Section Representations".
  <https://docs.openmc.org/en/stable/methods/cross_sections.html> **(free)**

## Resonance theory {.unnumbered}

- []{#ref-frohner}**[Fröhner]** F. H. Fröhner, *Evaluation and Analysis of
  Nuclear Resonance Data*, JEFF Report 18, OECD Nuclear Energy Agency (2000).
  §3: resonance theory (R-matrix, SLBW, MLBW, Reich–Moore, Doppler
  broadening); §4: statistical theory, unresolved range, Bondarenko factors.
  Catalogue page: <https://oecd-nea.org/jcms/pl_13418>; full text:
  <https://physics.mff.cuni.cz/kfnt/vyuka/statisticke_metody/Frohner.pdf>
  **(free)**
- []{#ref-lanethomas}**[Lane & Thomas]** A. M. Lane, R. G. Thomas,
  "R-Matrix Theory of Nuclear Reactions", *Reviews of Modern Physics* **30**,
  257–353 (1958). The standard reference for R-matrix theory.
  <https://doi.org/10.1103/RevModPhys.30.257> **(not free)**
- []{#ref-wignereisenbud}**[Wigner & Eisenbud]** E. P. Wigner, L. Eisenbud,
  "Higher Angular Momenta and Long Range Interaction in Resonance Reactions",
  *Physical Review* **72**, 29–41 (1947).
  <https://doi.org/10.1103/PhysRev.72.29> **(not free)**
- []{#ref-breitwigner}**[Breit & Wigner]** G. Breit, E. Wigner, "Capture of
  Slow Neutrons", *Physical Review* **49**, 519–531 (1936).
  <https://doi.org/10.1103/PhysRev.49.519> **(not free)**
- []{#ref-reichmoore}**[Reich & Moore]** C. W. Reich, M. S. Moore,
  "Multilevel Formula for the Fission Process", *Physical Review* **111**,
  929–933 (1958). <https://doi.org/10.1103/PhysRev.111.929> **(not free)**
- []{#ref-porterthomas}**[Porter & Thomas]** C. E. Porter, R. G. Thomas,
  "Fluctuations of Nuclear Reaction Widths", *Physical Review* **104**,
  483–491 (1956). <https://doi.org/10.1103/PhysRev.104.483> **(not free)**
- []{#ref-blattweisskopf}**[Blatt & Weisskopf]** J. M. Blatt,
  V. F. Weisskopf, *Theoretical Nuclear Physics*, Wiley (1952); Dover
  reprint (1991). Ch. VIII–X: compound nucleus, resonances. **(not free)**
- []{#ref-mughabghab}**[Mughabghab]** S. F. Mughabghab, *Atlas of Neutron
  Resonances: Resonance Properties and Thermal Cross Sections Z = 1–100*,
  6th ed., Elsevier (2018). Source of the thermal values and resonance
  parameters used as checks. **(not free)**
- []{#ref-kalbach}**[Kalbach]** C. Kalbach, "Systematics of Continuum
  Angular Distributions: Extensions to Higher Energies", *Physical Review C*
  **37**, 2350 (1988). The angular systematics used in ENDF MF6 (LANG = 2).
  <https://doi.org/10.1103/PhysRevC.37.2350> **(not free)**
- []{#ref-cullen}**[Cullen & Weisbin]** D. E. Cullen, C. R. Weisbin, "Exact
  Doppler Broadening of Tabulated Cross Sections", *Nuclear Science and
  Engineering* **60**, 199–229 (1976). <https://doi.org/10.13182/NSE76-1>
  **(not free)**
- []{#ref-solbrig}**[Solbrig]** A. W. Solbrig, Jr., "Doppler Broadening of
  Low-Energy Neutron Resonances", *Nuclear Science and Engineering* **10**,
  167–168 (1961). The exact free-gas kernel.
  <https://doi.org/10.13182/NSE61-A25954> **(not free)**
- []{#ref-lamb}**[Lamb]** W. E. Lamb, "Capture of Neutrons by Atoms in a
  Crystal", *Physical Review* **55**, 190–197 (1939). Lattice binding; the
  effective temperature.
  <https://doi.org/10.1103/PhysRev.55.190> **(not free)**

## Benchmarks {.unnumbered}

- []{#ref-vera}**[VERA]** A. T. Godfrey, *VERA Core Physics Benchmark
  Progression Problem Specifications*, CASL-U-2012-0131-004, Revision 4,
  Oak Ridge National Laboratory (2014). Problem 1: geometry (Table 1), states
  (Table P1-1), number densities (Table P1-4), CE KENO-VI reference
  eigenvalues (Table P1-5). <https://corephysics.com/docs/CASL-U-2012-0131-004.pdf>
  **(free)**. Inputs: <https://github.com/CASL/VERAin> **(free)**

- []{#ref-c5g7}**[C5G7]** OECD/NEA Expert Group on 3-D Radiation Transport
  Benchmarks (reference MCNP solution: M. A. Smith, N. Tsoulfanidis,
  R. N. Blomquist, E. E. Lewis, Argonne), *Benchmark on Deterministic
  Transport Calculations Without Spatial Homogenisation: A 2-D/3-D MOX Fuel
  Assembly Benchmark*, NEA/NSC/DOC(2003)16, OECD/NEA (2003). Geometry
  (Chapter 2, Figures 1–3), seven-group cross sections (Appendix A,
  Table 1), reference MCNP results (Tables 3–5).
  <https://www.oecd-nea.org/upload/docs/application/pdf/2019-12/nsc-doc2003-16.pdf>
  **(free)**

## Software {.unnumbered}

- []{#ref-armi}**[ARMI]** TerraPower, *ARMI: the Advanced Reactor Modeling
  Interface*, an open-source reactor analysis automation framework
  (Apache-2.0). <https://github.com/terrapower/armi>,
  documentation <https://terrapower.github.io/armi/> **(free)**

- []{#ref-openmoc}**[OpenMOC]** W. Boyd, S. Shaner, L. Li, B. Forget,
  K. Smith, "The OpenMOC method of characteristics neutral particle
  transport code", *Annals of Nuclear Energy* **68**, 43–52 (2014).
  <https://doi.org/10.1016/j.anucene.2013.12.012> **(not free)**
- []{#ref-openmocdocs}**[OpenMOC docs]** OpenMOC theory manual: §5 track
  generation (cyclic tracks, azimuthal weights),
  <https://mit-crpg.github.io/OpenMOC/methods/track_generation.html>; the
  method of characteristics,
  <https://mit-crpg.github.io/OpenMOC/methods/method_of_characteristics.html>
  **(free)**

## Reactor physics {.unnumbered}

- []{#ref-bell}**[Bell & Glasstone]** G. I. Bell, S. Glasstone, *Nuclear
  Reactor Theory*, Van Nostrand Reinhold (1970). Ch. 1: the transport
  equation; ch. 7: thermalization; ch. 8: resonance absorption (SLBW, Doppler
  broadening, unresolved parameters, narrow-resonance approximation).
  <https://www.osti.gov/biblio/4074688>
  **(free)**
- []{#ref-mit2205}**[MIT 22.05]** MIT OpenCourseWare, *22.05 Neutron Science
  and Reactor Physics* (Fall 2009).
  <https://ocw.mit.edu/courses/22-05-neutron-science-and-reactor-physics-fall-2009/>
  **(free)**; see also *22.106 Neutron Interactions and Applications*
  (Spring 2010),
  <https://ocw.mit.edu/courses/22-106-neutron-interactions-and-applications-spring-2010/>
  **(free)**
- []{#ref-lamarsh}**[Lamarsh]** J. R. Lamarsh, *Introduction to Nuclear
  Reactor Theory*, Addison-Wesley (1966). **(not free)**
- []{#ref-dh}**[Duderstadt & Hamilton]** J. J. Duderstadt, L. J. Hamilton,
  *Nuclear Reactor Analysis*, Wiley (1976). **(not free)**
- []{#ref-hebert}**[Hébert]** A. Hébert, *Applied Reactor Physics*, 3rd ed.,
  Presses internationales Polytechnique (2020). Lattice-code methods
  including resonance self-shielding. **(not free)**

## Mathematics and constants {.unnumbered}

- []{#ref-dlmf}**[DLMF]** *NIST Digital Library of Mathematical Functions*,
  §10.47, spherical Bessel functions. <https://dlmf.nist.gov/10.47>
  **(free)**
- []{#ref-codata}**[CODATA]** NIST, *Fundamental Physical Constants*
  (CODATA 2018 values used in the code).
  <https://physics.nist.gov/cuu/Constants/> **(free)**
