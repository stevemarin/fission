# Few-group constants

The pin-cell calculation of chapters 15–17 has 56 flat-source regions and
172 energy groups. A core has some 50,000 fuel rods, and a nodal diffusion
code describes it with one homogeneous node per assembly (or per quarter
assembly) and two energy groups. This chapter turns the lattice solution
into constants for such a node, in three steps:

1. **homogenize** the cell, keeping the 172 groups;
2. replace the infinite-lattice spectrum by the **B1 critical spectrum** of
   a leaking reactor;
3. **condense** to two groups, with the cut at 0.625 eV.

Each step preserves reaction rates by construction, and the tests check that
it does to rounding error. Here a "cell" is one pin; a 17×17 assembly works
the same way (step C4).

## Homogenization

Every reaction rate of the heterogeneous cell must be reproduced by one
homogeneous material with the same volume and the same average flux. For
reaction $x$ in fine group $g$,

$$
\Sigma_{x,g} = \frac{\sum_r V_r\,\Sigma_{x,r,g}\,\phi_{r,g}}{\sum_r V_r\,\phi_{r,g}},
$$

summed over the MOC regions $r$. This is **flux-volume weighting**. The
scattering matrices are weighted with the flux of the departure group,
including the P1 matrix. The net current is zero in an infinite lattice, so
nothing better is available.

$\chi$ is weighted with each region's fission production, so the
homogeneous cell emits the same fission neutrons. With these constants, the
infinite homogeneous medium has the MOC flux spectrum and the MOC
eigenvalue. For VERA 1A they agree to $4\times10^{-10}$.

The SPH factors of chapter 17 are already inside the region cross sections,
so the homogenized rates are the corrected ones.

Two remarks on what this does not preserve:

* **Surface fluxes and currents.** The homogeneous node has a flat flux,
  whereas the real flux dips in the fuel and peaks in the water. When
  neighbouring assemblies differ, the nodal code needs **discontinuity
  factors**, the ratio of the heterogeneous surface flux to the homogeneous
  one, to reproduce the interface currents. In a lattice of identical cells
  every interface sees the same factor on both sides, and it cancels. The
  constants written here therefore carry factors of 1. They will become
  meaningful with the assembly calculation (C4).
* **The transport correction.** MOC used the outflow-corrected total
  (chapter 16). The homogenization uses the plain total and both
  scattering moments, because the next step treats anisotropic scattering
  itself.

## The B1 fundamental mode

The lattice calculation assumes an infinite array: nothing leaks. In a core
the same assembly does leak, and leakage changes the spectrum. Fast
neutrons travel farther before they are absorbed, so they leak more, and
the spectrum seen by an assembly in a finite core is harder than the
infinite-lattice one. The standard correction is to solve the homogenized
medium in the **fundamental mode** of a large bare reactor [[Hébert](#ref-hebert)],

$$
\phi(x, E, \mu) = \psi(E, \mu)\, e^{iBx},
$$

and to choose the buckling $B^2$ that makes it exactly critical.

**Derivation.** Substituting into the one-dimensional transport equation
with linearly anisotropic scattering, the angular dependence separates.
Expanding $\psi$ in its Legendre moments, flux $\phi_g$ and current
$iJ_g$, and integrating the angular equation exactly instead of truncating
it at $P_1$ gives the **B1 equations**:

$$
\Sigma_{t,g}\phi_g + B J_g = \sum_{g'} \Sigma_{s0,g'\to g}\,\phi_{g'} + \frac{\chi_g}{k} \sum_{g'} \nu\Sigma_{f,g'}\phi_{g'},
$$

$$
3\,\alpha_g\,\Sigma_{t,g} J_g - B\phi_g = 3 \sum_{g'} \Sigma_{s1,g'\to g}\, J_{g'},
$$

with, for $x = B/\Sigma_{t,g}$,

$$
\alpha(x) = \frac{x^2 \arctan x}{3\,(x - \arctan x)} \quad (B^2>0), \qquad
\alpha(x) = \frac{x^2 \operatorname{artanh} x}{3\,(\operatorname{artanh} x - x)} \quad (B^2<0,\ x = |B|/\Sigma_{t,g}).
$$

$\alpha$ tends to 1 as $B \to 0$ (its series is $1 \pm 4x^2/15$). With
$\alpha = 1$ these are the $P_1$ equations, and for a diagonal P1 matrix
the second equation gives Fick's law, $J_g = D_g B\phi_g$ with
$D_g = 1/(3(\Sigma_{t,g} - \Sigma_{s1,g\to g}))$. **B1 computes the
diffusion coefficients**, $D_g = J_g/(B\phi_g)$, including the anisotropic
scattering between groups. For VERA 1A at the critical buckling, $B/\Sigma_t$
reaches 0.46 in the highest group (17–20 MeV, $\Sigma_t$ = 0.12/cm), where
$\alpha = 1.054$. It is much smaller wherever $\Sigma_t$ is larger, so the
B1 correction matters only for the fastest neutrons, which are also the
ones that leak.

**Solution** (`fewgroup/b1.py`). Eliminate the currents,
$J = B M\phi$ with $M = (3\,\mathrm{diag}(\alpha\Sigma_t) - 3S_1^T)^{-1}$.
That leaves $(\mathrm{diag}(\Sigma_t) - S_0^T + B^2 M)\,\phi = \chi F/k$, so
$k(B^2) = \nu\Sigma_f\cdot A(B^2)^{-1}\chi$: one dense 172×172 solve per
buckling. $k$ decreases with $B^2$. A secant search for $k = 1$ starts from
the one-group estimate $B^2 = (k_\infty-1)/M^2$ and converges in a few
steps.

## Condensation

With the critical spectrum $\phi_g$ and currents $J_g$, every few-group
constant is the rate-preserving average over its fine groups $g \in G$:

$$
\phi_G = \sum_{g\in G}\phi_g,\qquad
\Sigma_{x,G} = \frac{\sum_{g\in G}\Sigma_{x,g}\phi_g}{\phi_G},\qquad
\Sigma_{s,G\to G'} = \frac{\sum_{g\in G}\sum_{g'\in G'}\Sigma_{s0,g\to g'}\phi_g}{\phi_G},
$$

$$
\chi_G = \sum_{g\in G}\chi_g,\qquad
D_G = \frac{\sum_{g\in G} J_g/B}{\phi_G}.
$$

$D_G$ is weighted by the current, not the flux. Then the few-group leakage
$D_G B^2\phi_G$ equals the B1 leakage $\sum_g BJ_g$ exactly, and the
two-group diffusion equation is critical at the critical buckling. The test
suite checks this, and the converse: condensed with the infinite spectrum,
the two-group $k_\infty$ is the 172-group one.

$\Sigma_a$ is defined as total minus out-scattering, so that the nodal
code's removal $\Sigma_a + \Sigma_{s,1\to2}$ balances exactly. $(n,2n)$
counts there as negative absorption, as in MOC. $\kappa\Sigma_f$ uses
the **recoverable energy per fission** $E_R = E_T - E_\nu$ (all fission energy except the neutrinos'), from MF1/MT458
[[ENDF-102](#ref-endf102) §1.5]. ENDF/B-VIII.1 gives thermal values: 194.08
MeV for U-235 and 196.54 MeV for U-238. For U-235 and U-238 it also
tabulates the energy dependence of four components, which raises $E_R$ by
about 1 MeV at 1 MeV incident energy (`fewgroup/kappa.py`).

## Results

`make fewgroup STATES=1A PIN_ARGS="--rings --sph"` runs chapter 17's pin
cell, then all three steps (30 s). It writes the constants in the form the
nodal code reads to `endf-data/fewgroup/vera-p1-1A.json`.

| 1A | fast (> 0.625 eV) | thermal |
| --- | --- | --- |
| $D$ (cm) | 1.3517 | 0.4256 |
| $\Sigma_a$ (1/cm) | 0.010125 | 0.10037 |
| $\nu\Sigma_f$ (1/cm) | 0.0072015 | 0.14968 |
| $\kappa\Sigma_f$ (MeV/cm) | 0.5541 | 11.953 |
| $\Sigma_{s,1\to2}$, $\Sigma_{s,2\to1}$ (1/cm) | 0.015806 | 0.0016394 |
| flux fraction | 0.867 | 0.133 |

The critical buckling is $B^2 = 3.256\times10^{-3}$ cm⁻². The implied
migration area $M^2 = (k_\infty-1)/B^2 = 56.9$ cm² is typical of a PWR.
Leakage hardens the spectrum. The fast-to-thermal flux ratio rises from
6.19 (infinite lattice) to 6.54, and the two-group $k_\infty$ of the
B1-weighted constants is 1.18364, 154 pcm below the lattice's 1.18518.
That is the cost of weighting with a leaking spectrum: in a real core the
leakage is what removes the excess.

The states respond as they should:

| | 1A | 1B | 1C | 1D |
| --- | --- | --- | --- | --- |
| $B^2$ (10⁻³ cm⁻²) | 3.256 | 2.751 | 2.593 | 2.458 |
| $M^2$ (cm²) | 56.9 | 65.4 | 65.2 | 65.0 |
| $D_1$ (cm) | 1.3517 | 1.4216 | 1.4225 | 1.4232 |
| $\Sigma_{a,1}$ (1/cm) | 0.010125 | 0.009938 | 0.010122 | 0.010277 |
| $\Sigma_{s,1\to2}$ (1/cm) | 0.015806 | 0.013765 | 0.013711 | 0.013668 |

* **1A → 1B: less dense water.** The water density falls from 0.743 to
  0.661 g/cm³. Neutrons travel farther ($D_1$ +5%) and slow down less
  ($\Sigma_{s,1\to2}$ −13%), so $M^2$ grows by 15%.
* **1B → 1C → 1D: hotter fuel.** Only the fast absorption changes, +1.9%
  and +1.5%. That is Doppler broadening of the U-238 resonances. The
  thermal constants hardly move, because the water temperature is the
  same.

### Two checks

**Against OpenMC.** The OpenMC runs of chapter 17 already tally the rates
and fluxes needed to homogenize the cell their own way, condensed with
OpenMC's own infinite-lattice spectrum (`reference.few_group_reference`).
Against this code's constants with the infinite spectrum:

| 1A | this code | OpenMC | difference |
| --- | --- | --- | --- |
| $\phi_1/\phi_2$ | 6.186 | 6.172 | +0.23% |
| $\Sigma_{a,1}$ | 0.010354 | 0.010362 | −0.07% |
| $\Sigma_{a,2}$ | 0.10046 | 0.10001 | +0.45% |
| $\nu\Sigma_{f,1}$ | 0.0073010 | 0.0073392 | −0.52% |
| $\nu\Sigma_{f,2}$ | 0.14981 | 0.14901 | +0.54% |

The thermal constants are both about 0.5% high, and their ratio, which is
what $k$ depends on, agrees to 0.1%. The cause is spatial. This code's
thermal flux in the fuel, relative to the water, is 0.42% higher than
OpenMC's (0.8983 against 0.8946), so the fuel has a little more weight in
the homogenized thermal constants. Our thermal spectrum is also very
slightly cooler (mean energy −0.13%). Both point at thermal transport in the
water. The prime suspect is the outflow transport correction for
hydrogen's strongly forward-peaked scattering (chapter 16).

**Through the nodal code.** `make nodal-check` runs the nodal diffusion code
of branch `claude-nodal`, unchanged, on the 1A constants (in the tools
environment, which has its scipy dependency):

* an all-reflective core gives $k = 1.183641$: the two-group $k_\infty$;
* with $D_g B^2$ added to the absorption, $k = 1.000000$: the critical
  buckling carries through;
* a bare square of side $\pi\sqrt{2/B^2} = 77.9$ cm, with vacuum on its
  sides, gives $k = 1.0187$. Above 1, because the diffusion flux vanishes
  only at the extrapolated boundary, about $2D$ beyond the physical one,
  so the core acts larger than its geometric buckling says.

## Further reading

* [[Hébert](#ref-hebert)] (not free): homogenization, the B1 method and
  leakage models, discontinuity factors and SPH equivalence.
* [[Duderstadt & Hamilton](#ref-dh)] (not free): multigroup diffusion
  theory and few-group constants.
* [[ENDF-102](#ref-endf102)] §1.5: the components of fission energy release
  (free).
* Code: [`fewgroup/`](../../src/fission/fewgroup/) (`homogenize.py`,
  `b1.py`, `condense.py`, `kappa.py`),
  [`checks/nodal_fewgroup.py`](../../checks/nodal_fewgroup.py);
  `make fewgroup`, `make nodal-check`.
