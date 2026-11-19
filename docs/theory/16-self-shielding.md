# Self-shielding in a pin cell

Chapter 12 shielded the cross sections of a *homogenized* pin cell, and its
$k_\infty$ fell short of the reference by 4500 pcm: in the real cell the fuel
is a lump, and resonance neutrons are absorbed at its surface before they
reach the inside. This chapter keeps the regions apart. MOC (chapter 15)
solves the cell with fuel, gap, clad and water at their own temperatures,
using the ENDF/B-VIII.1 libraries of chapters 10 and 11. The resonance cross
sections of fuel and clad are shielded in four steps of increasing accuracy.
All four rest on the collision-probability treatment of heterogeneous systems
in [[Bell & Glasstone](#ref-bell) §2.8 and §8.4].

## The slowing-down balance of a fuel pin

Take two regions, the fuel F and everything else M, and neutrons at energy
$E$ inside a resonance. Assume a **flat source** in each region. Let $P(E)$ be
the probability that a neutron starting in F makes its next collision in M.
For neutrons at one energy every collision removes them (absorption or
scattering to another energy), so $P$ is a one-group, purely absorbing
quantity. It depends on $E$ only through the fuel's total cross section,
$P = P(\Sigma_{t,F}(E))$.

Collisions in F come from neutrons scattered in F that stay there, and from
neutrons scattered in M that reach F. By reciprocity the second term is
$P\,\Sigma_{t,F}\,\phi_M$. If the moderator has no resonances, its flux
per unit lethargy is smooth (the asymptotic value, 1). This gives
[[Bell & Glasstone](#ref-bell) eqs. 8.83–8.85]

$$
\Sigma_{t,F}(E)\,\phi_F(E) = \big(1 - P\big)\,S_F(E) + P\,\Sigma_{t,F}(E),
$$

where $S_F$ is the slowing-down source from scattering inside the fuel. The
methods below differ in how they approximate $P$ and $S_F$.

**Narrow resonances (NR).** If a resonance is narrow compared with the energy
a neutron loses in a collision with the fuel nuclides, the neutrons scattered
into it come from far above, where the flux is asymptotic: $S_F = \Sigma_p$,
the fuel's potential scattering ($\sigma_p = 4\pi a^2$ with the ENDF
scattering radius $a$). Then

$$
\phi_F(E) = P + (1 - P)\,\frac{\Sigma_p}{\Sigma_{t,F}(E)} .
$$

## Step 1: equivalence theory

The **rational approximation** writes $P = \Sigma_e/(\Sigma_{t,F} + \Sigma_e)$
with an energy-independent **escape cross section** $\Sigma_e$. For a fuel
rod in a lattice,

$$
\Sigma_e = \frac{1 - C}{\bar\ell},
\qquad \bar\ell = \frac{4V}{S}
$$

[[Bell & Glasstone](#ref-bell) eqs. 2.112, 2.122]. Here $\bar\ell$ is the
mean chord of the rod ($2r$ for a cylinder), and $C$ is the **Dancoff
correction**: the fraction of neutrons heading for a pin that are
intercepted by the neighbouring pins. With this form the NR flux becomes
$(\Sigma_p + \Sigma_e)/(\Sigma_t + \Sigma_e)$, the flux of a homogeneous
mixture with background cross section

$$
\sigma_{0,i} = \frac{\sum_{j\ne i} N_j\sigma_{p,j} + \Sigma_e}{N_i}
$$

per atom of the resonant nuclide $i$. This is the **equivalence principle**
[[Bell & Glasstone](#ref-bell) eqs. 8.87–8.88]: a heterogeneous cell
behaves like the homogeneous mixture of chapter 12 with this $\sigma_0$.
Each nuclide's group constants are averaged with
$w(E)/(\sigma_{t,i}(E) + \sigma_{0,i})$, which is the same calculation that
produces Bondarenko tables, done here directly at the cell's $\sigma_0$.

**The Dancoff correction from MOC.** $C$ is computed for the actual lattice
with a one-group fixed-source MOC calculation:

* the pellet is made black ($\Sigma = 1000$ cm⁻¹);
* every other region gets its potential-scattering cross section and a
  source $\Sigma/4\pi$, so that without the pellet the flux would be 1
  everywhere.

A pin alone in such a medium absorbs $S/4$ per unit height, the inward
partial current of an isotropic unit flux; in the lattice it absorbs
$(1-C)\,S/4$. The unshadowed current is computed with the same tracks, with
the black pin surrounded by an optically thick medium, so that the tracking
error cancels in the ratio. Checks: a transparent moderator gives $C = 1$; a
pin alone in a 12.6 cm water cell gives $|C| < 10^{-4}$.

For VERA 1A: $C = 0.342$, $\bar\ell = 0.819$ cm, $\Sigma_e = 0.804$ cm⁻¹,
and U-238 sees $\sigma_0 = 44$ b. For the clad, $C = 0.54$.

## Step 2: the exact escape probability

[[Bell & Glasstone](#ref-bell) §2.8b] note that the rational form
underestimates escape. With the same one-group MOC, $P(\Sigma)$ can be
computed exactly for any fuel cross section: put a flat source in the pellet
and count the collisions there. For the VERA pin the rational form is low by
up to 16% (0.446 against 0.516 at $\Sigma = 1$ cm⁻¹, left figure below).
Too little escape means too much shielding.

The method "cp" therefore uses the NR flux with the tabulated exact $P$. It
needs one more ingredient: $\Sigma_{t,F}(E)$, the pointwise total of the
*whole* fuel mixture. One flux shape then serves all fuel nuclides, and the
resonances of U-235 and U-238 shield each other where they overlap.

## Step 3: solving the slowing down

NR is exact only for narrow resonances. The low U-238 resonances (6.67, 20.9
and 36.7 eV) are not narrow for scattering on U-238 itself: a neutron loses
at most $1 - \alpha = 1.7\%$ of its energy in such a collision, so a neutron
scattered inside the resonance stays near it. The opposite limit, wide
resonances (WR), ignores the absorber's resonance scattering altogether. For
VERA 1A the two limits bracket the reference by thousands of pcm:

| U-238 resonance treatment | $k_\infty$ | vs KENO |
| --- | --- | --- |
| NR everywhere | 1.19701 | +997 pcm |
| WR below 200 eV | 1.16497 | −2207 pcm |
| WR everywhere | 1.14236 | −4468 pcm |

So neither is acceptable. The method "sd" solves the balance numerically
[[Bell & Glasstone](#ref-bell) §8.4c]. In lethargy $u$, with $\psi$ the flux
per unit lethargy,

$$
\Sigma_t(u)\,\psi(u) = (1 - P)\sum_i \frac{1}{1-\alpha_i}
\int_{u-\Delta_i}^{u}\Sigma_{s,i}(u')\,\psi(u')\,e^{-(u-u')}\,\mathrm du'
+ P\,\Sigma_t(u)\,\psi_M(u) .
$$

Each fuel nuclide $i$ contributes its own elastic kernel, isotropic in the
centre-of-mass system, with $\Delta_i = \ln(1/\alpha_i)$. The equation is
solved on the union of the pointwise grids between 1 eV and 24.8 keV:
272,000 points for the VERA fuel, about 2 s in Python. Above that range the
NR flux supplies the incoming neutrons. Some details matter:

* **Marching.** Each nuclide keeps the running integral
  $C_i(u) = \int \Sigma_{s,i}\psi\,e^{u'}\mathrm du'$ (trapezoid rule). The
  window is $C_i(u) - C_i(u - \Delta_i)$, and the unknown $\psi(u)$ enters
  linearly through the last half-step.
* **Consistency.** Pure scattering conserves neutrons, so nothing damps a
  small bias in the source. Interpolating $C$ linearly at the lower window
  edge left an error of $6\times10^{-4}$ in a constant-cross-section test.
  Integrating the partial step with the same trapezoid rule as $C$ itself
  reduces it to $7\times10^{-6}$.
* **Very heavy nuclides.** When a window is shorter than a grid step, the
  whole window lies in the last step and must be treated implicitly too.
* **Moderator level.** $\psi_M$ starts at 1 and is depleted by the
  absorption accumulated so far, relative to the cell's slowing-down power
  $\sum_r V_r\,\xi_r\Sigma_{p,r}$. For VERA this changes $k$ by 5 pcm.

The tests check the limits on synthetic data:

* pure scattering keeps $\psi = 1$;
* constant absorption gives the decay rate $\Sigma_a/(\xi\Sigma_t)$;
* an absorption resonance on a light absorber follows NR;
* a scattering resonance on a nucleus of mass $10^6$ shows no flux dip
  (WR);
* with mass $10^4$ the dip comes back. That nucleus loses little per
  collision, but its ~50 collisions per moderator collision drift a neutron
  through the resonance.

![Left: escape probability from the VERA fuel pin, computed by MOC and in the rational approximation. Right: the fuel flux around the 36.7 eV resonance of U-238, from the NR formula and from the slowing-down solution. Neutrons scattered in the resonance reappear just below it, and the absorption of the resonances above lowers the flux everywhere below.](../figures/self_shielding.svg)

## Results

VERA problem 1A, heterogeneous cell: 56 flat-source regions, 64 directions,
0.02 cm spacing, 3 polar angles, XMAS-172 groups, MOC with CMFD. Each run
takes 15–60 s.

| | $k_\infty$ | vs KENO | U-238 absorption per fission neutron |
| --- | --- | --- | --- |
| homogenized, NR (chapter 12) | 1.14194 | −4510 pcm | 0.303 |
| heterogeneous, unshielded | 0.89872 | −28832 pcm | 0.539 |
| heterogeneous, equivalence | 1.20189 | +1486 pcm | 0.246 |
| heterogeneous, exact $P$ with NR ("cp") | 1.19701 | +997 pcm | 0.250 |
| heterogeneous, slowing down solved ("sd") | 1.17991 | −713 pcm | 0.263 |

Even unshielded, the heterogeneous cell is much more reactive than the
homogenized one with dilute data (0.899 against 0.743). MOC resolves the
flux depression inside the pellet at the level of whole groups, which is
part of the self-shielding.

The **temperature trends** are less sensitive to the nuclear data than $k$
itself, and they agree much better:

| change | "sd" | KENO | chapter 12 (homogeneous NR) |
| --- | --- | --- | --- |
| 1A → 1B: hotter, less dense water | −530 pcm | −489 pcm | −501 pcm |
| 1B → 1C: fuel 600 → 900 K | −1150 pcm | −1043 pcm | −1492 pcm |
| 1C → 1D: fuel 900 → 1200 K | −978 pcm | −912 pcm | −1271 pcm |

The Doppler coefficient is now within about 10% of the reference; the
homogeneous NR calculation was 40% too strong. The absolute values are
1.17991, 1.17461, 1.16311 and 1.15333 for 1A–1D, between 713 and 927 pcm
below the reference.

**What the offset is not.** Varying the calculation for the "cp" case:

| change | effect on $k$ |
| --- | --- |
| no transport correction | −23 pcm |
| 128 directions, 0.01 cm, 5 fuel rings, 16 sectors | +70 pcm |
| 6 polar angles | +54 pcm |
| depleting the moderator flux in "sd" | −5 pcm |

**What it may be.** Several approximations remain, each plausibly worth
tens to a few hundred pcm:

* one flux shape for the whole pellet, although the rim is less shielded
  than the centre;
* infinitely dilute cross sections in the unresolved range;
* $\chi$ weighted by the library spectrum;
* the outflow transport correction for hydrogen;
* the data themselves: the reference used ENDF/B-VII.0.

Separating these needs an independent calculation with the same data:
continuous-energy Monte Carlo on ENDF/B-VIII.1. Chapter 17 does it. The
data turn out to account for only 24 pcm; almost all of the offset is
method, and most of it is the first two items: the multigroup treatment of
the shielded groups (corrected by SPH factors) and the pellet's radial
shielding.

## Further reading

* [[Bell & Glasstone](#ref-bell)] §2.8 (escape probabilities, the rational
  approximation, the Dancoff correction, Table 2.9) and §8.4 (resonance
  absorption in heterogeneous systems: collision probabilities, equivalence
  relations, numerical solution) (free).
* [[Hébert](#ref-hebert)] (not free): self-shielding in lattice codes,
  including subgroup methods and SPH factors.
* Code: [`shielding.py`](../../src/fission/pincell/shielding.py) (Dancoff
  correction, escape probabilities, background cross sections, the
  slowing-down march), [`heterogeneous.py`](../../src/fission/pincell/heterogeneous.py)
  (transport-corrected macroscopic data, the MOC pin cell, the four
  shielding options), [`moc.py`](../../src/fission/transport/moc.py)
  (`tracking`, `fixed_source`); `make pin STATES=1A SHIELD=sd`.
