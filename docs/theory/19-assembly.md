# The assembly, the clad and a second transport method

Chapters 15–18 built everything on one pin cell. This chapter takes three
steps away from it:

1. a **second transport method**, collision probabilities in a cylindrical
   cell, as an independent check of MOC;
2. the **clad in the slowing-down calculation**, the error that chapter 17
   suspected;
3. the **17×17 assembly** of VERA problem 2: guide tubes, an instrument
   tube, the water gap between assemblies, pin powers, and the
   discontinuity factors that chapter 18 left at 1.

## Collision probabilities

MOC follows tracks through the exact square cell. The classic alternative
for a pin cell [[Bell & Glasstone](#ref-bell) §2.8] makes two
approximations:

* **The Wigner–Seitz cell.** The square becomes a circle of the same area,
  $R = p/\sqrt\pi$.
* **A white boundary.** Every neutron that leaves comes back isotropically,
  instead of being reflected like a mirror.

In exchange, the transport problem becomes small and exact: once per group,
compute the first-flight collision probabilities $P_{ij}$ between the rings.
With flat sources $q_i$, the flux then follows from

$$
\Sigma_j\,\phi_j\,V_j = \sum_i q_i\,V_i\,P_{ij}.
$$

**Chords and Bickley functions.** In an infinite cylinder only the height
$y$ of a chord matters. A neutron moving at polar angle $\theta$ travels an
optical path $\tau/\sin\theta$ for a projected path $\tau$. Integrating
over the polar angle turns the exponential attenuation into the
Bickley–Naylor functions

$$
\mathrm{Ki}_n(x) = \int_0^{\pi/2}\cos^{n-1}u\;e^{-x/\cos u}\,du,\qquad
\frac{d\,\mathrm{Ki}_n}{dx} = -\mathrm{Ki}_{n-1}.
$$

Integrating $\mathrm{Ki}_1$ over a source segment $a$ and a target segment
$b$ further along the chord gives the number of first collisions in $b$ per
unit source density in $a$. Let $\tau_{ab}$ be the optical distance between
the two segments:

$$
N_{ab} = \frac{\mathrm{Ki}_3(\tau_{ab}) - \mathrm{Ki}_3(\tau_{ab}+\tau_a) - \mathrm{Ki}_3(\tau_{ab}+\tau_b) + \mathrm{Ki}_3(\tau_{ab}+\tau_a+\tau_b)}{\Sigma_a}.
$$

A segment with itself gives $[\tau_a - \mathrm{Ki}_3(0) + \mathrm{Ki}_3(\tau_a)]/\Sigma_a$,
and what passes the end of the chord escapes. Summed over the chords these
give $V_iP_{ij}$.

For the weights of the sum: averaged over directions, every direction of
travel sees the same chords. The chords at $-y$ repeat those at $y$, so the
weight is $2\int_0^R dy$. Along one direction only downstream pairs count.
The ring sequence $N\ldots k\ldots N$ is a palindrome, so those pairs
include both $i\to j$ and $j\to i$.

**The white boundary.** A neutron entering through the surface collides
first in ring $j$ with probability $P_{sj} = 4V_j\Sigma_jP_{j,\mathrm{esc}}/S$
(reciprocity), and crosses the cell without colliding with probability
$P_{ss} = 1-\sum_jP_{sj}$. Every escape returns as often as needed:
$P^w_{ij} = P_{ij} + P_{i,\mathrm{esc}}P_{sj}/(1-P_{ss})$.

**Tests** (`tests/test_cp.py`):

* $\mathrm{Ki}_n(0)$ and the derivative rule hold.
* Conservation, $\sum_jP_{ij}+P_{i,\mathrm{esc}} = 1$, holds exactly.
* Reciprocity, $V_i\Sigma_iP_{ij} = V_j\Sigma_jP_{ji}$, holds to quadrature
  accuracy. Nothing in the code imposes it, so it checks the $\mathrm{Ki}_3$
  combinations.
* A black cylinder escapes through its surface only:
  $P_{\mathrm{esc}}\to S/(4V\Sigma)$.
* A cell of one material gives its $k_\infty$, with a flat flux.

**On the VERA pin.** Solved on exactly the cross sections MOC uses
(shielded, transport-corrected, with SPH; `make pin PIN_ARGS="--rings --sph --cp"`),
the cylindrical cell gives $k = 1.18590$, **+71 pcm above MOC**. That is
the expected size and sign of the Wigner–Seitz approximation for a PWR
pin: the circle with a white boundary randomizes the directions of the
neutrons leaving the cell, and slightly overestimates the share of them
that reach the next fuel pin without a collision. Two different transport
methods agree to within the known error of the simpler one.

## The clad slows down too

Chapters 16 and 17 treat everything outside the fuel as a narrow-resonance
"outside". A neutron that collides there is replaced by one from the
asymptotic flux $\psi_M$. That is right for hydrogen, which takes a
neutron far from the resonance in one collision. It is wrong for
zirconium. A neutron scattered in the clad loses at most 4.4% of its
energy, so it stays near the resonance where it was. In a U-238 resonance
the clad's own flux dips too, because the clad sees the fuel next to it.
Treating the clad as narrow-resonance overstates the source returning to
the pellet, and the capture in it.

The cure is to march the clad with the pellet rings. Each region keeps its
own composition: the fuel rings have U and O, the clad Zr, Sn, Fe, Cr and
Hf. The balance at every energy point is

$$
\Sigma_{t,i}V_i\psi_i = \sum_jV_jP_{j\to i}S_j + \Sigma_{t,i}V_iE_i\psi_M,
$$

with each nuclide's elastic source in every region where it is present.
The first-collision probabilities now depend on two cross sections that
both vary with energy, the fuel's and the clad's. `coupled_collisions`
tabulates them on a 43 × 21 grid of $(\Sigma_{\mathrm{fuel}},\Sigma_{\mathrm{clad}})$,
one MOC fixed-source sweep per pair, and interpolates bilinearly in the
logarithms. The march is the one of chapter 17, now with a
cross section and a coefficient per region. Its Rust port agrees with
Python to $10^{-12}$ on a mixed-region test. SPH then corrects rings and
clad together.

Against the OpenMC runs of chapter 17 (`--rings --sph --clad`):

| | 1A | 1B | 1C | 1D |
| --- | --- | --- | --- | --- |
| chapter 17 − OpenMC | −162 pcm | −160 pcm | −245 pcm | −348 pcm |
| with the clad − OpenMC | **−77 pcm** | **−62 pcm** | −139 pcm | −234 pcm |

The excess U-238 capture between 1 eV and 25 keV falls from 0.0016 to
0.0009 per fission neutron (1A). Chapter 17 estimated the effect at +60 pcm
by crudely halving the clad's cross section in the escape model. It is
+85 to +115 pcm.

The Doppler steps do not change. 1B→1C is −1084 pcm against OpenMC's
−1007, and 1C→1D is −921 against −826. The clad was an error in the level,
not in the trend. What remains grows with the fuel temperature and sits in
the same scattering resonances (20.9, 36.7, 66, 102 and 208 eV). The
leading suspect is the scattering kinematics. The march scatters on nuclei
at rest, while below 400 kT (19.5 eV at 565 K, 41 eV at 1200 K) OpenMC
samples the thermal motion of the target. That motion grows with
temperature.

### Moving targets: chasing the Doppler trend

Before changing the method, rule out the data. The Doppler-broadened
U-238 tables of chapter 8 give a shielded capture resonance integral
($\sigma_0 = 50$ b, integrated over a union grid of both tables) of
16.0662, 16.6584 and 17.1702 b at 600, 900 and 1200 K. NJOY's broadened
tables of the same file give 16.0679, 16.6601 and 17.1719 b. The levels
differ by 0.01% and the temperature steps (+3.686% and +3.072%) are the
same to the last digit. The cross sections are not the cause.

What is left is what happens *in* a collision. A U-238 nucleus in the
fuel is not at rest. It moves with a Maxwellian velocity at the fuel
temperature, and the elastic scattering kernel depends on the *relative*
velocity of neutron and nucleus. There are three levels of physics:

1. **Target at rest.** A neutron of energy $E$ leaves with an energy
   between $\alpha E$ and $E$, $\alpha = ((A-1)/(A+1))^2 = 0.983$ for
   U-238. It can never gain energy. This is what the march of chapter 16
   does, at every temperature.
2. **Free gas with a constant cross section** (OpenMC's default below
   400 kT; [[OpenMC](#ref-openmc-physics) §5.10.1]). The target velocity
   is drawn from the Maxwellian, weighted by the relative speed only. The
   neutron can now gain energy as well as lose it, and the spread of the
   outgoing energies grows like $\sqrt{E\,kT/A}$. But the weighting ignores
   that $\sigma_s$ itself depends on the relative speed.
3. **Free gas with the energy-dependent cross section**
   ([[OpenMC](#ref-openmc-physics) §5.10.2]): the target velocity is
   weighted by $v_r\,\sigma_s(v_r)$, where $\sigma_s$ is the *0 K*
   cross section at the relative energy. Near a resonance, collisions with
   targets whose velocity puts the pair *on* the resonance peak are far
   more likely than the constant-cross-section model assumes. Becker,
   Dagan and Lohnert proved the rejection form of this
   ([[Becker](#ref-becker) DBRC]); Romano and Walsh sample the relative
   velocity directly ([[Romano & Walsh](#ref-rvs) RVS]). This is the exact
   free-gas kernel.

At 36.7 eV the Doppler width $\Gamma_D = \sqrt{4EkT/A}$ is 0.18 eV at
600 K and 0.25 eV at 1200 K. One scatter on a U-238 nucleus at rest moves a
neutron by up to $1-\alpha = 1.7\%$ of its energy, 0.62 eV. Kernel and resonance live on the same energy scale. Where a
neutron lands after a scatter in the wing of a resonance, inside it or
beside it, depends on the kernel, and so does how many neutrons are
captured there.

OpenMC can run all three models on the same data. A switch in the
reference spec selects them (`spec --free-gas-cutoff 4` puts targets at
rest above 4 eV; `--resonance-scattering` uses RVS for U-238 between 1 eV
and 1 keV; `--tag` keeps the runs apart). The pin cells 1B (600 K) and 1D
(1200 K), with this code's rings, SPH factors and clad:

| | 1B | 1D | 1B → 1D |
| --- | --- | --- | --- |
| OpenMC, constant-σ free gas below 400 kT (default) | 1.18162 | 1.16329 | −1833 pcm |
| OpenMC, target at rest above 4 eV | 1.18158 | 1.16188 | −1970 pcm |
| OpenMC, exact kernel (RVS) for U-238 | 1.18083 | 1.15982 | **−2101 pcm** |
| this code (target at rest) | 1.18100 | 1.16095 | −2005 pcm |
| KENO (ENDF/B-VII.0) | 1.182149 | 1.162603 | −1955 pcm |

The OpenMC statistics are ±13 to ±15 pcm per run. Three things follow:

* **The trend is the kinematics.** With targets at rest in both codes,
  our Doppler step and OpenMC's agree within 35 pcm (−2005 against −1970).
  Earlier, this code was compared with a reference that uses different
  kinematics. That difference grows with temperature: at 600 K the rest
  and default runs agree (−4 pcm), while at 1200 K they differ by −141 pcm.
* **The default reference was the weak one.** The constant-cross-section
  model puts the neutron's energy change in the wrong place in a
  resonance, and that makes the Doppler step *too weak*. The exact kernel
  makes it 268 pcm stronger than OpenMC's default (−2101 against −1833).
  So the "8–12% too strong" Doppler steps of chapters 17 and 19 were
  measured against an approximate reference.
* **Against the exact kernel, this code is now too weak.** We are +17 pcm
  at 1B and +113 pcm at 1D, and our step is 96 pcm (4.6%) too weak. The
  U-238 capture between 1 eV and 25 keV changes from +0.00019 per fission
  neutron at 1B to −0.00073 at 1D. With targets at rest we miss the extra
  capture that thermal motion causes.

The physically right reference for a Doppler coefficient is the exact
kernel, and the right fix is in the method: the march needs the free-gas
kernel with the 0 K cross section.

### The exact free-gas kernel in the march

**The kernel.** Measure speeds in the reduced unit of Doppler broadening,
$x = \sqrt{AE/kT}$ (chapter 8): $x$ and $x'$ for the neutron before and
after the collision, $z$ for the relative speed of neutron and nucleus
(relative energy $E_r = z^2 kT/A$). In the centre-of-mass frame the
scattering is isotropic, and the 0 K cross section $\sigma_0(E_r)$ decides
how often a given relative speed collides. Average over the Maxwellian of
the target and over the scattering angle. With $a = A/(A+1)$ and $w$ the
reduced speed of the centre of mass, the result is
[[Ouisloumen & Sanchez](#ref-ouisloumen)]

$$
\sigma_s(E \to E') = \frac{(A+1)^2}{2A^2E\sqrt{\pi}}
\int_0^\infty z\,\sigma_0(E_r)\,J(z)\,dz,
\qquad
J(z) = \int_{w_-}^{w_+} e^{\,x^2/A - z^2/(A+1) - w^2/a}\,dw,
$$

with $w_- = \max(|x - az|, |x' - az|)$ and $w_+ = \min(x, x') + az$. The
exponent is minus the target's squared speed, so it is never positive. It
is also at most $-(x - z)^2$: only relative speeds within a few units of
$x$ contribute. Three checks follow from the formula, and the code passes
all of them (`tests/test_moving.py`):

* Integrated over $E'$ it gives the Doppler-broadened $\sigma_s(E)$ of
  chapter 8, an independent calculation of the same average. Across the
  6.67 eV resonance the two agree to the mesh error. That error is 0.13%
  at the default spacing and falls as $\Delta x^2$; the stored kernel is
  normalized, so the march conserves neutrons anyway.
* It obeys detailed balance, $E e^{-E/kT}\sigma(E\to E') = E'
  e^{-E'/kT}\sigma(E'\to E)$, to $10^{-14}$.
* As $T \to 0$ it becomes the box $1/((1-\alpha)E)$ of the target at rest.
  For a constant $\sigma_0$, its integral is the closed-form free-gas
  $\sigma_s$ of chapter 8.

**Computing it.** For $x > 10$ (above 0.04 eV for U-238 at 1200 K) the upper
limit $w_+ \approx 2x$ contributes nothing. Then $J$ depends on $x'$ only
through which of $|x - az|$ and $|x' - az|$ is larger, and that switches
once, at $z^* = (x + x')/2a$. On each side $J$ depends on only one of the
two speeds:

$$
J(z) = D(x, z)\ \text{or}\ e^{(x^2 - x'^2)/A}\,D(x', z),\qquad
D(x, z) = e^{x^2/A - z^2/(A+1) - (x-az)^2/a}\,\tfrac{\sqrt{\pi a}}{2}\,
\mathrm{erfcx}\!\left(\tfrac{|x-az|}{\sqrt a}\right).
$$

The factor $e^{(x^2-x'^2)/A} = e^{(E-E')/kT}$ is detailed balance again.
The $z$ integral becomes a sum on a uniform grid with weights $W_j = \int
\mathrm{hat}_j(z)\,z\,\sigma_0\,dz$. The weights are exact for the 0 K
table, which is linear between its points, and every resonance enters only
through them. Each kernel entry is then a prefix sum of $W D$ for one node
plus a suffix sum for the other. Building the kernel costs one $D$ array
per node, not one integral per pair. The kernel lives on a mesh uniform in
$x$, with spacing $\Delta x = 1/8$ of a Doppler width. For U-238 between
1 eV and 1 keV at 600 K it has 4.7 million entries and takes 3 s to build
in Rust (`moving.py`, `moving.rs`). Halving $\Delta x$ moves $k$ of 1D by
1 pcm.

**In the march.** The source of a moving-target nuclide at $E$ is $E\sum_m
P(E_m \to E)\,\Phi_m$, where $\Phi_m$ is the hat moment of
$\Sigma_s\psi$ around mesh node $m$. The mesh energies are grid points of
the march. A moving target can also give the neutron energy, so the source
at $E$ needs the flux *below* $E$, which the march has not reached yet.
The march is therefore repeated (Gauss–Seidel iteration). The first pass
treats the targets at rest. Later passes take $\Phi$ from the current pass
where the march has already been and from the previous pass elsewhere,
until $\psi$ changes by less than $10^{-10}$. Neutrons under the bottom of
the march (1 eV) scatter up into it too. Their $\Sigma_s\psi$ is continued
flat below the last point. Without them the flux would dip by 10% at the
bottom.

A pure scatterer is a sharp test. The march keeps $\psi = 1$ to $+5\times
10^{-4}$, and that residual is physics: hot nuclei give energy back, about
$2kT/E$ per collision. There is one feature that is not a test failure. At
the top of the kernel's range, targets above move and targets below are at
rest. Neutrons can scatter up across this boundary, but none come back down
through the at-rest box. The flux shows a seam of several percent just
below it. OpenMC's resonance scattering has the same seam at its upper
energy, 1 keV here.

**Results.** U-238 scatters on moving targets between 1 eV and 1 keV
(`--moving`), as in OpenMC's RVS run. To compare like with like, OpenMC
was run once more with RVS for U-238 and every other nuclide at rest above
4 eV. This code treats them that way, and the RVS run of the table above
did not (constant-σ free gas below 400 kT):

| | 1B | 1D | 1B → 1D |
| --- | --- | --- | --- |
| OpenMC, all targets at rest above 4 eV | 1.18158 | 1.16188 | −1970 pcm |
| this code, targets at rest | 1.18100 | 1.16095 | −2005 pcm |
| OpenMC, RVS for U-238, others at rest | 1.18127 ± 0.00020 | 1.15955 ± 0.00019 | **−2172 ± 28 pcm** |
| this code, moving targets for U-238 | 1.18003 | 1.15795 | **−2208 pcm** |

The kinematics now match the reference, and so does the Doppler step: to
36 pcm (1.7%, 1.3σ), the same 35 pcm as with targets at rest. Moving
targets steepen our step by 203 pcm, and OpenMC's by 202. The level moved
too. Moving targets lower our $k$ by 97 pcm at 600 K and 300 pcm at
1200 K, OpenMC's by 31 ± 24 and 233 ± 24. This code now sits 124 and
160 pcm below the matched reference, against 58 and 94 with targets at
rest. A temperature-independent 66 ± 18 pcm of extra capture is not
explained yet. The excess U-238 capture between 1 eV and 25 keV,
+0.0009 to +0.0015 per fission neutron, is the place to look. The same
excess was the remaining error of chapter 17 with targets at rest.

### Chasing the level: the discretization was not converged

Two OpenMC eigenvalues each carry about ±17 pcm, so their difference
carries ±24. A quieter measure of the same change comes from the tallies:
$k_{\text{tally}} = \nu\Sigma_f / \Sigma_a$ summed over the cell. It sits a
constant 190 pcm below OpenMC's own eigenvalue in every run, but its
difference between two runs is much less noisy. On it, OpenMC's moving
targets lower $k$ by 70 pcm at 600 K and 267 pcm at 1200 K, and ours by 97
and 300. So the gap is about 30 pcm, not 66. Divided into energy bands,
the extra U-238 capture sits at 25–50 eV (the 36.7 eV resonance) and
10–25 eV (20.9 eV). These are the resonances with the largest neutron
widths, and they were also the worst bands with targets at rest
(+2.6% at 25–50 eV for 1D).

**Pellet rings.** In these resonances the flux dips sharply within a
fraction of a millimetre of the pellet surface. Three equal-area rings put
the outer 0.075 cm into one flat-flux region. With more rings (1D, and with
the tracking fixed as described next), we get:

| fuel rings | 3 | 6 | 10 |
| --- | --- | --- | --- |
| $k$ at rest − OpenMC at rest | +92 | +142 | +173 pcm |
| U-238 capture 25–50 eV vs OpenMC, at rest | +1.4% | +0.7% | +0.3% |
| moving − at rest | −294 | −283 | −277 pcm |

Moving targets broaden exactly these resonances. A coarse rim region
over-captures in them, and the broadening makes that error larger. The
"extra" moving-target capture was mostly this spatial error.

**Tracking.** Running the final 172-group MOC on finer tracks moved $k$
more than the rings did. To separate the two, the shielding, SPH factors and
libraries were computed once with the default tracking (64 azimuthal
directions over 2π, 0.02 cm spacing, 3 Gauss–Legendre polar angles), and
only the final eigenvalue solve was repeated. For 1D with three rings:

| azimuthal / spacing (cm) / polar | Δk vs default |
| --- | --- |
| 64 / 0.01 / 3 | +35 pcm |
| 64 / 0.02 / 6 | +46 pcm |
| 128 / 0.02 / 3 | +86 pcm |
| 128 / 0.01 / 6 | +134 pcm |
| 256 / 0.005 / 6 | +168 pcm |
| 256 / 0.005 / 12 | +186 pcm |

The same test with infinitely dilute libraries (no self-shielding at all)
moves $k$ by +560 pcm from the default to 512 / 0.005 / 12, still rising
by 10–25 pcm per doubling of any one parameter. This is MOC itself, not
the shielding. The test also shows why the shift is larger without
shielding: unshielded resonance groups are optically thicker.

The reason is the balance of an optically thick region. If
$\Sigma_t \ell \gg 1$ along every segment, the sweep's region flux becomes

$$
\bar\phi_i \approx \frac{q_i}{\Sigma_t} + \frac{1}{\Sigma_t A_i}
\sum_{\text{entering segments}} W\,\left(\psi_\text{in} - \frac{q_i}{\Sigma_t}\right),
\qquad W = 2\pi\, w_a\, \omega_p\, \delta_a \sin\theta_p .
$$

The second term is the inflow through the region's surface, and the
tracks evaluate it as a quadrature. For an isotropic inflow its polar part
is $\int_0^1 \sin\theta\,d\mu = \pi/4$. The integrand $\sqrt{1-\mu^2}$
has an infinite slope at $\mu = 1$, so Gauss–Legendre in $\mu$ converges
slowly: it overestimates the integral by 0.46% with 3 points, 0.069% with
6 and 0.0095% with 12. A larger estimated inflow means more neutrons enter
the thick resonance region of the pellet rim, so there is too much capture
and too low a $k$. The sign matches. At intermediate optical thickness the
polar integrand becomes $\sin\theta\, e^{-\tau/\sin\theta}$ (a Bickley
function), which is harder still. Lattice codes use polar sets built for
these integrals for this reason (Tabuchi–Yamamoto, chapter 15 and
[[OpenMOC docs](#ref-openmocdocs)]). Spacing matters in the same term:
the number of tracks crossing a small ring is right only on average.

The default tracking was verified on C5G7 (chapter 15). Its seven groups
have no optically thick resonance groups. For
the 172-group pin cell it is about 190 pcm low. The agreement with OpenMC
quoted in chapters 17–19 was partly this error cancelling the opposite
one in the 6.67 eV resonance (below).

**Shielding inputs.** Collision probabilities, escape tables and SPH
factors are MOC calculations too. Making every one of them use 128 / 0.01 /
6 instead of only the final solve adds another 30 pcm, mostly from 50 to
300 eV. This is smaller, and it is left out of the table below.

**With both converged.** The final solve on 256 / 0.005 / 12 with 10
rings gives, compared with the matched OpenMC runs:

| | U-238 capture added by moving targets (10⁻⁵ per fission neutron) | moving − rest, $k$ |
| --- | --- | --- |
| 1B, this code | +63 | −86 pcm |
| 1B, OpenMC | +50 ± 8 | −70 pcm (tally) |
| 1D, this code | +209 | −277 pcm |
| 1D, OpenMC | +201 ± 8 | −267 pcm (tally) |

The moving-target kernel agrees with OpenMC's RVS to 1–1.6σ. At 600 K a
small 36.7 eV residual of +9 ± 3 remains. So nothing was wrong in the
kinematics. The 66 pcm came from OpenMC noise (about half) and from the
spatial and angular discretization, which moving targets amplify.

**What the converged level shows.** With the discretization fixed, this
code is +169 pcm (1B) and +173 pcm (1D) above OpenMC at rest. All of that
is one band: U-238 capture at 4–10 eV, the 6.67 eV resonance, is 2.6% too
low (−0.0012 per fission neutron) at both temperatures and for every ring
count. Every other band is now within 0.5%. Before, the polar error added
capture at the pellet rim and hid most of this. The 6.67 eV resonance
is capture-dominated ($\Gamma_n$ = 1.5 meV, $\Gamma_\gamma$ = 22.7 meV in
ENDF/B-VIII.1), the opposite of the 36.7 eV one ($\Gamma_n$ = 33.6 meV,
$\Gamma_\gamma$ = 22.3 meV).

### The 6.67 eV deficit: hydrogen moves too

A capture deficit can mean two things: the cross section the flux sees is
wrong (self-shielding), or the flux itself is wrong. OpenMC also tallies the
172-group flux in every region, so the two can be told apart. Group by group,
the effective cross section (capture over flux) agrees within 0.1–0.3% outside
the resonance groups. The **flux** is low, by the same amount in every
region, the moderator included (1D, 3 rings, converged tracking):

| group (eV) | 22.6–25.0 | 9.2–9.9 | 7.5–8.3 | 4.1–5.0 | 4.0–4.1 | 3.4–4.0 | 372–454 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| moderator flux vs OpenMC | −0.6% | −1.5% | −2.0% | −2.6% | −1.2% | −0.8% | −0.1% |

A deficit that grows as $1/E$ (deficit $\times E \approx$ 12–16 %·eV, a few
$kT/E$) and jumps at 4 eV points at the moderator, not the fuel. Below 4 eV
the library used the S(α,β) law of H in water. Above 4 eV hydrogen was a
target **at rest**. OpenMC puts a target at rest only if it is heavier than
a neutron. Its source has the test
`E >= free_gas_threshold * kT && nuc.awr_ > 1.0` in `sample_target_velocity`
(`src/physics.cpp`, v0.16.0), and H-1 has $A$ = 0.99917. So in OpenMC
hydrogen is bound (S(α,β)) up to the end of the water data at 10 eV, and a
free gas at all energies above.

The physics is simple. On a proton at rest a neutron keeps on average $E/2$.
On a proton moving thermally it keeps $E/2 + 0.75\,kT$, which this code's
free-gas kernel gives at 5, 10 and 20 eV (565 K). Each collision takes a
little less energy, so neutrons spend longer at each energy, and the
epithermal flux rises by a few $kT/E$: 2.6% at 4.5 eV. Bound hydrogen in
water moves even faster than the water's temperature suggests. The water
law's short-collision-time temperature at 565 K is 1285 K. That is why the
law matters up to 10 eV.

The library now follows the same rule (chapter 11): a scattering law up to
its own end (9.906 eV for water in XMAS-172), then a free gas for a nuclide
lighter than the neutron up to 1 keV, then at rest. Each step was tested
on its own. A free gas from 4 eV to 1 keV takes the 1D level from +93 to
−26 pcm (3 rings). Moving the law's end to 10 eV adds −96 pcm and brings
the flux between 4 and 25 eV to within ±0.5% everywhere. The whole
reaction balance changes. The deficit at 6.67 eV and the excess of thermal
U-235 fission (+55 → +3 × 10⁻⁵ per fission neutron) were cancelling each
other.

With 10 rings and the final solve on 256 / 0.005 / 12:

| | $k$ − OpenMC at rest | U-238 capture 1–4 eV | 4–10 eV | 25–50 eV | 50–100 eV |
| --- | --- | --- | --- | --- | --- |
| 1B, H at rest above 4 eV | +169 pcm | −1.1% | −2.7% | +0.2% | +0.1% |
| 1B, H moving | −44 pcm | −0.0% | −0.6% | +0.5% | +0.4% |
| 1D, H at rest above 4 eV | +173 pcm | −1.1% | −2.6% | +0.3% | +0.4% |
| 1D, H moving | −42 pcm | −0.0% | −0.5% | +0.7% | +0.7% |

Every band of U-238 capture is now within 1%. The rest of the level is the
shielding inputs' tracking (+30 pcm, see above) and the remaining capture
excess in the large-$\Gamma_n$ resonances. The moving-target effect of
U-238 does not change: 1B −86 pcm, 1D −278 pcm. The march still treats
hydrogen at rest above 1 eV. It only shapes the flux inside groups, so this
matters little; it is a candidate for the last −0.5% at 4–10 eV.

**The default settings now show their error.** With the default tracking
(64 / 0.02 / 3) and three rings, the pin compared with OpenMC's default runs
moves from −77 / −62 / −139 / −234 pcm (1A–1D) to −271 / −274 / −351 /
−447 pcm. The two errors no longer cancel: about 200 pcm of missing
epithermal flux is gone, and about 270 pcm of discretization (tracking and
rings) is left in plain view. The numbers quoted earlier in chapters 12–19
were computed with hydrogen at rest above 4 eV. The assembly results above
have not been rerun.

### Converged defaults

Two changes make the default pin cell converged.

**Polar angles.** Tabuchi–Yamamoto replaces Gauss–Legendre (chapter 15). Its
3 angles reproduce $\mathrm{Ki}_3$ to 3.4×10⁻⁵, better than 12 Gauss–Legendre
angles. On the unshielded 1D pin at 512 / 0.005, 3 Tabuchi–Yamamoto angles
are within 5 pcm of 12 Gauss–Legendre angles.

**Azimuthal angles and rings.** With the polar error gone, what remains is
the number of directions in the plane. In a resonance group the pellet rim
is nearly black, so its capture equals the current of neutrons entering the
pellet. The tracks count that current by direction, and the inflow into a
square cell's pin depends on direction (the corners and the flat sides
differ). For the shielded 1D pin, with the self-shielding done once and the
final solve repeated:

| azimuthal directions | 64 | 128 | 256 | 512 |
| --- | --- | --- | --- | --- |
| $k$ − 512 / 0.005 cm, spacing 0.01 cm | −99 | −58 | −14 | |

Spacing below 0.02 cm changes $k$ by about 10 pcm. Running everything
(shielding, SPH and the final solve) on each candidate, 1D at rest:

| tracking, rings | $k$ − OpenMC at rest | time |
| --- | --- | --- |
| 64 / 0.02, 3 rings (old default) | −220 pcm | 30 s |
| 128 / 0.02, 10 rings | −74 pcm | 2 min |
| **256 / 0.02, 10 rings** (new default) | **−23 pcm** | 3.5 min |
| 256 / 0.01, 10 rings | −12 pcm | 7 min |

The times are with the libraries already cached; a new state adds about
2 minutes to build them. Most of the time is in the fixed-source solves of
the self-shielding: the ring and clad collision probabilities on a 43 × 21
grid of cross sections, and the SPH iterations. Two changes halved it.
Each solve stops when the region fluxes stop changing; checking every
boundary angular flux after every sweep cost more than the sweep. SPH
iterations also start from the previous iteration's boundary fluxes.

Against OpenMC with the new defaults (`--rings --sph --clad`):

| | 1A | 1B | 1C | 1D |
| --- | --- | --- | --- | --- |
| this code | 1.18619 | 1.18069 | 1.17069 | 1.16165 |
| vs OpenMC default | −62 ± 14 | −92 ± 15 | −85 ± 15 | −164 ± 15 pcm |
| vs OpenMC, targets at rest above 4 eV | | −89 ± 13 | | −23 ± 15 pcm |

Against the matched runs (targets at rest), every reaction rate agrees within
10⁻⁴ per fission neutron. Summed, those rates put $k$ within 10 pcm at both
temperatures. The spread of the eigenvalue differences is mostly OpenMC's:
two eigenvalues of ±13–15 pcm each, and its eigenvalue sits 173–207 pcm
above its own tally ratio, varying from run to run. The default OpenMC
runs also sample free-gas motion of every nuclide below 400 kT, which
this code (targets at rest) leaves out. That is why 1D is lower against
them. C5G7 with Tabuchi–Yamamoto: +12 pcm, pin powers within 0.26%.

The assembly keeps 64 / 0.03 and 3 rings (`--fuel-rings`): 256 directions
on 17×17 would take hours. It gets the Tabuchi–Yamamoto polar set like
everything else.

Moving targets are therefore an option, not the default. They make the
temperature dependence right. For a Doppler coefficient (the M2 branches)
they belong in.

## The 17×17 assembly

VERA problem 2 is one Westinghouse 17×17 lattice in 2D with reflective
sides [[VERA](#ref-vera) Tables 2, P2-1 to P2-5]. It has 264 fuel rods, 24
guide tubes and one instrument tube (Zircaloy-4, filled with water in
cases 2A–2D), on a 1.26 cm pitch. The model file is `models/vera-p2.toml`.
The guide-tube map is Figure 3 of the specification.

**The gap.** Assemblies sit 21.50 cm apart, 0.08 cm more than 17 × 1.26.
Each lattice is therefore surrounded by 0.04 cm of water, half of the gap
between neighbours. The model format gains a lattice key `gap`. The
transport mesh adds no regions for it. A point in the gap belongs to the
outermost region of the nearest edge cell, in the same sector, by angle
from that cell's centre. The areas of those regions grow by their share of
the gap strip, computed exactly as wedges of an off-centre rectangle. The
ray tracer (Python and Rust) shifts its cell lines by the gap. CMFD uses
the true widths of the edge cells. Its corrected couplings make the
converged answer independent of that choice, but the right widths keep the
acceleration effective. The area test (tracks against exact areas) and the
Rust parity test cover a gapped lattice.

**Water-only cells and CMFD.** The transport correction makes hydrogen's
within-group scattering negative in the fast groups: its P1 moment exceeds
its P0 within-group term. With the first, flat-flux guess, a cell of pure
water (a guide tube) can then get a negative source and a negative coarse
flux, and CMFD divides by it. The pin cell never showed this, because its
only coarse cell contains fuel. The converged flux is positive, so the
solver now runs plain power iterations until every coarse cell flux is
positive, then accelerates. In practice this takes a few iterations.

**Shielding.** Every fuel rod shares one shielded library per pellet ring.
The escape probabilities are computed over the whole assembly, so they
average over rods next to guide tubes, at the edge and in the middle.
Rod-by-rod Dancoff factors would be the next refinement.

### Results

`make assembly STATES=2A PIN_ARGS="--rings --sph --clad"` runs one state:
15,384 flat-source regions, 172 groups, about 11 minutes. The OpenMC
references use the same chain as chapter 17, with the lattice and the
gap: 10⁸ neutron histories for 2A (±9 pcm, pins ±0.14%) and 4×10⁷ for
2B–2D.

| state | this code | OpenMC (VIII.1) | KENO (VII.0) | − OpenMC | − KENO |
| --- | --- | --- | --- | --- | --- |
| 2A | 1.18191 | 1.18202 ± 0.00009 | 1.182175 | −11 pcm | −26 pcm |
| 2B | 1.18298 | 1.18344 ± 0.00014 | 1.183360 | −46 pcm | −38 pcm |
| 2C | 1.17313 | 1.17392 ± 0.00015 | 1.173751 | −80 pcm | −63 pcm |
| 2D | 1.16471 | 1.16588 ± 0.00015 | 1.165591 | −117 pcm | −88 pcm |

* **The data.** The two Monte Carlo codes agree within 30 pcm, so for the
  assembly too, ENDF/B-VII.0 and VIII.1 are interchangeable.
* **The method.** The assembly sits closer to the references than the pin
  cell does (−77 pcm for 1A). The Doppler steps are again stronger than
  OpenMC's default ones, 2B→2C by 3.5% and 2C→2D by 4.7%. This is the
  kinematics trend of the pin cell, smaller: see "Moving targets" above.
* **The mesh matters more here.** With 32 directions and 0.05 cm between
  tracks, 2A comes out at 1.18075, −127 pcm. With 64 directions and 0.03
  cm it is −11 pcm. The 0.04 cm gap and the 0.057 cm clad are thinner than
  the coarse track spacing, and they enter every self-shielding table. The
  assembly default is therefore the finer mesh.

**Pin powers.** Fission rate per rod, normalized to a fuel average of 1:

| | 2A | 2B | 2C | 2D |
| --- | --- | --- | --- | --- |
| vs KENO: RMS (max) | 0.07% (0.22%) | 0.06% (0.15%) | 0.07% (0.17%) | 0.05% (0.11%) |
| vs OpenMC: RMS | 0.13% | 0.19% | 0.20% | 0.20% |
| OpenMC's own uncertainty per rod | 0.14% | 0.21% | 0.21% | 0.21% |

Against KENO, whose pin powers carry about 0.02% uncertainty, the shape is
right to better than a tenth of a percent. The highest rod in 2A, next to
a guide tube, is 1.0519 against 1.0516. Against OpenMC the
differences are its statistical noise. That holds although every rod
shares one shielded library: the radial power shape is set mostly by the
water around the guide tubes, and the transport calculation resolves it.

### Discontinuity factors

A homogeneous node has a flat flux. The real assembly's flux on its surface
differs from its average, and the nodal code needs the ratio
$f_g = \phi_{s,g}/\bar\phi_g$ to reproduce the currents between unlike
neighbours [[Hébert](#ref-hebert)]. With reflective sides the
homogeneous flux is exactly the average, so $f$ needs only the surface
flux.

MOC has the surface flux almost for free. Every track starts and ends on
the boundary, and the solver keeps the angular flux entering there. A track
of 2D direction $\mathbf d$ crosses a side of normal $\mathbf n$ once per
$|\mathbf d\cdot\mathbf n|/\delta$ of length. Each entering crossing
therefore adds $W\psi/(\sin\theta\,|\mathbf d\cdot\mathbf n|)$ to the
side's integrated scalar flux, with the weights $W$ of the face currents.
On a reflective side the outgoing flux mirrors the incoming one, so the
incoming half of the sphere carries half of $\phi$: a factor 2. A test
checks the normalization. In a uniform medium the surface flux equals the
volume flux, to the tracking error.

For 2A, condensed to two groups (B1 spectrum: $B^2 = 3.23\times10^{-3}$
cm⁻²; $D$ = 1.372 / 0.410 cm, $\Sigma_a$ = 0.00941 / 0.0943 cm⁻¹,
$\nu\Sigma_f$ = 0.00666 / 0.1360 cm⁻¹, $\Sigma_{s,1\to2}$ = 0.01730 cm⁻¹):

| 2A | fast | thermal |
| --- | --- | --- |
| discontinuity factor, this code | 0.9962 | 0.9995 |
| flux in the 0.04 cm gap / average, OpenMC | 0.9937 | 1.0013 |

The two columns are not the same quantity. MOC gives the flux exactly on
the surface. OpenMC gives the average over the thin water strip next to
it, which comes close because the normal derivative vanishes at a
reflective side. They agree to 0.2–0.3%, and both say the same thing:
for a lattice of only UO₂ rods and empty guide tubes, the edge looks like
the average, and $f \approx 1$.

The factors are sensitive to the tracking. With 32 directions and 0.05 cm,
the fast factor is 0.9995. Discontinuity factors will matter where
neighbours differ: next to a reflector, a rodded assembly, or MOX fuel.
These are the cases for the M track.

## Further reading

* [[Bell & Glasstone](#ref-bell)] §2.8 (escape and collision probabilities
  in cylinders, Bickley functions) (free).
* [[Hébert](#ref-hebert)] (not free): collision probabilities in
  cylindrical cells, white boundaries, discontinuity factors.
* [[VERA](#ref-vera)] problem 2: geometry, states, pin powers (free).
* [[Ouisloumen & Sanchez](#ref-ouisloumen)] (not free): the free-gas
  kernel with the 0 K cross section; [[Becker](#ref-becker)] and
  [[Romano & Walsh](#ref-rvs)] (not free) for its Monte Carlo forms.
* Code: [`transport/cp.py`](../../src/fission/transport/cp.py),
  [`pincell/shielding.py`](../../src/fission/pincell/shielding.py)
  (`coupled_collisions`, `coupled_flux`),
  [`pincell/assembly.py`](../../src/fission/pincell/assembly.py),
  [`pincell/moving.py`](../../src/fission/pincell/moving.py) (free-gas kernel),
  [`transport/moc.py`](../../src/fission/transport/moc.py) (`surface_flux`),
  [`models/vera-p2.toml`](../../models/vera-p2.toml);
  `make assembly STATES=2A PIN_ARGS="--rings --sph --clad"`.
