# Resonance physics: the common ingredients

All three resonance formalisms in chapters 4–6 are built from the same pieces:
the compound nucleus and its widths, angular momentum and spin statistics,
the neutron wave number, and the "hard-sphere" functions that describe a
neutron outside the nucleus. This chapter introduces them once.

## The compound nucleus

Bohr's picture (1936): a slow neutron that enters a nucleus shares its energy
among all nucleons almost immediately, forming a **compound nucleus**
$^{A+1}X^*$ that "forgets" how it was made. The excitation energy is

$$
E^* = S_n + \frac{A}{A+1}\,E ,
$$

where $S_n$ is the neutron separation energy (4.806 MeV for U-239, the Q value
of U-238 capture) and $\frac{A}{A+1}E$ is the neutron's energy in the
centre-of-mass frame. Like an atom, the compound nucleus has discrete quantum
states. When $E^*$ matches one of them, the neutron is captured into that state
far more readily, and the cross section has a **resonance**.

The state does not live forever. It decays through **channels**: re-emitting
the neutron (elastic scattering), emitting gamma rays (capture), or splitting
(fission). Each channel $c$ has a **partial width** $\Gamma_c$, related to the
decay rate $\lambda_c$ by $\Gamma_c = \hbar \lambda_c$. The total width
$\Gamma = \sum_c \Gamma_c$ sets the lifetime

$$
\tau = \frac{\hbar}{\Gamma}, \qquad \hbar = 6.582\times10^{-16}\ \text{eV s}.
$$

For the U-238 resonance at 6.67 eV, $\Gamma = 24$ meV, so $\tau \approx
2.7\times10^{-14}$ s: about $10^8$ times longer than the time a neutron needs to
cross the nucleus. That large ratio is why the compound nucleus "forgets" its
formation. The *branching ratio* $\Gamma_c/\Gamma$ is the probability of
decaying through channel $c$; this independence of formation and decay is the
**Bohr hypothesis**, and every single-level formula below has the structure
(formation probability) $\times$ (branching ratio).

### What the parameters look like

Each resolved resonance in ENDF is a row of numbers: $E_r$, $J$, $\Gamma_n$,
$\Gamma_\gamma$, and for fissile nuclides $\Gamma_f$. The first U-238
resonances:

| $E_r$ (eV) | $l$ | $J$ | $\Gamma_n$ (meV) | $\Gamma_\gamma$ (meV) |
| --- | --- | --- | --- | --- |
| 4.41 | 1 | 3/2 | 0.0000556 | 22.5 |
| **6.674** | 0 | 1/2 | **1.492** | **22.7** |
| 10.24 | 1 | 3/2 | 0.00085 | 22.5 |
| **20.87** | 0 | 1/2 | **10.08** | **22.8** |
| **36.68** | 0 | 1/2 | **33.59** | **22.3** |

Two patterns hold across all 3345 resolved U-238 resonances:

* **$\Gamma_\gamma$ is nearly constant** (≈ 22.5 meV). Gamma decay proceeds
  through thousands of paths to lower states; the sum of many random
  contributions has little relative spread.
* **$\Gamma_n$ varies over many orders of magnitude.** Neutron emission from a
  given state goes through essentially *one* channel, whose amplitude is a
  random-looking projection of a complicated wave function. Its square follows
  the **Porter–Thomas** distribution (a $\chi^2$ distribution with one degree of
  freedom) [[Porter & Thomas](#ref-porterthomas); [Fröhner](#ref-frohner)
  §4.1.1], so very small widths are common. The p-wave resonances are smaller
  still for a separate reason, the centrifugal barrier (section 3.4).

Resonances with **negative energy** ($E_r < 0$) also appear: 22 of them for
U-238. They are compound states just *below* the neutron separation energy.
They cannot be reached by a real neutron, but their tails extend into the
positive-energy region and shape the thermal cross section. Evaluators
often adjust them so that the thermal cross sections come out right
[[Fröhner](#ref-frohner) §3.4.3].

## Angular momentum and the spin factor $g_J$

A neutron with orbital angular momentum $l$ (in units of $\hbar$; $l = 0, 1, 2$
are called s-, p-, d-waves) meets a target of spin $I$. The neutron spin
$i = \tfrac12$ couples with $I$ into the **channel spin** $s = |I - \tfrac12|$
or $I + \tfrac12$, and $s$ couples with $l$ into the total angular momentum
$J$, running from $|l - s|$ to $l + s$. A resonance is a state of definite
$J$ (and parity $(-1)^l$), so only the matching combinations can form it.

For unpolarized neutrons and targets, all $(2i+1)(2I+1)$ spin orientations are
equally likely; a state of spin $J$ has $2J+1$ substates. The fraction able to
form it is the **statistical spin factor**

$$
g_J = \frac{2J+1}{(2i+1)(2I+1)} = \frac{2J+1}{2(2I+1)} .
$$

Summed over every $(s, J)$ combination allowed for one $l$, the weights add up
to $2l+1$, the number of orientations of the orbital angular momentum:

$$
\sum_{s}\sum_{J=|l-s|}^{l+s} g_J = 2l+1 .
$$

Examples:

* U-238 ($I = 0$): s-wave resonances have $J = \tfrac12$, $g = 1$; p-wave
  resonances have $J = \tfrac12$ ($g = 1$) or $\tfrac32$ ($g = 2$), total
  $3 = 2l+1$.
* U-235 ($I = \tfrac72$): s-wave resonances have $J = 3$ ($g = \tfrac{7}{16}$)
  or $J = 4$ ($g = \tfrac{9}{16}$).

When a spin combination has no resonances, it still contributes hard-sphere
scattering (section 3.4). Keeping track of these "missing" weights matters in
chapters 5 and 6.

## The wave number

In the centre-of-mass frame the neutron–target system has reduced mass
$\mu = m_n M/(m_n + M)$ and kinetic energy $E_{cm} = E\, M/(m_n + M)$. The
wave number is $k = \sqrt{2\mu E_{cm}}/\hbar$. Writing
$\text{AWRI} = M/m_n$,

$$
k = \frac{\sqrt{2 m_n}}{\hbar}\,\frac{\text{AWRI}}{\text{AWRI} + 1}\,\sqrt{|E|}
\qquad\text{[ENDF-102, eq. D.10]},
$$

and with $m_n c^2 = 939.565$ MeV and $\hbar c = 197.327$ MeV fm (CODATA 2018
[[CODATA](#ref-codata)]),

$$
k = 2.196\,81\times10^{-3}\ \frac{\text{AWRI}}{\text{AWRI}+1}\sqrt{E\,[\text{eV}]}
\quad\text{in }(10^{-12}\,\text{cm})^{-1}.
$$

(The ENDF manual quotes $2.196771\times10^{-3}$, from older constants; the
difference, 17 ppm, is far below the accuracy of any resonance parameter.)

The **reduced wavelength** $1/k$ is the quantum "size" of the neutron. At
6.67 eV on U-238, $1/k = 177\times10^{-12}$ cm $\approx 1800$ fm, about two
hundred times the nuclear radius ($a = 8.4$ fm). The natural cross-section
scale is

$$
\frac{\pi}{k^2} = 98\,400\ \text{b at 6.67 eV},
$$

which is how a nucleus of geometric cross section $\pi a^2 \approx 2$ b can
capture with a cross section of 22,800 b.

Because the wave is so much larger than the nucleus, only small $l$ matter.
Semi-classically, a neutron with angular momentum $l\hbar$ passes the nucleus at
impact parameter $b \approx l/k$; it can only touch it if $b \lesssim a$, so
$l \lesssim ka$. For U-238 at 1 keV, $ka \approx 0.06$: s-waves dominate and
p-waves are strongly suppressed. Section 3.4 makes that quantitative.

## Outside the nucleus: the hard-sphere functions

Outside the nuclear radius $a$, the neutron moves freely. Its radial wave
function $u_l(r)$ satisfies the free radial Schrödinger equation

$$
\frac{\mathrm{d}^2 u_l}{\mathrm{d}r^2} + \left[k^2 - \frac{l(l+1)}{r^2}\right] u_l = 0 ,
$$

whose solutions are the **Riccati–Bessel functions** of $\rho = kr$
[[DLMF](#ref-dlmf) §10.47; [Lane & Thomas](#ref-lanethomas)]:

$$
F_l(\rho) = \rho\, j_l(\rho) \ (\text{regular}),\qquad
G_l(\rho) = -\rho\, y_l(\rho) \ (\text{irregular}),
$$

for example $F_0 = \sin\rho$, $G_0 = \cos\rho$,
$F_1 = \sin\rho/\rho - \cos\rho$, $G_1 = \cos\rho/\rho + \sin\rho$. The term
$l(l+1)/r^2$ is the **centrifugal barrier**. Evaluated at the channel radius,
$\rho = ka$, three combinations describe everything resonance theory needs from
the outside region:

$$
\begin{aligned}
P_l &= \frac{\rho}{F_l^2 + G_l^2} &&\text{penetrability (penetration factor)}\\[2pt]
S_l &= \frac{\rho\,(F_l F_l' + G_l G_l')}{F_l^2 + G_l^2} &&\text{shift factor}\\[2pt]
\varphi_l &= \arctan\frac{F_l}{G_l} &&\text{hard-sphere phase shift}
\end{aligned}
$$

For $l = 0 \dots 3$ they have closed forms [[ENDF-102](#ref-endf102) Table D.1,
eqs. D.11–D.13]:

| $l$ | $P_l(\rho)$ | $S_l(\rho)$ | $\varphi_l(\hat\rho)$ |
| --- | --- | --- | --- |
| 0 | $\rho$ | $0$ | $\hat\rho$ |
| 1 | $\dfrac{\rho^3}{1+\rho^2}$ | $-\dfrac{1}{1+\rho^2}$ | $\hat\rho - \arctan\hat\rho$ |
| 2 | $\dfrac{\rho^5}{9+3\rho^2+\rho^4}$ | $-\dfrac{18+3\rho^2}{9+3\rho^2+\rho^4}$ | $\hat\rho - \arctan\dfrac{3\hat\rho}{3-\hat\rho^2}$ |
| 3 | $\dfrac{\rho^7}{225+45\rho^2+6\rho^4+\rho^6}$ | $-\dfrac{675+90\rho^2+6\rho^4}{225+45\rho^2+6\rho^4+\rho^6}$ | $\hat\rho - \arctan\dfrac{\hat\rho(15-\hat\rho^2)}{15-6\hat\rho^2}$ |

These closed forms are checked in the test suite against the Bessel-function
definitions above, an independent route to the same numbers.

![Left: penetrabilities. At small $\rho$, $P_l \approx \rho^{2l+1}/[(2l-1)!!]^2$: each extra unit of $l$ costs a factor $\rho^2$. Right: $\sin^2\varphi_l$, which controls potential scattering.](../figures/hard_sphere.svg)

**Penetrability.** $P_l$ measures how easily the neutron reaches the nuclear
surface through the centrifugal barrier; it multiplies the neutron width. For
U-238 at 1 keV, $\rho = ka = 0.058$ and $P_1/P_0 = \rho^2/(1+\rho^2) = 0.0034$:
p-wave neutrons are about 300 times less likely to get in. This is why most U-238
p-wave resonances are tiny, and why only s-waves matter in the thermal range.

**Phase shift and potential scattering.** Even without resonances a nucleus
scatters neutrons, like an impenetrable sphere. That **potential scattering**
cross section is

$$
\sigma_{\text{pot}} = \frac{4\pi}{k^2}\sum_l (2l+1)\sin^2\varphi_l
\;\xrightarrow{\;k \to 0\;}\; 4\pi\, AP^2 ,
$$

four times the geometric area. For U-238, $AP = 0.948\times10^{-12}$ cm gives
$4\pi\,AP^2 = 11.3$ b, the plateau visible between resonances.

**Shift factor.** $S_l$ makes the apparent resonance energy drift slightly
with $E$ for $l > 0$ (section 4.1). It vanishes for s-waves.

## Channel radius versus scattering radius

Two radii appear in the formulas [[ENDF-102](#ref-endf102) eq. D.14]:

* the **channel radius** $a$, used in $P_l$ and $S_l$ ($\rho = ka$);
* the **scattering radius** AP, used in the phase shift ($\hat\rho = k\,AP$).
  It is fitted by the evaluator to reproduce the measured potential
  scattering.

The flag NAPS selects $a$: for NAPS = 0, $a = 0.123\,\text{AWRI}^{1/3} + 0.08$
(in $10^{-12}$ cm, a nuclear-size formula); for NAPS = 1, $a = AP$.
Reich–Moore evaluations may give a separate AP for each $l$ (called APL).

| Nuclide | Formalism | NAPS | $a$ ($10^{-12}$ cm) | AP |
| --- | --- | --- | --- | --- |
| U-238 | Reich–Moore | 0 | 0.840 | 0.948 (APL) |
| U-235 | Reich–Moore | 1 | $= AP$ | 0.960 (APL, $l = 0$) |
| Zr-90 | MLBW | 0 | 0.630 | 0.716 |

## Energy dependence of the widths

ENDF gives $\Gamma_n$ *at the resonance energy*. Since the neutron must
penetrate the barrier, the width at another energy scales with the
penetrability [[ENDF-102](#ref-endf102) eq. D.7]:

$$
\Gamma_n(E) = \Gamma_n(|E_r|)\,\frac{P_l(E)}{P_l(|E_r|)} .
$$

R-matrix theory writes $\Gamma_n = 2P_l\gamma_n^2$ with an energy-independent
**reduced width** $\gamma_n^2$; the formula above is the same statement. For a
negative-energy level the penetrability is evaluated at $|E_r|$ (an ENDF
convention). The gamma and fission widths are taken as constant in energy.

For s-waves, $P_0 = ka \propto \sqrt{E}$, so $\Gamma_n \propto \sqrt E$. This
is one of the two ingredients of the $1/v$ law.

## Formalisms in ENDF

Resolved resonance parameters can be given in several formalisms (the LRF flag
in MF2) [[ENDF-102](#ref-endf102) §2.2]:

| LRF | Formalism | Level interference | Used here by |
| --- | --- | --- | --- |
| 1 | Single-Level Breit–Wigner (SLBW) | none | (reference, chapter 4) |
| 2 | Multi-Level Breit–Wigner (MLBW) | elastic only | U-234, Zr-90…96 |
| 3 | Reich–Moore | full (R-matrix) | U-235, U-238 |
| 7 | R-Matrix Limited | full, general channels | (not used by our nuclides) |

At higher energies each evaluation switches to an **unresolved range**
(LRU = 2) where only averages are given. H-1, O-16 and the boron isotopes have
no resonance parameters at all (LRP = 0); for them MF3 is complete.

## Further reading

* [[Fröhner](#ref-frohner)] §3 (resonance theory, hard-sphere functions) and
  §4.1 (width and level statistics) (free, NEA/JEFF Report 18). The best single
  free text on these topics.
* [[ENDF-102](#ref-endf102)] Appendix D.1 and Table D.1 (free).
* [[Lane & Thomas](#ref-lanethomas)] §III–IV: the original derivations (not
  free).
* [[Blatt & Weisskopf](#ref-blattweisskopf)] ch. VIII–X: compound nucleus and
  resonance reactions (not free; inexpensive reprint).
* Code: [`channel.py`](../../src/nodal/endf/channel.py) (wave number,
  $P_l$, $S_l$, $\varphi_l$, $g_J$), [`resonances.py`](../../src/nodal/endf/resonances.py)
  (reading MF2), and the Rust port [`channel.rs`](../../rust/src/channel.rs).
