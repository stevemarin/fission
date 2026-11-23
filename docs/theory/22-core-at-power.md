# The core at power

Chapter 21 solved a core whose constants were given. A real core at power
sets its own constants. The power heats the fuel, the fuel temperature
changes the U-238 resonance absorption (Doppler), the heat goes into the
water, which gets hotter and lighter as it rises through the core, and less
dense water moderates less. Every node's constants depend on its state
(chapter 20), and its state depends on the power. This chapter closes that
loop for an illustrative core built entirely from this code's own
constants.

## Reflector constants

The core is surrounded by water. Its nodes need two-group constants, and
they cannot come from water alone: a reflector has no source, and its
spectrum is the one the fuel puts into it. Fast neutrons leak in from the
fuel, slow down in the water, and the thermal flux peaks a few centimetres
out, where they are thermal but the fuel can no longer absorb them. That
peak is what makes the reflector return neutrons to the core.

The constants therefore come from a calculation with fuel next to water:
one row of 17 VERA fuel cells (one assembly width, 21.42 cm) followed by
17 cells of water, solved with the heterogeneous MOC of chapters 16–19
(`lattices.reflector` in `models/vera-p1.toml`). The model format now
allows such rectangular lattices. With reflective sides the slab repeats
as fuel | water | water | fuel, so the water between two fuel slabs is
42.84 cm thick. Only the water cells are homogenized, with the slab's own
flux; there is no B1 spectrum, since there is nothing to make critical.
The diffusion coefficient is the transport one, $D_g = 1/(3\Sigma_{tr,g})$
with $\Sigma_{tr} = \Sigma_t - \sum_{g'}\Sigma_{s1,g\to g'}$, condensed so
that each fine group's current $D_g\phi_g$ is kept. CMFD does not converge
on this long slab of fuel and water, so it runs with plain power
iteration.

With fuel at 600 K and water at 600 K and 0.661 g/cm³ (state 1B):

| boron | slab $k$ | $D_1$, $D_2$ (cm) | $\Sigma_{a,1}$, $\Sigma_{a,2}$ (1/cm) | $\Sigma_{1\to2}$ (1/cm) |
| --- | --- | --- | --- | --- |
| 0 ppm | 1.14160 | 2.177, 0.304 | 0.00037, 0.00898 | 0.0382 |
| 1300 ppm | 0.99131 | 2.175, 0.308 | 0.00101, 0.03044 | 0.0380 |
| 2000 ppm | 0.93650 | 2.175, 0.310 | 0.00136, 0.04155 | 0.0377 |

The thermal absorption is the water's own (hydrogen capture, 0.009 cm⁻¹)
plus the boron, which more than triples it at 1300 ppm. These form an M3
table with boron as the only axis (`make reflector`).

## Water: IAPWS-IF97

The coolant enters at 565 K and 15.5 MPa and leaves some 35 K hotter. Its
density falls by about 10% on the way, and the constants follow. The
industrial standard for water properties is IAPWS-IF97
[[IAPWS-IF97](#ref-if97)]. Its region 1, compressed liquid up to 623.15 K,
covers a PWR core. It gives the dimensionless Gibbs free energy as a sum of
34 terms,

$$
\gamma(\pi, \tau) = \frac{g}{RT} = \sum_{i=1}^{34} n_i\,(7.1 - \pi)^{I_i}\,(\tau - 1.222)^{J_i},
\qquad \pi = \frac{p}{16.53\ \text{MPa}},\quad \tau = \frac{1386\ \text{K}}{T},
$$

and every property follows from its derivatives: the specific volume
$v = \pi\gamma_\pi RT/p$, the enthalpy $h = \tau\gamma_\tau RT$ and the heat
capacity $c_p = -\tau^2\gamma_{\tau\tau}R$, with $R$ = 461.526 J/(kg K). The
implementation (`feedback/water.py`) reproduces the release's verification
values (its Table 5) to all nine printed digits. At 2250 psia (15.513 MPa) it
gives 742.7 kg/m³ at 565 K and 661.2 kg/m³ at 600 K, the 0.743 and 0.661
g/cm³ of VERA states 1A and 1B. A temperature for a given enthalpy comes
from Newton's method, whose slope is $c_p$.

## The fuel rod

A rod producing $q'$ watts per metre, with uniform heat generation in the
pellet, has (cylindrical conduction)

$$
\int_{T_s}^{T_c} k(T)\,\mathrm{d}T = \frac{q'}{4\pi}\ \text{(pellet)},\qquad
T_s - T_{ci} = \frac{q'}{2\pi r_p h_\text{gap}}\ \text{(gap)},\qquad
T_{ci} - T_{co} = \frac{q'\ln(r_{co}/r_{ci})}{2\pi k_\text{clad}}\ \text{(clad)}.
$$

The pellet integral handles a conductivity that falls with temperature:
UO2 conducts worse when hot, so the centre runs hotter than a constant $k$
would give. The properties are the simplified ones of the OECD/NEA PWR
MOX/UO2 core transient benchmark specification [[MOX benchmark](#ref-moxbench)
§7]: $k_{UO_2} = 1.05 + 2150/(T - 73.15)$ W/(m K),
$k_{Zr} = 7.51 + 2.09\times10^{-2}T - 1.45\times10^{-5}T^2 + 7.67\times10^{-9}T^3$
W/(m K), $h_\text{gap}$ = 10⁴ W/(m² K). The temperature that sets the
Doppler broadening is the specification's weighted one,
$T_\text{eff} = 0.3\,T_c + 0.7\,T_s$. The heat transfer from the clad to the
water, the last step, is the subject of the next section; the
specification leaves it to code-specific correlations.

For the average Watts Bar rod, 17.8 kW/m (17.67 MW × 97.4% over 264 rods
of 3.6576 m) in 580 K water, this gives a centre at 1040 K, a pellet
surface at 671 K and $T_\text{eff}$ = 782 K.

## From the clad to the water

The heat leaves the clad through a thin layer of slow water at its surface,
where it flows by conduction. Its resistance is lumped into a film
coefficient $h$: the clad surface is $q'/(2\pi r_{co} h)$ hotter than the
bulk water. For turbulent flow along a heated wall the classic correlation
is Dittus–Boelter, in McAdams' form [[MOOSE](#ref-moose)]:

$$
\mathrm{Nu} = \frac{h D_h}{\lambda} = 0.023\,\mathrm{Re}^{0.8}\,\mathrm{Pr}^{0.4},\qquad
\mathrm{Re} = \frac{G D_h}{\mu},\qquad \mathrm{Pr} = \frac{c_p\,\mu}{\lambda}.
$$

The Reynolds number compares inertia with viscosity: the higher it is, the
thinner the slow layer and the better the transfer. The Prandtl number
compares how fast momentum and heat diffuse. In a square rod bundle the
water between the rods behaves somewhat better than in a round tube, by
Presser's factor $\psi = 0.9217 + 0.1478\,P/D - 0.1130\,e^{-7(P/D-1)}$,
which is 1.106 for the VERA lattice ($P/D$ = 1.26/0.95).

Two more water properties are needed, each with its own IAPWS release: the
viscosity [[IAPWS viscosity](#ref-iapwsvisc)] and the thermal conductivity
[[IAPWS conductivity](#ref-iapwscond)], in their form for industrial use
(densities and heat capacities from IF97). Both are a dilute-gas term
times a density term, $\bar\mu = \bar\mu_0(\bar T)\,\bar\mu_1(\bar T,\bar\rho)$.
The conductivity adds a critical enhancement $\lambda_2$, which grows near
the critical point (647 K, 22 MPa). It is still 2.6% of the total at 620 K
and 20 MPa, so it is kept; it needs the isochoric heat capacity and
$(\partial\rho/\partial p)_T$ from IF97 as well. The implementation
reproduces the releases' verification tables to every printed digit,
including all the intermediate quantities of the conductivity (its Table 7).

The flow sees the 21.50 cm square of an assembly minus its 264 rods, 24
guide tubes and instrument tube: 246.65 cm² (VERA Tables 1–2), with a wetted
perimeter that gives a hydraulic diameter $D_h = 4A/P$ = 1.118 cm. The
rated flow of 85.97 kg/s gives a mass flux of 3485 kg/(m² s). From 565 to
600 K this means Re = 4.2–5.0 × 10⁵, Pr = 0.84–0.98 and
$h$ = 39–42 kW/(m² K), so the film adds about 15 K at the average rod.
The Reynolds number is beyond the range of the correlation
(10⁴–1.2 × 10⁵): Dittus–Boelter is used for PWR channels all the same, but
here it is extrapolated.

## The coupled iteration

`feedback.solve` alternates between neutronics and thermal-hydraulics
(Picard iteration):

1. evaluate every node's constants from its material's table at its state
   (fuel temperature, water temperature and density, boron);
2. solve the core with NEM, starting from the previous solution;
3. power per node from $\kappa\Sigma_f\phi$, scaled so that every assembly
   produces its rated power on average;
4. per channel (one column of nodes, with an assembly's flow scaled by the
   column's share of an assembly's area), the enthalpy rises with the power
   below each node centre, and IF97 gives the water's temperature and
   density;
5. per fuel node, the rod temperatures from its linear heat rate;
6. repeat until $k$ and every temperature stop changing.

Reflector nodes take the water of their channel: the radial reflector the
inlet water, the top reflector the water leaving the fuel. There is no
cross flow between channels and no boiling (IF97 region 1 ends at
623.15 K, below which the code stops).

## The core

`models/wbn1-core.toml` has the shape and operating conditions of Watts
Bar Unit 1 from the VERA specification [[VERA](#ref-vera)]: 193
assemblies on a 21.50 cm pitch (Figure 9, quarter symmetry), a 365.76 cm
fuel stack, 17.67 MW and 85.97 kg/s per assembly (0.6823 Mlb/hr), inlet
water at 565 K and 2250 psia, and 1300 ppm of boron (Table P6-1). It is
not the real core. Every assembly is the VERA problem 1 fuel (3.1% UO2,
without guide tubes or burnable absorbers), the reflector is water only
(21.5 cm thick, as in VERA problem 5, without the steel baffle), and there
are 21.5 cm of water below and above the fuel.

## Fuel constants from MOC

The fuel nodes use the heterogeneous pin of chapter 19, not the homogenized
cell of chapter 20, whose Doppler feedback is too strong. The core runs at
fixed boron, so a smaller matrix suffices (`[branches.core]`): fuel at 600,
900 and 1200 K, water at 565 and 600 K, densities 0.661 and 0.743 g/cm³, all
at 1300 ppm. That is 12 points at 128 directions and 10 rings, 85 minutes.
The MOC constants have clearly weaker feedback at 1B than the homogenized
cell:

| coefficient | MOC pin | homogenized cell |
| --- | --- | --- |
| fuel temperature | −2.57 pcm/K | −3.96 pcm/K |
| moderator temperature | −0.76 pcm/K | −1.73 pcm/K |
| moderator density | +28.1 pcm per 0.01 g/cm³ | +22.5 pcm per 0.01 g/cm³ |

For this table the fit keeps 8 terms. Interpolation in $\sqrt{T_f}$ misses
the left-out 900 K node by 6.9 pcm, against 66 pcm in $T_f$.

## Results

The iteration converges in 12 steps at half power (7 s, NEM starting each
time from the previous solution, with its sweeps in Rust).

| | zero power | 50% power | 50%, no film |
| --- | --- | --- | --- |
| $k$ | 1.16722 | 1.15793 | 1.15840 |
| largest assembly power | 2.176 | 1.959 | 1.965 |
| axial offset (top − bottom) | 0.0% | −7.5% | −7.7% |
| Doppler temperature, average / max | 565 / 565 K | 680 / 912 K | 672 / 890 K |
| peak fuel centreline | 565 K | 1331 K | 1305 K |
| water outlet | 565 K | 569–600 K | 569–600 K |

The film adds 15–20 K to every temperature in the rod. The fuel runs hotter,
and the Doppler feedback takes 47 pcm more.

At zero power everything is at the inlet temperature, and the axial shape
is symmetric, as it must be with identical reflectors above and below. At
power three feedbacks show:

* **The power defect.** Going from zero to half power costs
  $\Delta\rho = 1/1.15793 - 1/1.16722$ = 687 pcm. The fuel heats up (Doppler)
  and the water thins out. This is the reactivity the boron or the control
  rods must give back to raise the power.
* **Axial offset.** The water enters cool and dense at the bottom and
  leaves hotter and lighter at the top, so the bottom moderates better and
  the power shifts down (−7.5%).
* **Radial flattening.** The hottest assemblies have the hottest fuel and
  the lightest water, so they lose more reactivity than their neighbours:
  the peak assembly falls from 2.18 to 1.96.

**Why only half power.** With a single fuel type, the power peaks in the
centre: 2.18 times the average at zero power. At full power the hottest
channel would heat its water by about twice the average rise, past 623 K,
the end of IF97 region 1, and beyond saturation at 15.5 MPa. The code stops
there (`make feedback` at `--power 1` reports water at 624.3 K). A real core
avoids this with its loading pattern: three enrichments, the most enriched
fuel at the edge, and burnable absorbers in the centre (VERA Figure 9). The
next section does the enrichments.

## The loading pattern

Watts Bar Unit 1 started its first cycle with three fuel types
[[VERA](#ref-vera) §2.1, Figure 9]: 2.11%, 2.619% and 3.10% U-235. The two
lower enrichments alternate like a chess board in the centre, and the
3.10% fuel sits on the edge. A fresh core has no burnt fuel to place, so
the enrichments do that job. The outer assemblies leak neutrons into the
reflector, so they get the most reactive fuel. In the centre, each
assembly's neighbours are less reactive than itself, which flattens the
power. `models/wbn1-zoned.toml` is that core:

```text
a b a b a b a c R
b a b a b a c c R        a  2.11%    b  2.619%    c  3.10%
a b a b a b a c R        R  baffle and water
b a b a b a c c R
a b a b b b c R R        quarter core, centre at top left
b a b a b c c R .
a c a c c c R R .
c c c c R R R . .
R R R R R . . . .
```

**Three lattice calculations.** Every fuel type is the 17×17 VERA problem
2 assembly with guide tubes (chapter 19) and its own enrichment. The model
format lets a file *extend* another: `models/wbn1-fuel-211.toml` takes
everything from `vera-p2.toml` and changes only the uranium. The U-234 and
U-236 content follows the enrichment $w$ by the formulas of VERA Table 17,
$0.007731\,w^{1.0837}$ and $0.0046\,w$ (weight percent; the test
`test_wbn1_enrichment_zones` checks the files against them). Each type gets
a small state matrix around 2B at 1300 ppm: fuel at 600 and 900 K, water at
565 and 600 K, densities 0.661 and 0.743 g/cm³. That is eight assembly
calculations per type, at 64 directions, 0.03 cm spacing and 3 rings (the
assembly settings of chapter 19), about 13 minutes each, 5.6 hours for all
three.

| fuel | $k_\infty$ (600 K, 565 K, 0.661) | fuel temperature | moderator temperature | moderator density |
| --- | --- | --- | --- | --- |
| 2.11% | 1.06368 | −2.70 pcm/K | −0.96 pcm/K | −88.2 pcm per 0.01 g/cm³ |
| 2.619% | 1.13223 | −2.50 pcm/K | −0.68 pcm/K | −47.9 pcm per 0.01 g/cm³ |
| 3.10% | 1.18170 | −2.36 pcm/K | −0.51 pcm/K | −22.1 pcm per 0.01 g/cm³ |

Two trends are worth understanding.

* **Doppler grows as the enrichment falls.** Resonance capture in U-238
  takes a larger share of the neutron economy when there is less U-235 to
  compete with it.
* **The density coefficient is negative.** Denser water should moderate
  better and raise $k$, and in the VERA pin cell it does (+28 pcm per
  0.01 g/cm³, previous table). But the boron dissolved at 1300 ppm by
  weight comes with the water: denser water carries more boron. The
  assembly also has more water than the pin cell (25 water-filled tubes and
  the gap between assemblies). With more moderator already present, extra
  moderation gains less than the extra boron costs. The lattice is
  **overmoderated**. Lower enrichment makes it worse, because the boron
  then competes with less U-235.

What a real core feels is the two together. When the water heats up, its
density falls by 2.33 kg/m³ per kelvin between 565 and 600 K (IF97). That
turns the coefficients into a total moderator temperature coefficient of
$-0.96 + 88.2 \times 0.233$ = **+19.6 pcm/K** for the 2.11% fuel, +10.5 for
2.619% and +4.6 for 3.10%. A positive moderator temperature coefficient is
what a PWR at power must avoid: hotter water would raise the power, which
heats the water further. Real cores avoid it with burnable absorbers.
The Pyrex rods of VERA Figure 9 hold down the fresh fuel's reactivity, so
less soluble boron is needed, and they take water out of the guide tubes.
They are the next step. Here the Doppler effect is still strong enough to
keep the total power coefficient negative, as the results below show.

**The baffle.** Between the outer assemblies and the water sits a
stainless-steel baffle, 2.85 cm thick and 0.19 cm from the fuel
[[VERA](#ref-vera) Table 16]. Its constants come from the slab of
"Reflector constants", rebuilt as `models/wbn1-baffle.toml`: 17 cells of
fuel, 2 cells of steel and 15 of water. The slab's cells are 1.26 cm wide,
so the steel is 2.52 cm thick with its density raised by 2.85/2.52 to keep
the same steel per unit area, and the 0.19 cm water gap is left out. The
steel is SS304 from VERA Table P3-3, all 19 nuclides listed there. Four
of them the code cannot read yet and leaves out: Mn-55 and Cr-53
(Reich–Moore with negative channel spins), Fe-57 (R-Matrix Limited with an
inelastic channel) and Fe-58 (an unresolved format). Silicon and Fe-54
came with the Pyrex work of the next section. A first version without
them, and without C and P, had $D_1$ = 1.400 cm, $\Sigma_{a,2}$ = 0.03663/cm and
$k$ 8 pcm lower in the core below. The water and steel cells are
homogenized together into one reflector material at 1300 ppm:

| | water only | baffle + water |
| --- | --- | --- |
| $D_1$, $D_2$ (cm) | 2.175, 0.308 | 1.338, 0.307 |
| $\Sigma_{a,1}$, $\Sigma_{a,2}$ (1/cm) | 0.00101, 0.03044 | 0.00226, 0.03691 |
| $\Sigma_{1\to2}$ (1/cm) | 0.0380 | 0.0270 |

The steel scatters fast neutrons without slowing them much (iron is 56
times heavier than a neutron). That shortens the fast diffusion length, and
it absorbs in both groups. It also displaces water, so less moderation
happens in the reflector. The water above and below the fuel keeps the
water-only constants. The baffle reflector also carries discontinuity
factors from the same slab, 1.269 (fast) and 0.258 (thermal). Their
derivation is in "Control rods" below, where the comparison with KENO's
power map called for them. All the results that follow include them.

**Results.** All the constants are in place, so `make feedback
FEEDBACK_MODEL=wbn1-zoned` runs at any power:

| | zero power | 50% power | 100% power |
| --- | --- | --- | --- |
| $k$ | 1.08451 | 1.08066 | 1.07621 |
| power defect | | 329 pcm | 711 pcm |
| largest assembly power | 1.482 | 1.454 | 1.422 |
| axial offset (top − bottom) | 0.0% | +1.2% | 0.0% |
| Doppler temperature, average / max | 565 / 565 K | 678 / 813 K | 810 / 1099 K |
| peak fuel centreline | 565 K | 1096 K | 1786 K |
| water outlet | 565 K | 573–592 K | 581–613 K |
| feedback iterations | 3 | 7 | 9 |

The assembly powers at full power (quarter core, centre at top left):

```text
1.282  1.422  1.249  1.346  1.138  1.168  0.884  0.637
1.422  1.266  1.388  1.199  1.266  1.036  1.061  0.624
1.249  1.388  1.221  1.323  1.116  1.127  0.821  0.568
1.346  1.199  1.323  1.160  1.234  0.973  0.889  0.418
1.138  1.266  1.116  1.234  1.184  1.017  0.691
1.168  1.036  1.127  0.973  1.017  0.876  0.438
0.884  1.061  0.821  0.889  0.691  0.438
0.637  0.624  0.568  0.418
```

What the loading pattern does:

* **The peak falls from 2.18 to 1.48** at zero power. The centre no longer
  dominates, and the chess board shows: every 2.619% assembly has more
  power than its 2.11% neighbours (1.482 against 1.330 and 1.307 next to the
  centre). The hottest channel now heats its water to 613 K, so the core
  runs at full power below the 623.15 K limit of IF97 region 1. (The code
  does not test for local boiling.)
* **The edge still loses.** Even with 3.10% fuel, the outer assemblies make
  0.42–0.64 of the average. The baffle reflects fewer thermal neutrons than
  water would.
* **The axial shape barely moves with power.** In the uniform core of the
  previous section the power slid to the bottom (−7.5% at half power),
  because cool, dense water there moderated better. Here denser water is
  *less* reactive (overmoderated), which pushes the power up. The Doppler
  effect pushes it down, because the fuel at the top is hotter: it sits in
  hotter water. The two nearly cancel.
* **The power defect is 711 pcm at full power** and comes from the Doppler
  effect. The water's contribution is positive, so the defect is smaller
  than the fuel temperature alone would give.

$k = 1.08$ at 1300 ppm is far from critical. The measured critical boron
of this core at zero power, with its Pyrex rods in place and all control
rods out, is 1291 ppm [[VERA](#ref-vera) Table 18]. The burnable
absorbers take reactivity out of the fresh fuel and pull the moderator
coefficient back towards negative. They are the next section.

## Burnable absorbers

A fresh core has too much reactivity. Soluble boron can hold it down, but
the more boron the water carries, the more positive the moderator
coefficient becomes (previous section). A burnable absorber is a solid
neutron absorber placed in the core for the first cycle. It holds down
part of the excess, so less boron is needed. As the cycle goes on, its
B-10 is consumed together with the U-235. (The burning itself needs
depletion, which this code does not have yet. Here the absorbers are as
loaded.)

**The Pyrex rod.** Watts Bar's first core uses Pyrex: borosilicate glass,
12.5 wt% B₂O₃ in SiO₂, density 2.25 g/cm³ [[VERA](#ref-vera) Table 6]. The
glass is a tube between two stainless-steel walls, with helium gaps on
both sides, and the rod sits in a guide tube (radii in cm):

| region | outer radius |
| --- | --- |
| helium (centre) | 0.214 |
| SS304 inner tube | 0.231 |
| helium | 0.241 |
| Pyrex | 0.427 |
| helium | 0.437 |
| SS304 clad | 0.484 |
| water | 0.561 |
| Zircaloy guide tube | 0.602 |

The number densities come from VERA Table P2-4. One check on them: the
B-10 in the glass annulus is
$N_{B10}\,\pi(0.427^2 - 0.241^2)\,m_{B10}$ = 6.25 mg per cm of rod.
VERA states 6.24 mg/cm. The hole in the middle keeps the boron where the
thermal flux is: B-10 absorbs so strongly that the glass would shield its
own centre.

The glass is the reason for the R-Matrix Limited format of chapter 6.
Natural silicon (Si-28, 29, 30) is evaluated only in LRF = 7. The same
work let the full SS304 of Table P2-4 into the rod walls and the baffle:
Fe-54 is LRF = 7 with only the gamma and neutron channels. Fe-57 (it has an
inelastic channel), Mn-55 and Cr-53 (negative channel spins) and Fe-58
(unresolved LRU = 2, LRF = 1) are still left out (`pincell.run.UNSUPPORTED`).

**Where the rods go.** A 17×17 assembly has 24 guide tubes. VERA Figure 5
lists the patterns that fill 4, 8, 12, 16, 20 or all 24 of them, keeping
the assembly's eighth symmetry. `models/vera-p2-pyrex.toml` builds each
pattern as a lattice `pyrex_N` and adds the two VERA states that test
them. Both are at 600 K with 1300 ppm boron. 2E uses Figure 5's 12-rod pattern
and 2F all 24.

**Against KENO.** The assembly solver is the one of chapter 19 (64
directions, 0.03 cm tracks, 3 pellet rings, SPH, coupled clad), run with
`make assembly`. The Pyrex rod gets no shielding of its own: its nuclides
are infinitely dilute, which suits B-10, whose absorption falls as 1/v
through the thermal and epithermal range. The rod's own self-shielding (fewer thermal neutrons
reach the inside of the glass) is done by MOC across its regions.

| state | Pyrex rods | $k$ (this code) | KENO | difference | pin powers, RMS / max |
| --- | --- | --- | --- | --- | --- |
| 2B | 0 | 1.18143 | 1.183360 | −193 pcm | 0.06% / 0.14% |
| 2E | 12 | 1.06927 | 1.069627 | −36 pcm | 0.08% / 0.18% |
| 2F | 24 | 0.97657 | 0.976018 | +56 pcm | 0.13% / 0.30% |

KENO from VERA Table P2-5, pin powers from Figures P2-7 and P2-8. The
offset in 2B comes from the fuel lattice (chapter 19). Measured against it,
$\rho_{2E} - \rho_{2B}$ is −8879 pcm here and −8985 pcm in KENO, and with
24 rods −17756 against −17952. (VERA has no Pyrex-free state at 2E's water
density, so the difference also contains the denser water, 0.743 against
0.661 g/cm³. At 1300 ppm that change is worth about +185 pcm in the pin
cell of chapter 20, small next to the rods.) The Pyrex is
**1.1–1.2% too weak**, consistently for both patterns. That is about 100 pcm too little
holddown for 12 rods. The pin powers stay close to KENO with the rods
in: the lowest pin falls from 0.941 (2B) to 0.931 (2E) and 0.908 (2F), and
no pin differs by more than 0.3%.

**State matrices.** The core needs constants for every lattice in it.
`models/wbn1-fuel-262.toml` and `wbn1-fuel-310.toml` now extend
`vera-p2-pyrex`. Each has one `[branches.pN]` set per Pyrex pattern of
VERA Figure 9: 16, 20 and 24 rods in 2.619% fuel, and 8, 12, 16 and 24
in 3.10% fuel. They use the same axes as the unpoisoned `[branches.core]`
(fuel temperature, moderator temperature, water density, 2×2×2 points
around 2B). The 2.11% fuel never carries Pyrex. `models/wbn1-pyrex.toml`
places them as in Figure 9 (1712 rods). It also splits each poisoned
assembly axially, because the glass covers 15.761 to 376.441 cm while the
fuel runs from 11.951 to 377.711 cm (Tables 4 and 6). Each assembly then has 3.81
cm of unpoisoned fuel at the bottom, 360.68 cm with Pyrex, and 1.27 cm
unpoisoned at the top. Each of the 56 points ran in 16–30 min (two
matrices in parallel on 16 cores).

**The rods turn the moderator coefficient negative.** At the base state
(lattice $k_\infty$, from `make branches`):

| lattice | moderator T (pcm/K) | density (pcm per 0.01 g/cm³) | total MTC (pcm/K) |
| --- | --- | --- | --- |
| 2.619%, no Pyrex | — | −48 | +10.5 |
| 2.619%, 16 rods | −4.08 | +20.0 | −8.7 |
| 2.619%, 20 rods | −4.80 | +36.3 | −13.3 |
| 2.619%, 24 rods | −5.39 | +51.9 | −17.5 |
| 3.10%, no Pyrex | — | −22 | +4.6 |
| 3.10%, 8 rods | −1.93 | +7.8 | −3.7 |
| 3.10%, 12 rods | −2.61 | +22.5 | −7.8 |
| 3.10%, 16 rods | −3.24 | +36.7 | −11.8 |
| 3.10%, 24 rods | −4.29 | +64.3 | −19.3 |

(The total is the moderator-temperature term plus the density term times
0.233, the IF97 density drop per kelvin, as in the previous section.)
The density coefficient changes sign, and more rods make it more
positive. Less dense water has two effects with opposite signs. It
moderates less, which costs reactivity. It also carries less boron, which
gains reactivity. Without Pyrex at 1300 ppm the boron term wins (the
previous section). The rods shift the balance in two ways. Each rod takes
the water out of its guide tube, so the lattice has less moderator and
sits closer to the undermoderated side. And the B-10 in the glass takes a
share of the thermal absorptions that does not leave when the water
expands, so the soluble boron's share, the term that gains reactivity, is
smaller. The lattice ends on the side where losing water loses reactivity.
That is what a PWR at power needs.

**The core.** `make feedback FEEDBACK_MODEL=wbn1-pyrex`:

| | zero power | 50% power | 100% power |
| --- | --- | --- | --- |
| $k$ | 0.99915 | 0.99461 | 0.98951 |
| power defect | | 457 pcm | 975 pcm |
| largest assembly power | 1.324 | 1.300 | 1.277 |
| axial offset (top − bottom) | −0.2% | −3.5% | −7.3% |
| Doppler temperature, average / max | 565 / 565 K | 677 / 784 K | 808 / 1039 K |
| peak fuel centreline | 565 K | 1029 K | 1640 K |
| water outlet | 565 K | 577–589 K | 588–609 K |
| feedback iterations | 3 | 7 | 20 |

The assembly powers at full power (quarter core, centre at top left; the
Pyrex map is in `models/wbn1-pyrex.toml`):

```text
1.058  0.977  1.062  1.023  1.140  1.052  1.063  0.772
0.977  1.051  0.936  1.112  1.063  1.155  1.032  0.872
1.062  0.936  1.093  1.049  1.175  1.127  1.079  0.783
1.023  1.112  1.049  1.171  1.088  1.140  1.019  0.639
1.140  1.063  1.175  1.088  1.277  0.904  0.915
1.052  1.155  1.127  1.140  0.904  0.913  0.608
1.063  1.032  1.079  1.019  0.915  0.608
0.772  0.872  0.783  0.639
```

* **Critical at 1300 ppm.** The rods are worth −7878 pcm in the core at
  zero power, and $k$ falls from 1.08451 to 0.99915. Our tables have no
  boron axis, so the code cannot search the critical boron itself. With the
  measured differential boron worth, −10.77 pcm/ppm [[VERA](#ref-vera)
  Table P5-7], the −85 pcm at 1300 ppm is made up by 8 ppm less boron: a
  predicted critical boron of 1292 ppm against the measured 1291. That is
  closer than the model deserves. It has no spacer grids, no nozzles, a
  smeared baffle and assembly constants that are −193 pcm off in 2B and
  1.1% weak in the Pyrex, and some of these errors cancel. (Before the
  reflector discontinuity factors of "Control rods" it was 1.00025, 11 ppm
  high.)
* **The isothermal temperature coefficient.** At zero power the start-up
  tests heat fuel and water together, from 560 to 570 K, and measure the
  reactivity change. `python -m fission.feedback wbn1-pyrex --itc` does
  the same (`feedback.core.isothermal`: one temperature everywhere, IF97
  density, fixed boron, NEM). It gives $k$ = 0.999420, 0.999146 and
  0.998837 at 560, 565 and 570 K: **−5.84 pcm/K = −3.24 pcm/°F**.
  Table P5-7 gives −2.17 pcm/°F measured and −3.18 ± 0.04 pcm/°F from CE
  KENO-VI. We agree with KENO within 2%, which is the closer comparison
  (same idealized core, no thermal expansion). 560 K and 0.752 g/cm³ are just outside
  the tables (565–600 K, 0.661–0.743 g/cm³), so that point is extrapolated
  linearly, and the reflector keeps its 565 K constants. Without
  Pyrex (`wbn1-zoned`) the same calculation gives **+3.43 pcm/°F**.
* **Flatter, and the peak moves.** The largest assembly power falls from
  1.482 to 1.324 at zero power. The peak is now the one 2.619% assembly
  without Pyrex (row 5, column 5 of the quarter map). Every other 2.619%
  assembly carries 16 to 24 rods and gives up power to the 2.11% fuel.
  The centre (2.11%) falls from 1.330 to 0.995, but its rodded 2.619%
  neighbours fall further (1.482 to 0.918), so the chess board is now
  inverted: the 2.11% assemblies make more power than the 2.619%
  assemblies around them.
* **The power goes to the bottom.** The axial offset goes from −0.2% to
  −7.3% at full power. Like the uniform core of "Results" and unlike the
  zoned core without Pyrex, denser water is now *more* reactive, so the
  cool bottom half takes more power. Doppler pushes the same way.
* **The power defect grows** from 711 to 975 pcm, because now the water
  adds to the Doppler effect instead of offsetting it.

## Control rods

Soluble boron changes slowly: the plant dilutes or borates the whole
primary circuit. Fast changes, and shutdown, are the job of the control
rods. A rod cluster control assembly (RCCA) is a spider of 24 rodlets that
enters the 24 guide tubes of one assembly from above. The rods move in
groups, called banks. Watts Bar's cycle 1 has four control banks A–D and
four shutdown banks SA–SD [[VERA](#ref-vera) §1.6, Figure 10].

**Two absorbers in one rod.** The Watts Bar rodlets are hybrids (Table 8).
The lower 101.6 cm is silver-indium-cadmium (Ag-In-Cd, 80/15/5 wt%,
10.2 g/cm³, radius 0.382 cm). The upper 259.08 cm is boron carbide (B₄C,
1.76 g/cm³, radius 0.373 cm). Both sit in an SS304 tube (0.386–0.484 cm)
with helium between, and fully inserted the poison starts 5.08 cm above
the bottom of the fuel. One step moves the rod 1.5875 cm, and 230 steps
take it out. The number densities of Table P2-4 give back the
specification exactly: 10.20 and 1.76 g/cm³, and 80.0/15.0/5.0 wt%.

The two absorbers capture in different ways. Boron absorbs through
B-10(n,α), 3844 b at 0.0253 eV and falling as $1/v$, with no resonances:
in B₄C the mean free path of a thermal neutron is 0.17 mm. Ag-In-Cd
absorbs mostly in resonances. Our 0 K tables give:

| nuclide | capture peak | at 0.0253 eV | in the rod |
| --- | --- | --- | --- |
| In-115 | 38 200 b at 1.457 eV | 202 b | mean free path 34 μm |
| Ag-109 | 31 100 b at 5.19 eV | 90 b | mean free path 15 μm |
| Ag-107 | 2 970 b at 16.3 eV | 37 b | |
| Cd-113 | 61 800 b at 0.174 eV | 19 900 b | |

At the peak of the Ag-109 resonance a neutron entering the rod travels
15 μm. The rod is 3.8 mm thick, so it is black there. A group constant
that averages the capture cross section with a smooth weighting spectrum
(infinite dilution) assumes the flux inside the rod is the same at the
peak as in the wings of the resonance. In reality the flux at the peak is
almost zero a few micrometres in. That is resonance self-shielding, the
problem of chapters 16, 17 and 19, now in an absorber instead of the fuel.

**Shielding the absorber.** The code already has the tools. For the fuel
it computes, by MOC, the escape probability $P(\Sigma)$ of the region for
any total cross section, and turns it into a flux shape with the
narrow-resonance balance (chapter 16):

$$
\phi(E) = P(\Sigma_t(E)) + \bigl(1 - P(\Sigma_t(E))\bigr)\,\frac{\Sigma_p}{\Sigma_t(E)}.
$$

The first term is the neutrons that arrive from outside, attenuated. The
second is the ones that scatter inside. For a black absorber the first
term carries almost everything. A material marked `shielded = true` in
the model file gets this treatment: its own escape probability table
from the assembly geometry (all 24 rodlets together), its own flux shape
and libraries, and SPH factors (chapter 17) so that the multigroup MOC
reproduces the shielded reaction rates. The fuel's slowing-down march is
not used here. The narrow-resonance balance treats the rod's own
scattering simply, and in a nearly black absorber the inflow term
dominates anyway.

**Against KENO.** VERA problem 2 has both rods in all 24 guide tubes, at
600 K and 0.743 g/cm³ (states 2G and 2H; KENO from Table P2-5, pin powers
from Figures P2-9 and P2-10). `models/vera-p2-rods.toml` holds them:

| state | $k$ (this code) | KENO | difference | pin powers, RMS / max |
| --- | --- | --- | --- | --- |
| 2G, Ag-In-Cd shielded | 0.84881 | 0.847695 | +111 pcm | 0.22% / 0.60% |
| 2G, Ag-In-Cd infinitely dilute | 0.82223 | 0.847695 | −2546 pcm | 0.74% / 1.41% |
| 2H, B₄C | 0.79064 | 0.788221 | +242 pcm | 0.19% / 0.49% |

Without shielding the rod absorbs far too much: 2.5% in $k$, and
$\rho_{2G} - \rho_{2B}$ = −36977 pcm, 11% more than KENO's. With it, measured from the 2B offset like the Pyrex,
$\rho_{2G} - \rho_{2B}$ is −33169 pcm against KENO's −33462 (0.9% weak), and
$\rho_{2H} - \rho_{2B}$ −41837 against −42363 (1.2% weak). The B₄C rod is
as weak as the Pyrex glass (1.1–1.2%). Both get their worth from B-10
(n,α), which suggests a common cause in the boron data or its thermal
treatment. KENO's reference uses ENDF/B-VII.0, this code ENDF/B-VIII.1.
The pin powers next to the rods agree to 0.5–0.6%. The infinitely dilute
rod depresses its neighbours too much: the largest error, 1.4%, is at
the assembly corner, the pin farthest from any rod, which takes the power
the rodded region loses.

**Rodded lattices for the core.** The core needs constants for fuel with
rods in, at every state the feedback can reach. Each fuel file
(`wbn1-fuel-{211,262,310}.toml`) gains `[branches.aic]` and
`[branches.b4c]`: the 24-rod lattice with the Ag-In-Cd or the B₄C
rodlets, on the same 2×2×2 axes as `[branches.core]`. That is 48 assembly
calculations. With two running side by side they took 25–30 min each for
Ag-In-Cd, whose SPH adds a fixed-source MOC per group, and about 20 min for
B₄C. The fuel temperature coefficient of the rodded 3.10% lattice is
−3.54 pcm/K, against −2.36 without the rods. The rods take most of the
thermal neutrons, so resonance capture in U-238 is a larger share of
what is left.

**Banks in the model.** `models/wbn1-rods.toml` extends `wbn1-pyrex`. It
gives each fuel lattice its rodded versions (`rodded = { aic = ..., b4c =
... }`), and it describes the rod and the banks as data (chapter 13). The
banks sit on the quarter map of VERA Figure 10:

```text
D . A . D . C .          quarter core, centre at top left
. . . . . SB . .
A . C . . . B .          57 RCCAs in the full core, all in
. . . A . SC . .         assemblies without Pyrex
D . . . D . SA
. SB . SD . .
C . B . SA
. . . .
```

A position at $s$ steps withdrawn has its absorber column start at
$26.58 + 1.5875\,s$ cm, measured from the bottom of the lower reflector.
The code cuts that assembly's axial stack at the edges of the Ag-In-Cd
and B₄C zones (`nodal.mesh.stacks`). The cut nodes take the rodded
lattice, and the nodal mesh is cut at the union of all the stacks
(chapter 13). Fully inserted, the Ag-In-Cd occupies the bottom 101.6 cm
of the fuel and the B₄C the rest.

**The core's power map against KENO.** VERA gives the KENO-VI assembly
powers of this core at its initial criticality: 1285 ppm, bank D at 167
steps, everything else out (Figure P5-6). `make feedback
FEEDBACK_MODEL=wbn1-rods` at zero power with `--bank D=167` (at 1300
ppm, the tables have no boron axis) gives the same map with a clear tilt.
The centre assemblies are 8% low, the edge ones up to 13% high, RMS 6.0%.
The tilt is the same at 2×2 nodes per assembly, so it is not the nodal
mesh. Nor is it the assembly constants, which match KENO within 0.2% in
every lattice (above, and chapter 19). What is left is the interface
between the fuel and the baffle reflector.

**Reflector discontinuity factors.** A homogenized node reproduces its
heterogeneous reaction rates and face currents only if the flux is
allowed to jump at its faces, $f_A\phi_A = f_B\phi_B$ (chapter 19
computes the factors $f$ of an assembly, chapter 21's NEM applies them).
Between two fuel assemblies the jump is small (ADFs of 0.996–1.0). At
the reflector it is not. Smith's equivalence theory [[Smith](#ref-smith1980)] gets the factors
from the same fuel | baffle | water slab that gives the reflector its
constants.

1. The MOC solution of the slab gives the scalar flux $\phi_s$ on the
   fuel-reflector face (`moc.face_flux`: one more sweep from the
   converged solution, each crossing contributing $W\psi/(\sin\theta
   |\cos\varphi|)$).
2. Each side, homogenized with its own MOC flux, is then solved as a
   homogeneous two-group diffusion problem with the far side reflective,
   as in the slab:
   $$
   -D\phi'' + \Bigl(A - \tfrac{1}{k}F\Bigr)\phi = 0,
   $$
   with the removal matrix $A$ (absorption, out-scattering, minus
   in-scattering) and the fission matrix $F_{gh} = \chi_g\nu\Sigma_{f,h}$,
   zero in the reflector. Write
   $D^{-1}(A - F/k) = V\,\mathrm{diag}(\kappa_i^2)\,V^{-1}$. The solution is a sum of
   $\cosh\kappa_i(W-x)$ modes, and on the face
   $$
   \phi(0) = V\,\mathrm{diag}\bigl(\kappa_i W \coth \kappa_i W\bigr)\,V^{-1}\,\bar\phi,
   $$
   where $\bar\phi$ is the side's average flux. Fixing $\bar\phi$ is the
   same as fixing the net current through the face: by the side's
   balance, $J = (A - F/k)\bar\phi\,W$. In the fuel the fundamental mode
   has $\kappa^2 < 0$, and the cosh becomes a cosine; complex arithmetic
   covers both (`fewgroup.reflector.homogeneous_surface_flux`, tested
   against fine-mesh finite differences).
3. Each side's factor is $f = \phi_s / \phi(0)$. The nodal interface
   condition sees only the ratio $f_{\mathrm{refl}}/f_{\mathrm{fuel}}$, and in the
   core the fuel keeps the ADF of its own infinite lattice (about 1). So
   the reflector table gets $f_{\mathrm{refl}}/f_{\mathrm{fuel}}$.

The slab gives:

| | fast | thermal |
| --- | --- | --- |
| reflector side, $\phi_s/\phi_{\mathrm{refl}}(0)$ | 1.153 | 0.248 |
| fuel side, $\phi_s/\phi_{\mathrm{fuel}}(W)$ | 0.909 | 0.962 |
| reflector ADF | 1.269 | 0.258 |

The thermal factor is far from 1, and for a clear reason. In the slab the
thermal flux on the face is low: the steel next to the fuel absorbs, and
the thermal flux peaks several centimetres out, in the water behind the
baffle. The homogenized reflector smears steel and water together and
puts its thermal flux where its source is, at the face. Without the
factor, continuity of flux at the face makes the reflector send far too
many thermal neutrons back into the fuel. The fuel side matters too. A
first version took only the reflector side, $f = 1.153 / 0.248$, and
overcorrected: the centre went 3% high and the edge 9% low, RMS 3.0%.

With both sides, the assembly powers against KENO (per cent; Figure P5-6,
octant, centre at top left):

```text
 -1.4
 -1.9  -1.1
 -1.0  -2.2  -0.4
 -1.1  -0.4  -1.1  +0.5
 +0.1  -0.7  +0.5  +0.0  +0.7
 +0.0  +0.8  +0.1  +1.3  +1.1  +1.4
 +1.2  +0.2  +1.7  +0.8  -0.4  -3.0
 +0.2  -0.0  +0.5  -0.8
```

RMS 1.1%, largest −3.0%, from 6.0% and +13%. The eigenvalue at this state
is 0.99815 at 1300 ppm. With the measured boron worth (−10.77 pcm/ppm) it
becomes about −23 pcm at the measured 1285 ppm. All the core results of
this chapter's baffle sections are with these factors (the uniform core
of "Results" has a water reflector and keeps the older constants without
them).

**Bank worths.** `python -m fission.feedback wbn1-rods --bank-worths`
inserts each bank fully, alone, at zero power. This is how VERA's KENO
reference is defined: each bank with the others withdrawn
[[VERA](#ref-vera) Problem 5, "Reference solution results and measured data"]. The measured worths come
from rod swap against bank D, inferred with KENO shadow factors, at
1170 ppm. VERA calls the boron a secondary effect (2 pcm on bank D):

| bank | this code (pcm) | KENO | difference | measured | difference | no reflector ADF |
| --- | --- | --- | --- | --- | --- | --- |
| D | 1368 | 1386 | −1.3% | 1342 | +1.9% | 1320 (−4.8%) |
| A | 861 | 898 | −4.1% | 843 | +2.1% | 797 (−11.3%) |
| C | 958 | 984 | −2.6% | 951 | +0.8% | 911 (−7.4%) |
| SB | 1055 | 1066 | −1.0% | 1056 | −0.1% | 1033 (−3.1%) |
| B | 897 | 875 | +2.5% | 879 | +2.0% | 940 (+7.4%) |
| SC | 466 | 499 | −6.6% | 480 | −2.9% | 477 (−4.5%) |
| SA | 469 | 447 | +5.0% | 435 | +7.9% | 530 (+18.6%) |
| SD | 466 | 499 | −6.6% | 480 | −2.9% | 477 (−4.5%) |
| total | 6540 | 6654 | −1.7% | 6467 | +1.1% | 6484 (−2.6%) |

(KENO and measured values from Table P5-7; the last column against KENO.)

* **The reflector made the pattern.** Without the factors, central banks
  (A, C, D) were too weak and peripheral ones (B, SA) too strong. A
  bank's worth grows roughly with the square of the flux where it sits
  (first-order perturbation theory weights an absorber by $\phi\phi^\dagger$),
  so a radial tilt shows up twice in it. With the factors the spread
  falls from −11..+19% to −7..+5% against KENO.
* **The total is 1.7% weak against KENO**, close to the 0.9% (Ag-In-Cd)
  and 1.2% (B₄C) of the single assembly. Most of the rod's worth comes from
  the lower Ag-In-Cd zone at zero power, because the axial flux peaks
  mid-core and the Ag-In-Cd covers the bottom 102 cm.
* **SC and SD.** In the real core these banks sit at positions that are
  related by a quarter rotation, not by reflection. A reflective quarter
  core cannot reproduce that exactly. Here SC and SD come out equal (466
  pcm). KENO, on the full core, also has them equal (499). They are 6.6%
  weak, more than the other banks, perhaps from this approximation.
* **Against the plant** the total is within 1.1%, and every bank but SA
  is within 3%. KENO itself is 2.9% high on the total. Part of our
  apparent agreement with the measurement is the 1.7% our rods are weak.

## Further reading

* [[IAPWS-IF97](#ref-if97)] (free): the water formulation, region 1.
* [[MOX benchmark](#ref-moxbench) §7] (free): the simplified fuel-rod
  properties.
* [[IAPWS viscosity](#ref-iapwsvisc)], [[IAPWS conductivity](#ref-iapwscond)]
  (free): transport properties of water.
* [[MOOSE](#ref-moose)] (free): the Dittus–Boelter correlation and Presser's
  bundle factor, with their ranges.
* [[VERA](#ref-vera)] (free): Watts Bar core data (Figure 9, Tables 1, 18,
  P6-1); the Pyrex rod and its patterns (Table 6, Figure 5, Table P2-4),
  KENO references (Table P2-5, Figures P2-7, P2-8) and the start-up tests
  (Table P5-7); the control rods (Table 8, Figure 10), KENO for 2G/2H
  (Figures P2-9, P2-10) and the core's power map (Figure P5-6).
* [[Smith](#ref-smith1980)] (free to read): discontinuity factors, the
  equivalence behind the reflector factors.
* Code: [`feedback/water.py`](../../src/fission/feedback/water.py),
  [`feedback/convection.py`](../../src/fission/feedback/convection.py),
  [`feedback/fuel.py`](../../src/fission/feedback/fuel.py),
  [`feedback/core.py`](../../src/fission/feedback/core.py) (`make feedback`),
  [`fewgroup/reflector.py`](../../src/fission/fewgroup/reflector.py)
  (`make reflector`), `models/wbn1-core.toml`, `models/wbn1-zoned.toml`
  (with `wbn1-fuel-*.toml` and `wbn1-baffle.toml`), `models/wbn1-pyrex.toml`,
  `models/vera-p2-rods.toml`, `models/wbn1-rods.toml`;
  [`nodal/mesh.py`](../../src/fission/nodal/mesh.py) (`stacks`: the banks
  cut the axial stacks).
