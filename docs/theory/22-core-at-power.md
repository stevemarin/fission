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
water-only constants.

**Results.** All the constants are in place, so `make feedback
FEEDBACK_MODEL=wbn1-zoned` runs at any power:

| | zero power | 50% power | 100% power |
| --- | --- | --- | --- |
| $k$ | 1.08495 | 1.08115 | 1.07678 |
| power defect | | 323 pcm | 699 pcm |
| largest assembly power | 1.436 | 1.414 | 1.387 |
| axial offset (top − bottom) | 0.0% | +1.1% | 0.0% |
| Doppler temperature, average / max | 565 / 565 K | 677 / 806 K | 809 / 1084 K |
| peak fuel centreline | 565 K | 1079 K | 1748 K |
| water outlet | 565 K | 574–591 K | 583–612 K |
| feedback iterations | 3 | 7 | 10 |

The assembly powers at full power (quarter core, centre at top left):

```text
1.250  1.387  1.220  1.319  1.119  1.155  0.882  0.667
1.387  1.235  1.357  1.175  1.246  1.025  1.059  0.655
1.220  1.357  1.196  1.300  1.101  1.120  0.823  0.598
1.319  1.175  1.300  1.143  1.223  0.972  0.903  0.465
1.119  1.246  1.101  1.223  1.181  1.023  0.734
1.155  1.025  1.120  0.972  1.023  0.893  0.489
0.882  1.059  0.823  0.903  0.734  0.489
0.667  0.655  0.598  0.465
```

What the loading pattern does:

* **The peak falls from 2.18 to 1.44** at zero power. The centre no longer
  dominates, and the chess board shows: every 2.619% assembly has more
  power than its 2.11% neighbours (1.436 against 1.288 and 1.267 next to the
  centre). The hottest channel now heats its water to 612 K, so the core
  runs at full power below the 623.15 K limit of IF97 region 1. (The code
  does not test for local boiling.)
* **The edge still loses.** Even with 3.10% fuel, the outer assemblies make
  0.47–0.67 of the average. The baffle reflects fewer thermal neutrons than
  water would.
* **The axial shape barely moves with power.** In the uniform core of the
  previous section the power slid to the bottom (−7.5% at half power),
  because cool, dense water there moderated better. Here denser water is
  *less* reactive (overmoderated), which pushes the power up. The Doppler
  effect pushes it down, because the fuel at the top is hotter: it sits in
  hotter water. The two nearly cancel.
* **The power defect is 699 pcm at full power** and comes from the Doppler
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
the rods' worth $\rho_{2E} - \rho_{2B}$ is −8879 pcm here and −8985 pcm in
KENO, and with 24 rods −17756 against −17952. The Pyrex is **1.1–1.2% too
weak**, consistently for both patterns. That is about 100 pcm too little
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
| $k$ | 1.00025 | 0.99572 | 0.99067 |
| power defect | | 454 pcm | 967 pcm |
| largest assembly power | 1.322 | 1.296 | 1.271 |
| axial offset (top − bottom) | −0.2% | −3.6% | −7.4% |
| Doppler temperature, average / max | 565 / 565 K | 677 / 783 K | 808 / 1036 K |
| peak fuel centreline | 565 K | 1027 K | 1634 K |
| water outlet | 565 K | 578–589 K | 591–609 K |
| feedback iterations | 3 | 7 | 19 |

The assembly powers at full power (quarter core, centre at top left; the
Pyrex map is in `models/wbn1-pyrex.toml`):

```text
1.005  0.929  1.014  0.982  1.105  1.031  1.054  0.807
0.929  1.000  0.895  1.070  1.032  1.133  1.026  0.918
1.014  0.895  1.050  1.014  1.147  1.113  1.081  0.827
0.982  1.070  1.014  1.140  1.070  1.137  1.040  0.722
1.105  1.032  1.147  1.070  1.271  0.915  0.993
1.031  1.133  1.113  1.137  0.915  0.946  0.699
1.054  1.026  1.081  1.040  0.993  0.699
0.807  0.918  0.827  0.722
```

* **Critical at 1300 ppm.** The rods are worth −7805 pcm in the core at
  zero power, and $k$ falls from 1.08495 to 1.00025. Our tables have no
  boron axis, so the code cannot search the critical boron itself. With the
  measured differential boron worth, −10.77 pcm/ppm [[VERA](#ref-vera)
  Table P5-7], the +25 pcm at 1300 ppm becomes +122 pcm at the measured
  1291 ppm: a predicted critical boron about 11 ppm too high. That is
  closer than the model deserves. It has no spacer grids, no nozzles, a
  smeared baffle and assembly constants that are −193 pcm off in 2B and
  1.1% weak in the Pyrex, and some of these errors cancel.
* **The isothermal temperature coefficient.** At zero power the start-up
  tests heat fuel and water together, from 560 to 570 K, and measure the
  reactivity change. `python -m fission.feedback wbn1-pyrex --itc` does
  the same (`feedback.core.isothermal`: one temperature everywhere, IF97
  density, fixed boron, NEM). It gives $k$ = 1.000528, 1.000249 and
  0.999934 at 560, 565 and 570 K: **−5.94 pcm/K = −3.30 pcm/°F**.
  Table P5-7 gives −2.17 pcm/°F measured and −3.18 ± 0.04 pcm/°F from CE
  KENO-VI. We agree with KENO within 4%, which is the closer comparison
  (same idealized core, no thermal expansion). 560 K and 0.752 g/cm³ are just outside
  the tables (565–600 K, 0.661–0.743 g/cm³), so that point is extrapolated
  linearly, and the reflector keeps its 565 K constants. Without
  Pyrex (`wbn1-zoned`) the same calculation gives **+3.42 pcm/°F**.
* **Flatter, and the peak moves.** The largest assembly power falls from
  1.436 to 1.322 at zero power. The peak is now the one 2.619% assembly
  without Pyrex (row 5, column 5 of the quarter map). Every other 2.619%
  assembly carries 16 to 24 rods and gives up power to the 2.11% fuel.
  The centre (2.11%) falls from 1.288 to 0.927, but its rodded 2.619%
  neighbours fall further (1.436 to 0.857), so the chess board is now
  inverted: the 2.11% assemblies make more power than the 2.619%
  assemblies around them.
* **The power goes to the bottom.** The axial offset goes from −0.2% to
  −7.4% at full power. Like the uniform core of "Results" and unlike the
  zoned core without Pyrex, denser water is now *more* reactive, so the
  cool bottom half takes more power. Doppler pushes the same way.
* **The power defect grows** from 699 to 967 pcm, because now the water
  adds to the Doppler effect instead of offsetting it.

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
  (Table P5-7).
* Code: [`feedback/water.py`](../../src/fission/feedback/water.py),
  [`feedback/convection.py`](../../src/fission/feedback/convection.py),
  [`feedback/fuel.py`](../../src/fission/feedback/fuel.py),
  [`feedback/core.py`](../../src/fission/feedback/core.py) (`make feedback`),
  [`fewgroup/reflector.py`](../../src/fission/fewgroup/reflector.py)
  (`make reflector`), `models/wbn1-core.toml`, `models/wbn1-zoned.toml`
  (with `wbn1-fuel-*.toml` and `wbn1-baffle.toml`).
