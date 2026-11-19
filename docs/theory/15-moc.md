# The method of characteristics

Chapter 14 divided the geometry into flat-source regions and covered it with
cyclic tracks. This chapter solves the multigroup transport equation on them,
for the multiplication factor $k$ and the flux in every region. It then checks
the solver on the standard benchmark for this kind of code, C5G7.

## The transport equation along a track

For energy group $g$, the angular flux $\psi_g$ along a straight flight path
obeys

$$
\frac{\mathrm{d}\psi_g}{\mathrm{d}s} + \Sigma_{t,g}\,\psi_g = Q_g ,
\qquad
Q_{i,g} = \frac{1}{4\pi}\Big[\sum_{g'} \Sigma_{s,g'\to g}\,\phi_{i,g'}
+ \frac{\chi_{i,g}}{k}\sum_{g'}\nu\Sigma_{f,i,g'}\,\phi_{i,g'}\Big].
$$

Scattering is taken as **isotropic**, and $Q_{i,g}$ is constant in each flat-source
region $i$. Anisotropic scattering is approximated with the **transport
correction**: the total cross section is replaced by the transport cross
section $\Sigma_{tr} = \Sigma_t - \bar\mu\,\Sigma_s$, which keeps the
diffusion length right [[Bell & Glasstone](#ref-bell)]. The C5G7 data
are given this way.

A 2D track carries neutrons moving at every **polar angle** $\theta$ (measured
from the axial direction $z$). Such a neutron travels $\ell/\sin\theta$ while
crossing a segment of 2D length $\ell$. With
$\tau = \Sigma_{t,g}\,\ell/\sin\theta$, the solution of chapter 14 gives the
change of angular flux over one segment:

$$
\Delta\psi = \Big(\psi_{\text{in}} - \frac{Q}{\Sigma_t}\Big)\left(1 - e^{-\tau}\right),
\qquad
\psi_{\text{out}} = \psi_{\text{in}} - \Delta\psi .
$$

$1 - e^{-\tau}$ is computed as $-\mathrm{expm1}(-\tau)$, which stays accurate
for the very thin segments that tracks cut near corners.

## From segments to region fluxes (derivation)

Integrate the transport equation over one segment. The neutrons that leave
it, minus those that enter, equal the source inside minus the collisions:
$\psi_{\text{out}} - \psi_{\text{in}} = \frac{\ell}{\sin\theta}\,Q -
\Sigma_t\int\psi\,\mathrm d s$. So the integral of $\psi$ along the segment's
3D path is $\big(\frac{\ell}{\sin\theta}Q + \Delta\psi\big)/\Sigma_t$.
Measured along the 2D track, multiply by $\sin\theta$:

$$
\int_{\text{segment}} \psi \,\mathrm d s_{2D}
= \frac{\ell\,Q + \sin\theta\,\Delta\psi}{\Sigma_t}.
$$

A track stands for a strip of width $\delta$. Summing $\delta\times$ this
integral over all segments in region $i$ integrates $\psi$ over the region's
area. Since $\sum\delta\,\ell = A_i$ (chapter 14), the average angular flux
in direction $(\varphi, \theta)$ is

$$
\bar\psi_i = \frac{Q_i}{\Sigma_{t,i}}
+ \frac{1}{\Sigma_{t,i} A_i}\sum_{\text{segments in } i}\delta\,\sin\theta\,\Delta\psi .
$$

The scalar flux integrates over all directions, $\phi = \int_{4\pi}\psi\,\mathrm
d\Omega$. In 2D the flux is symmetric above and below the plane, so the
polar angles on $(0, \pi/2)$ are counted twice. Each track is swept both ways,
which covers the azimuths in $(0,\pi)$ and $(\pi, 2\pi)$. With azimuthal
weights $w_a$ (summing to 1 over $(0,\pi)$, chapter 14) and polar weights
$\omega_p$ (summing to 1 over $(0,\pi/2)$), each sweep direction carries the
solid angle $2\pi\,w_a\,\omega_p$, and all directions together add up to
$4\pi$. Therefore

$$
\phi_{i,g} = \frac{4\pi\,Q_{i,g}}{\Sigma_{t,i,g}}
+ \frac{1}{\Sigma_{t,i,g}\,A_i}\sum_{\text{segments in } i} W\,\Delta\psi_g ,
\qquad
W = 2\pi\,w_a\,\omega_p\,\delta_a\,\sin\theta_p .
$$

The first term is the flux of an infinite medium with that source. The
second term corrects it with what the tracks see: neutrons streaming in
from elsewhere, or escaping.

**Polar quadrature.** This code uses Gauss–Legendre points in
$\mu = \cos\theta$ on $(0, 1)$. Production MOC codes often use the
Tabuchi–Yamamoto set, which is optimized for the integrals that appear here
and is more accurate with two or three angles. Gauss–Legendre is easy to
derive and converges as the number of angles grows.

**Boundaries.** A sweep's outgoing flux becomes, in the next iteration, the
incoming flux of the sweep it is linked to (reflective sides). At a vacuum
side it is counted as **leakage**, $\sum W\,\psi_{\text{out}}$.

**Area correction.** Before solving, every segment length in region $i$ is
multiplied by $A_i/\tilde A_i$, with $\tilde A_i$ the area estimated from the
tracks (chapter 14). This makes the weighted track areas exact.

## Power iteration

$k$ is the largest eigenvalue of the transport problem. The solver uses
**power iteration** with a single transport sweep per iteration:

1. From the current $\phi$ and $k$, form the source $Q$ (scattering plus
   fission over $k$).
2. Sweep every track in both directions, and form the new $\phi$.
3. Update $k \leftarrow k\cdot\sum_i A_i\,\nu\Sigma_f\phi^{\text{new}} /
   \sum_i A_i\,\nu\Sigma_f\phi^{\text{old}}$.

At convergence, neutron production divided by (absorption + leakage) equals
$k$; the tests check this. "Absorption" here must be
$\Sigma_t - \sum_{g'}\Sigma_{s,g\to g'}$, the removal the transport equation
actually uses. The benchmark's tabulated absorption agrees with it only to
the six printed digits.

**Slow convergence.** One sweep moves angular flux one track length across
the geometry. In regions that are thin and scatter strongly, the incoming
fluxes approach their final values only slowly. Successive changes shrink
by a nearly constant ratio $\rho$, close to 1 (0.995 per iteration in the
one-material test). A small change per iteration therefore does not mean
convergence: after a change $d$, the remaining error is
$d(\rho + \rho^2 + \dots) = d\,\rho/(1-\rho)$. The solver measures $\rho$ over
the last 10 iterations and stops when this *estimated remaining error* of
$k$ and of the fission source is below the tolerance. The cure for the
slowness is acceleration: a coarse-mesh diffusion problem that moves
neutrons across the whole core at once (CMFD), which is step B3.

**Rust.** The sweep is the hot loop. Its Rust port (`moc.rs`) splits the
tracks over threads: within one sweep, tracks are independent, because their
incoming fluxes come from the previous iteration. Each thread keeps its own
tally, and the tallies are added at the end. The iteration around it stays in
Python. On a pin cell, Python and Rust agree to $10^{-15}$ in $k$; one C5G7
sweep (5.9 million segments, 7 groups, 3 polar angles) takes 0.17 s.

## Checks

* **One material everywhere.** If all regions hold the same material, the
  flat flux is exact whatever the tracking. MOC must then return the
  infinite-medium $k_\infty$ of that material (chapter 12). It does, to
  $10^{-9}$, with the flux flat to $10^{-9}$.
* **Neutron balance**, with reflective and with vacuum sides: production
  over absorption plus leakage equals $k$ to $10^{-7}$.
* **Discretization.** Refining the tracks (8→64 directions, 0.1→0.01 cm)
  converges $k$ of a pin cell.

## C5G7

The C5G7 benchmark [[C5G7](#ref-c5g7)] is a quarter of a small PWR core:

* two UO$_2$ and two MOX assemblies of 17×17 pins, surrounded by water;
* reflective sides on the top and left, vacuum on the bottom and right;
* seven-group, transport-corrected cross sections for seven materials:
  UO$_2$, three MOX enrichments, guide tube, fission chamber and moderator.

The pins are cylinders of radius 0.54 cm, with no separate clad: fuel and
clad are mixed into one material, so there is no fuel-coolant
homogenization. The reference is a multigroup Monte Carlo (MCNP) solution
with the same cross sections: $k = 1.18655$, to 0.003%. Any difference is
therefore an error of the deterministic method, not of the data.

The model file `models/c5g7.toml` holds the cross sections (Table 1 of the
specification, extracted from its text by script) and the layout. The
assembly maps were read from the specification's Figure 3 by classifying
the grey level of each cell. They match the stated counts: 24 guide tubes
and one fission chamber per assembly, and MOX zones of 64, 100 and 100
pins.

Results, each converged to an estimated error below $10^{-6}$ in $k$ and
$10^{-5}$ in the fission source:

| tracking | regions | segments | sweeps | time | $k$ | vs MCNP |
| --- | --- | --- | --- | --- | --- | --- |
| 16 directions, 0.1 cm, 2 polar, 1 ring, 4 sectors | 15,028 | 0.87 M | 1241 | 82 s | 1.185432 | −112 pcm |
| 32 directions, 0.05 cm, 3 polar, 2 rings, 8 sectors | 39,304 | 5.9 M | 1242 | 314 s | 1.186933 | +38 pcm |
| as above, 6 polar | 39,304 | 5.9 M | 1242 | 705 s | 1.186981 | +43 pcm |

| | reference (MCNP) | coarse | fine (3 polar) |
| --- | --- | --- | --- |
| inner UO$_2$ assembly power | 492.8 | 490.0 (−0.58%) | 491.1 (−0.34%) |
| MOX assembly power | 211.7 | 212.5 (+0.36%) | 212.2 (+0.22%) |
| outer UO$_2$ assembly power | 139.8 | 141.1 (+0.96%) | 140.5 (+0.53%) |
| maximum pin power | 2.498 | 2.482 (−0.63%) | 2.488 (−0.40%) |

The pin powers are normalized to an average of 1 over the 1056 fuel pins,
and an assembly's power is the sum over its 264 pins. Doubling the polar
angles changes $k$ by only 5 pcm. The remaining difference therefore comes
mostly from the spatial and azimuthal discretization: flat sources in
large regions (the reflector cells are 1.26 cm squares cut only into
sectors) and the track density. For comparison, the benchmark's
participants differ from the reference by 0.1 to about 950 pcm
[[C5G7](#ref-c5g7), Table 17]. Among the characteristics codes, DeCART
is at +4 pcm and APOLLO2 at −31 pcm. Refining the regions of this
calculation, especially in the reflector, is the next lever, and a cheap
one once CMFD cuts the sweep count.

Each run needed about 1240 sweeps. That number depends on the core, not on
the tracking, and it is what CMFD acceleration (B3) is meant to cut.

![C5G7 with coarse tracking (16 directions, 0.1 cm, 2 polar angles). Left: pin powers of the four fuel assemblies. Right: how slowly power iteration approaches its final eigenvalue.](../figures/c5g7.svg)

## Further reading

* [[C5G7](#ref-c5g7)] (free): the specification (Chapter 2, Appendix A) and
  the results of the participating codes.
* [[OpenMOC](#ref-openmoc)], [[OpenMOC docs](#ref-openmocdocs)] (free):
  the same equations, the Tabuchi–Yamamoto quadrature, and CMFD.
* [[Hébert](#ref-hebert)] (not free): the method of characteristics and its
  acceleration.
* Code: [`moc.py`](../../src/fission/transport/moc.py) (sweep, power
  iteration, convergence estimate), [`moc.rs`](../../rust/src/moc.rs),
  [`run.py`](../../src/fission/transport/run.py) (C5G7 driver and reference
  values), [`c5g7.toml`](../../models/c5g7.toml),
  [`c5g7-xs.toml`](../../models/c5g7-xs.toml); `make moc`.
