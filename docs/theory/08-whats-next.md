# What comes next

The tables of chapter 7 describe a nucleus at rest (0 K) in the resolved
range, and averaged ("infinitely dilute") cross sections in the unresolved
range. Real fuel is hot, and in a fuel pin the neutron flux itself has
resonance structure. This chapter previews the next processing steps; each
will get its own chapter once implemented.

## Doppler broadening

Target nuclei are in thermal motion, so a neutron of laboratory energy $E$
meets them at a spread of *relative* energies. The cross section an
experiment or a reactor sees is the 0 K cross section averaged over that
spread. For a free-gas (Maxwellian) target and a linearly interpolable 0 K
table, the result is exact in closed form [[Cullen & Weisbin](#ref-cullen);
[NJOY](#ref-njoy) ch. 4]. With $\alpha = \text{AWR}/(k_B T)$, $y^2 = \alpha E$
and $x^2 = \alpha E'$:

$$
\sigma(y, T) = \frac{1}{y^2\sqrt{\pi}}\int_0^\infty x^2\,\sigma(x, 0)
\left[e^{-(x-y)^2} - e^{-(x+y)^2}\right]\mathrm{d}x .
$$

The scale of the effect is the **Doppler width**

$$
\Gamma_D = 2\sqrt{\frac{E\,k_B T}{\text{AWR}}} .
$$

For the 6.67 eV U-238 resonance in fuel at 900 K, $\Gamma_D = 0.094$ eV, four
times the natural width $\Gamma = 0.024$ eV. The peak is lowered and the wings
are raised. For a single resonance with $\xi = \Gamma/\Gamma_D$, the peak
falls by the factor
$\psi(0, \xi) = \tfrac{\sqrt\pi}{2}\,\xi\, e^{\xi^2/4}\operatorname{erfc}(\xi/2)
\approx 0.20$, from 22,765 b to about 4,500 b. The *area* stays essentially
the same.

Doppler broadening does not change the infinitely dilute resonance integral,
but it changes the *self-shielded* one. This is the origin of the **Doppler
coefficient**: hotter fuel captures more neutrons in U-238. It is the most
important prompt negative feedback in a light-water reactor
[[Bell & Glasstone](#ref-bell) §8.1d, §8.3e; [Fröhner](#ref-frohner) §3.5].

## Self-shielding

In a fuel pin the flux is not smooth across a resonance. Where U-238 captures
strongly, few neutrons survive, and the flux dips, so the reaction rate is
much smaller than $\sigma \times$ (average flux) suggests. The quantity a
group calculation needs is the flux-weighted **effective cross section**

$$
\sigma_{\text{eff}} = \frac{\int \sigma(E)\,\phi(E)\,\mathrm{d}E}{\int \phi(E)\,\mathrm{d}E} .
$$

The classic estimate of the flux is the **narrow-resonance approximation**
[[Bell & Glasstone](#ref-bell) §8.3]: for a resonance narrow compared with the
energy lost per collision,

$$
\phi(E) \propto \frac{1}{E\,\left[\sigma_t(E) + \sigma_0\right]},
$$

where $\sigma_0$ is the **background cross section** per absorber atom (all
other scattering in the mixture, plus a geometric term for the escape of
neutrons from the fuel). For $\sigma_0 \to \infty$ the flux is $1/E$ and
$\sigma_{\text{eff}}$ is the infinitely dilute value of chapter 7. For small
$\sigma_0$, the big resonances are strongly shielded. Tables of
$\sigma_{\text{eff}}/\sigma_\infty$ versus $\sigma_0$ and temperature are
**Bondarenko factors** [[Fröhner](#ref-frohner) §4.3.1]. The geometric part
of $\sigma_0$ (Dancoff factors, escape probabilities, equivalence relations
[[Bell & Glasstone](#ref-bell) §8.4]) is where the transport solver of the
pin-cell model comes in.

## The unresolved range

Between ~20 keV and ~150 keV (U-238), only average resonance parameters are
known: mean level spacing, average widths, and the number of degrees of
freedom of their Porter–Thomas-type distributions [[ENDF-102](#ref-endf102)
§D.2]. Averaged cross sections follow from statistical (Hauser–Feshbach type)
formulas [[Fröhner](#ref-frohner) §4.2]. Self-shielding there needs either
Bondarenko factors (NJOY's UNRESR, [[NJOY](#ref-njoy) ch. 5]) or
**probability tables**, distributions of the cross section sampled from
statistical resonance ladders (NJOY's PURR, ch. 23). The unresolved range
contributes little to U-238's dilute resonance integral (0.9 b of 275 b), but
it matters for fast spectra and for the Doppler coefficient.

## Thermal scattering

Below a few eV, the hydrogen in water is not a free atom. Chemical binding,
molecular vibration and rotation change how neutrons exchange energy with it.
The free-atom H-1 cross section of chapter 1 (20.4 b) is replaced by the
**thermal scattering law** $S(\alpha, \beta)$ from the thermal sublibrary
(MF7, e.g. `tsl-HinH2O`) [[ENDF-102](#ref-endf102) §7; [NJOY](#ref-njoy)
ch. 7 THERMR, ch. 24 LEAPR]. It sets the shape of the thermal neutron
spectrum, and with it every thermal reaction rate.

## Further reading

* [[NJOY](#ref-njoy)] ch. 4 BROADR, ch. 5 UNRESR, ch. 7 THERMR, ch. 23 PURR
  (free).
* [[Fröhner](#ref-frohner)] §3.5 (Doppler broadening, including point data),
  §4 (unresolved range, Bondarenko factors) (free).
* [[Bell & Glasstone](#ref-bell)] ch. 7 (neutron thermalization) and ch. 8
  (resonance absorption: Doppler broadening §8.1d, unresolved parameters §8.2,
  the narrow-resonance approximation §8.3) (free).
* [[OpenMC](#ref-openmc)] methods documentation: how a modern Monte Carlo
  code represents and broadens cross sections (free).
* [[Cullen & Weisbin](#ref-cullen)], the exact broadening of tabulated data
  (not free); [[Hébert](#ref-hebert)], resonance self-shielding in lattice
  codes (not free).
