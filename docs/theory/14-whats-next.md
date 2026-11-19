# What comes next

Chapters 7 and 8 give pointwise tables at any temperature, with averaged
("infinitely dilute") cross sections in the unresolved range; chapter 9 adds
the data on the neutrons that reactions emit, chapter 10 condenses
everything into infinitely dilute multigroup constants, chapter 11 adds
thermal scattering, chapter 12 shields them in a homogeneous mixture, and
chapter 13 describes the reactor model they are applied to. In
a real fuel pin, though, fuel and moderator are separated, and the flux has
spatial as well as energy structure. This chapter previews the next processing
steps; each will get its own chapter once implemented.

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
pin-cell model comes in. Section 8.7 already used this flux to show the
Doppler effect: the shielded U-238 capture grows with temperature while the
dilute value does not.

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
