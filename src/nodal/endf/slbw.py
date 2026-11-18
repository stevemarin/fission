"""Single-Level Breit-Wigner (SLBW) cross sections from resolved resonance parameters.

SLBW treats every resonance as if it were alone: the cross section is a sum of
independent Lorentzian-like peaks (ENDF-102, Appendix D.1.1). For one resonance
r with spin J and orbital momentum l:

    capture   sigma_g = (pi/k^2) g_J  Gn Gg / ((E - E'r)^2 + G^2/4)
    fission   sigma_f = (pi/k^2) g_J  Gn Gf / ((E - E'r)^2 + G^2/4)
    elastic   sigma_n = (4 pi/k^2)(2l+1) sin^2(phi_l)                       potential
                      + (pi/k^2) g_J [Gn^2 - 2 Gn G sin^2(phi_l)
                                      + 2 (E - E'r) Gn sin(2 phi_l)]
                                    / ((E - E'r)^2 + G^2/4)                   resonance

with G = Gn + Gg + Gf the total width. The ingredients:

  * pi/k^2   the neutron's "size" (reduced wavelength squared): huge at low E,
             which is why slow-neutron cross sections can be thousands of barns.
  * g_J      the share of spin orientations able to form the compound state.
  * Gn(E)    the neutron width varies with energy as the penetrability:
             Gn(E) = Gn(|Er|) P_l(E) / P_l(|Er|). For s-waves P ~ k ~ sqrt(E), so
             far below a resonance capture goes like (1/k^2) sqrt(E) ~ 1/sqrt(E) = 1/v.
  * E'r      the resonance energy shifted by the shift factor S_l (zero for l=0).
  * The sin(2 phi) term is interference between resonance scattering and
             potential scattering: it is negative below E'r (a dip in the
             elastic cross section) and positive above.

The weakness of SLBW is that resonances never talk to each other: no level-level
interference, and nothing guarantees the elastic cross section stays positive.
Multi-level formalisms (MLBW, Reich-Moore) fix that.
"""

from __future__ import annotations

import math
from collections.abc import Iterable
from dataclasses import dataclass

from .channel import channel_radius, penetrability_shift_phase, spin_factor, wave_number
from .resonances import LGroup, ResolvedRange, Resonance


@dataclass(frozen=True)
class ResonanceXS:
    """Resonance-range cross sections at one energy (barns)."""

    elastic: float
    capture: float
    fission: float

    @property
    def total(self) -> float:
        return self.elastic + self.capture + self.fission


@dataclass(frozen=True)
class SingleResonance:
    """The derived quantities of one resonance at one energy (for inspection)."""

    k: float  # wave number, 1/(10^-12 cm)
    gn: float  # neutron width at this energy (eV)
    width: float  # total width Gamma (eV)
    shifted_energy: float  # E'r (eV)
    g: float  # spin factor
    phi: float  # hard-sphere phase shift
    capture: float
    fission: float
    elastic_resonance: float  # resonance + interference part of elastic (b)


def _radii(rr: ResolvedRange, group: LGroup) -> tuple[float, float]:
    """(channel radius a, scattering radius AP) for one l group."""
    ap = (
        group.apl if group.apl > 0 else rr.ap
    )  # Reich-Moore may give an l-dependent radius
    return channel_radius(group.awri, ap, rr.naps), ap


def single_resonance(
    rr: ResolvedRange, group: LGroup, res: Resonance, energy: float
) -> SingleResonance:
    """SLBW contribution of one resonance at one energy, with its ingredients."""
    a, ap = _radii(rr, group)
    k = wave_number(energy, group.awri)
    p, s, phi = penetrability_shift_phase(group.l, k * a, k * ap)

    # penetrability and shift at the resonance energy itself
    kr = wave_number(res.energy, group.awri)
    pr, sr, _ = penetrability_shift_phase(group.l, kr * a, kr * ap)

    gn = res.gn * p / pr
    shifted = res.energy + (sr - s) * res.gn / (2.0 * pr)
    gf = abs(res.gfa) + abs(res.gfb)
    width = gn + res.gg + gf
    g = spin_factor(res.j, rr.spi)

    de = energy - shifted
    denom = de * de + 0.25 * width * width
    scale = math.pi / (k * k) * g / denom
    sin_phi = math.sin(phi)
    elastic = scale * (
        gn * gn - 2.0 * gn * width * sin_phi**2 + 2.0 * de * gn * math.sin(2 * phi)
    )
    return SingleResonance(
        k=k,
        gn=gn,
        width=width,
        shifted_energy=shifted,
        g=g,
        phi=phi,
        capture=scale * gn * res.gg,
        fission=scale * gn * gf,
        elastic_resonance=elastic,
    )


def potential_scattering(rr: ResolvedRange, group: LGroup, energy: float) -> float:
    """Hard-sphere scattering of one l group: (4 pi / k^2)(2l + 1) sin^2(phi_l)."""
    a, ap = _radii(rr, group)
    k = wave_number(energy, group.awri)
    _, _, phi = penetrability_shift_phase(group.l, k * a, k * ap)
    return 4.0 * math.pi / (k * k) * (2 * group.l + 1) * math.sin(phi) ** 2


def slbw(
    rr: ResolvedRange, energy: float, only: Iterable[Resonance] | None = None
) -> ResonanceXS:
    """SLBW cross sections at one energy.

    `only` restricts the sum to selected resonances (e.g. one isolated level);
    potential scattering is always included for every l group.
    """
    keep = None if only is None else {id(r) for r in only}
    elastic = capture = fission = 0.0
    for group in rr.groups:
        elastic += potential_scattering(rr, group, energy)
        for res in group.resonances:
            if keep is not None and id(res) not in keep:
                continue
            c = single_resonance(rr, group, res, energy)
            elastic += c.elastic_resonance
            capture += c.capture
            fission += c.fission
    return ResonanceXS(elastic=elastic, capture=capture, fission=fission)


def main(argv: list[str] | None = None) -> None:
    """Scan SLBW cross sections across one resonance.

    uv run python -m nodal.endf.slbw U238 6.674
    """
    import argparse

    from ..data import EndfLibrary
    from .evaluation import Evaluation

    ap = argparse.ArgumentParser(
        prog="python -m nodal.endf.slbw", description=main.__doc__
    )
    ap.add_argument("nuclide")
    ap.add_argument("energy", type=float, help="approximate resonance energy (eV)")
    ap.add_argument(
        "--widths", type=float, default=6.0, help="half-range in total widths"
    )
    args = ap.parse_args(argv)

    rr = (
        Evaluation.from_file(EndfLibrary().neutron(args.nuclide)).resonances[0].resolved
    )
    if rr is None:
        raise SystemExit(f"{args.nuclide} has no resolved resonances")
    group, res = min(
        ((g, r) for g in rr.groups for r in g.resonances),
        key=lambda gr: abs(gr[1].energy - args.energy),
    )
    c = single_resonance(rr, group, res, res.energy)
    print(
        f"{args.nuclide} resonance at {res.energy:g} eV: l={group.l}  J={res.j}  g={c.g:.4g}"
    )
    print(
        f"  Gn={res.gn:.4g} eV  Gg={res.gg:.4g} eV  Gf={abs(res.gfa) + abs(res.gfb):.3g} eV"
        f"  G={c.width:.4g} eV"
    )
    print(
        f"  k={c.k:.5g} /1e-12cm   pi/k^2={math.pi / c.k**2:.5g} b   phi={c.phi:.4g}\n"
    )
    print(
        f"  {'E (eV)':>10} {'capture':>11} {'elastic':>11} {'fission':>11}"
        f" {'all-level el.':>14}   (barns; first three: this resonance alone)"
    )
    for i in range(-10, 11):
        e = res.energy + args.widths * c.width * i / 10
        one, every = slbw(rr, e, only=[res]), slbw(rr, e)
        print(
            f"  {e:10.5f} {one.capture:11.5g} {one.elastic:11.5g} {one.fission:11.3g}"
            f" {every.elastic:14.5g}"
        )


if __name__ == "__main__":
    main()
