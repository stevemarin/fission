"""Summary of a nuclide's secondary-neutron data.

uv run python -m fission.endf.distributions U235
"""

from __future__ import annotations

import argparse
import math

from ...data import EndfLibrary
from ..evaluation import Evaluation
from .angular import angular_distribution
from .energy import energy_distribution
from .kinematics import elastic_alpha, mean_lethargy_gain, two_body
from .mf6 import Continuum, TwoBody, energy_angle
from .nu import delayed_groups, multiplicity_distribution, nu_bar


def _mean_lab_cosine(ad, energy: float, awr: float, n: int = 400) -> float:
    h = 2.0 / n
    return h * math.fsum(
        ad(energy, mu) * two_body(energy, mu, awr)[1]
        for mu in (-1 + (k + 0.5) * h for k in range(n))
    )


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(prog="python -m fission.endf.distributions")
    ap.add_argument("nuclide")
    args = ap.parse_args(argv)
    ev = Evaluation.from_file(EndfLibrary().neutron(args.nuclide))
    awr = ev.header.awr
    print(f"{args.nuclide}: AWR = {awr:.4f}")

    if (1, 452) in ev.sections:
        nt, npr, nd = (nu_bar(ev, mt) for mt in (452, 456, 455))
        print("\nneutrons per fission       E (eV)   nu_t     nu_p     nu_d     beta")
        for e in (0.0253, 1e6, 1.4e7):
            print(
                f"{'':25s} {e:8.3g} {nt(e):8.4f} {npr(e):8.4f} {nd(e):8.5f} {nd(e) / nt(e):8.5f}"
            )
        g = delayed_groups(ev)
        print(
            "\ndelayed groups: half-life (s), abundance, mean energy (keV) at 0.0253 eV"
        )
        for t, a, m in zip(
            g.half_lives, g.abundances(0.0253), g.mean_energies(0.0253), strict=True
        ):
            print(f"  {t:8.3f} {a:8.4f} {m / 1e3:8.1f}")
        md = multiplicity_distribution(ev)
        if md:
            p = md.probabilities(0.0253)
            print("\nP(nu) at 0.0253 eV:", " ".join(f"{v:.3f}" for v in p[:8]))
        chi = energy_distribution(ev, 18)
        print(
            f"\nprompt fission spectrum: <E'> = {chi.mean(0.0253) / 1e6:.3f} MeV at thermal, "
            f"{chi.mean(1.4e7) / 1e6:.3f} MeV at 14 MeV"
        )

    print(
        f"\nelastic: alpha = {elastic_alpha(awr):.4f}, xi = {mean_lethargy_gain(awr):.4f}, "
        f"2/(3A) = {2 / (3 * awr):.4f}"
    )
    ad = angular_distribution(ev, 2)
    for e in (1.0, 1e5, 1e6, 1e7):
        print(
            f"  E = {e:8.0e} eV: <mu_cm> = {ad.mean_cosine(e):+.4f}, "
            f"<mu_lab> = {_mean_lab_cosine(ad, e, awr):+.4f}"
        )

    print("\nneutron-emitting reactions with MF6 data:")
    levels = []
    for mf, mt in sorted(ev.sections):
        if mf != 6:
            continue
        products = energy_angle(ev, mt).products
        if not products:
            continue
        d = products[0].distribution
        if isinstance(d, TwoBody):
            levels.append(mt)
        elif isinstance(d, Continuum):
            print(
                f"  MT{mt:<4d} continuum (LANG={d.lang}), yield {products[0].yield_(1.4e7):.0f}"
            )
        else:
            print(
                f"  MT{mt:<4d} LAW={products[0].law} ({len(products)} neutron products)"
            )
    if levels:
        print(f"  MT{levels[0]}-{levels[-1]}: {len(levels)} discrete levels (two-body)")


if __name__ == "__main__":
    main()
