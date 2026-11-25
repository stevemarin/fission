"""Zero-power physics tests of Watts Bar 1 cycle 1 (VERA problem 5).

    uv run python -m fission.feedback.zppt wbn1-rods [CASE ...]

The 32 cases of VERA's Table P5-2 (boron, temperature, bank positions) are
the start-up tests: the initial criticality, the all-rods-out critical,
the rod-swap criticals and the bank worths at the dilution endpoint
(1170 ppm), the isothermal temperature coefficient (560 / 570 K) and the
bank D integral worth curve (1230 ppm). Table P5-6 gives the CE KENO-VI
eigenvalue of each case. Every case runs here at its own boron and
temperature (``core.isothermal``: fuel and water at one temperature, zero
power), and the eigenvalue is compared with KENO's. Bank worths, the
differential boron worth and the ITC follow from differences of cases,
the way VERA derives its own (Table P5-7).
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass

from fission.model import Model, load_model

from .core import isothermal


@dataclass(frozen=True)
class Case:
    boron: float  # ppm
    temperature: float  # K, fuel and water
    banks: dict[str, float]  # bank -> steps withdrawn (others out)
    keno: float  # CE KENO-VI k-effective, Table P5-6 (+-1.0-1.5e-5)


# VERA CASL-U-2012-0131-004 Rev. 4, Tables P5-2 and P5-6 (cases 1-32)
CASES: dict[int, Case] = {
    1: Case(1285, 565, {"D": 167}, 0.999899),
    2: Case(1291, 565, {}, 1.000321),
    3: Case(1170, 565, {"A": 0, "D": 97}, 0.998797),
    4: Case(1170, 565, {"B": 0, "D": 113}, 0.999358),
    5: Case(1170, 565, {"C": 0, "D": 119}, 0.999039),
    6: Case(1170, 565, {"D": 18}, 0.999084),
    7: Case(1170, 565, {"SA": 0, "D": 69}, 0.999022),
    8: Case(1170, 565, {"SB": 0, "D": 134}, 0.999324),
    9: Case(1170, 565, {"SC": 0, "D": 71}, 0.998983),
    10: Case(1170, 565, {"SD": 0, "D": 71}, 0.998976),
    11: Case(1170, 565, {}, 1.012841),
    12: Case(1170, 565, {"A": 0}, 1.003716),
    13: Case(1170, 565, {"B": 0}, 1.003941),
    14: Case(1170, 565, {"C": 0}, 1.002843),
    15: Case(1170, 565, {"D": 0}, 0.998815),
    16: Case(1170, 565, {"SA": 0}, 1.008281),
    17: Case(1170, 565, {"SB": 0}, 1.002018),
    18: Case(1170, 565, {"SC": 0}, 1.007749),
    19: Case(1170, 565, {"SD": 0}, 1.007745),
    20: Case(1291, 560, {}, 1.000608),
    21: Case(1291, 570, {}, 1.000034),
    **{
        22 + i: Case(1230, 565, {} if steps is None else {"D": steps}, keno)
        for i, (steps, keno) in enumerate(
            zip(
                (0, 23, 46, 69, 92, 115, 138, 161, 184, 207, None),
                (
                    0.992755,
                    0.993162,
                    0.994555,
                    0.997369,
                    1.000279,
                    1.002542,
                    1.004163,
                    1.005300,
                    1.006073,
                    1.006468,
                    1.006584,
                ),
                strict=True,
            )
        )
    },
}
ROD_WORTH_CASES = {
    "A": 12,
    "B": 13,
    "C": 14,
    "D": 15,
    "SA": 16,
    "SB": 17,
    "SC": 18,
    "SD": 19,
}


def pcm(k_from: float, k_to: float) -> float:
    """Reactivity change rho(k_to) - rho(k_from) in pcm."""
    return 1e5 * (1 / k_from - 1 / k_to)


def run(model: Model, cases: dict[int, Case] | None = None, refine: int = 1,
        axial: float | None = None) -> dict[int, float]:  # fmt: skip
    """k-effective of each case."""
    out = {}
    for n, c in (CASES if cases is None else cases).items():
        out[n] = isothermal(model, [float(c.temperature)], refine, axial,
                            c.banks or None, float(c.boron))[0]  # fmt: skip
    return out


def derived(k: dict[int, float]) -> dict[str, float]:
    """Bank worths at 1170 ppm (each bank alone against case 11, pcm), the
    differential boron worth (cases 2 and 11, pcm/ppm) and the ITC (cases
    20 and 21, pcm/F), from the eigenvalues of the cases."""
    out = {b: pcm(k[n], k[11]) for b, n in ROD_WORTH_CASES.items()}
    out["DBW"] = pcm(k[11], k[2]) / (CASES[2].boron - CASES[11].boron)
    dt = CASES[21].temperature - CASES[20].temperature
    out["ITC"] = pcm(k[20], k[21]) / (1.8 * dt)
    return out


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(prog="python -m fission.feedback.zppt", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)  # fmt: skip
    ap.add_argument("model", nargs="?", default="wbn1-rods")
    ap.add_argument("cases", nargs="*", type=int, help="case numbers (default all)")
    ap.add_argument("--refine", type=int, default=1)
    ap.add_argument("--axial", type=float, default=None, help="largest axial node (cm)")
    args = ap.parse_args(argv)
    model = load_model(args.model)
    cases = {n: CASES[n] for n in args.cases} if args.cases else CASES
    k: dict[int, float] = {}
    for n, c in cases.items():
        k |= run(model, {n: c}, args.refine, args.axial)
        banks = " ".join(f"{b}={s:g}" for b, s in c.banks.items()) or "ARO"
        print(f"case {n:2d}: {c.boron:g} ppm, {c.temperature:g} K, {banks:9s} "
              f"k = {k[n]:.6f}, KENO {c.keno:.6f}, {pcm(c.keno, k[n]):+4.0f} pcm",
              flush=True)  # fmt: skip
    if set(CASES) <= set(k):
        ours, keno = derived(k), derived({n: c.keno for n, c in CASES.items()})
        for key in ours:
            unit = {"DBW": "pcm/ppm", "ITC": "pcm/F"}.get(key, "pcm")
            print(f"  {key:3s} {ours[key]:9.2f} {unit:7s} KENO {keno[key]:9.2f} "
                  f"({100 * (ours[key] / keno[key] - 1):+5.1f}%)")  # fmt: skip


if __name__ == "__main__":
    main()
