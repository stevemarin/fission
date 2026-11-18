"""Doppler broadening of linear-linear cross-section tables (free-gas target).

Clean Python reference; slow (a full U-238 table needs the Rust port).

Theory
------
The 0 K table gives sigma_0 for a target nucleus at rest. Real nuclei move:
with a Maxwellian (free-gas) velocity distribution at temperature T, the
reaction rate a neutron of lab speed v sees is the average of
|v - V| sigma_0(|v - V|) over target velocities V. Doing the angular part of
that average analytically leaves a one-dimensional integral (the "exact
kernel" of Solbrig, and Cullen & Weisbin's SIGMA1). In reduced speeds

    y^2 = alpha E,   x^2 = alpha E',   alpha = AWR / (k_B T),

(y is the neutron speed in units of the most probable target speed) it reads

    sigma(y, T) = 1 / (y^2 sqrt(pi)) * Int_0^inf x^2 sigma_0(x) K(x, y) dx,
    K(x, y)     = exp(-(x - y)^2) - exp(-(x + y)^2).

Only |x - y| < ~6 contributes (exp(-36) ~ 2e-16), i.e. relative energies
within a few Doppler widths of E.

Because the integral is linear in sigma_0, and sigma_0 is linear between
table points, the broadened value at one energy is a fixed weighted sum of
the table values: sigma(y) = sum_i w_i(y) sigma_i. The same weights broaden
every reaction (`broadening_weights`).

On each segment the integrand is a polynomial times a Gaussian. Cullen and
Weisbin integrate it in closed form with error functions
(`cullen_weisbin`, kept as an independent check). The main routine uses
Gauss-Legendre quadrature on panels at most 0.5 wide in x instead (2, 4 or
8 points depending on the width): the integrand is so smooth there that this
reaches machine precision, and it
avoids the cancellation of the closed form on very narrow segments (the
difference of two nearly equal error functions) as well as any need for
erfc (the Rust standard library has none).

Outside the table: below the first point sigma_0 is extended as 1/v, above
the last point as a constant.
"""

from __future__ import annotations

import argparse
import bisect
import math
from collections import defaultdict
from collections.abc import Callable, Sequence

from fission.endf.pointwise import PointwiseTable
from fission.endf.reconstruct import adaptive_grid

BOLTZMANN_EV_PER_K = 8.617333262e-5  # k_B (eV/K), exact in SI 2019 / CODATA 2018
Z_CUTOFF = 6.0  # kernel neglected beyond |x - y| > 6: exp(-36) ~ 2e-16
PANEL_WIDTH = 0.5  # widest quadrature panel, in reduced speed x
# Gauss-Legendre rules on [-1, 1] (= numpy.polynomial.legendre.leggauss(n)),
# written out so that Python and Rust use bit-identical values
GAUSS_RULES: dict[int, tuple[tuple[float, ...], tuple[float, ...]]] = {
    2: ((-0.5773502691896257, 0.5773502691896257), (1.0, 1.0)),
    4: (
        (-0.8611363115940526, -0.33998104358485626, 0.33998104358485626, 0.8611363115940526),
        (0.34785484513745357, 0.6521451548625464, 0.6521451548625464, 0.34785484513745357),
    ),
    8: (
        (-0.9602898564975362, -0.7966664774136267, -0.525532409916329, -0.18343464249564978,
         0.18343464249564978, 0.525532409916329, 0.7966664774136267, 0.9602898564975362),
        (0.10122853629037706, 0.22238103445337443, 0.3137066458778869, 0.36268378337836166,
         0.36268378337836166, 0.3137066458778869, 0.22238103445337443, 0.10122853629037706),
    ),
}  # fmt: skip


def gauss_rule(h: float) -> tuple[tuple[float, ...], tuple[float, ...]]:
    """Fewest points that integrate a panel of width h (in x) to ~1e-14.

    The n-point error scales as h^(2n+1): most segments of a resonance table
    are narrow (h ~ 0.001-0.01), and 2 points are already exact there.
    """
    return GAUSS_RULES[2 if h <= 0.01 else 4 if h <= 0.1 else 8]


def alpha(awr: float, temperature: float) -> float:
    """alpha = AWR / (k_B T) in 1/eV, so that y^2 = alpha * E."""
    return awr / (BOLTZMANN_EV_PER_K * temperature)


def doppler_width(energy: float, awr: float, temperature: float) -> float:
    """Gamma_D = 2 sqrt(E k_B T / AWR) (eV): the energy scale of the broadening."""
    return 2.0 * math.sqrt(energy * BOLTZMANN_EV_PER_K * temperature / awr)


def kernel(x: float, y: float) -> float:
    """K(x, y) = exp(-(x-y)^2) - exp(-(x+y)^2), without cancellation at small x*y."""
    return math.exp(-((x - y) ** 2)) * -math.expm1(-4.0 * x * y)


def _integrate(f: Callable[[float], float], a: float, b: float, y: float) -> float:
    """Int_a^b f(x) K(x, y) dx by Gauss-Legendre on panels <= PANEL_WIDTH."""
    n = max(1, math.ceil((b - a) / PANEL_WIDTH))
    h = (b - a) / n
    total = 0.0
    for p in range(n):
        mid = a + (p + 0.5) * h
        for t, w in zip(*gauss_rule(h), strict=True):
            x = mid + 0.5 * h * t
            total += w * f(x) * kernel(x, y)
    return 0.5 * h * total


def _segment_weights(
    xl: float, xr: float, a: float, b: float, y: float
) -> tuple[float, float]:
    """Contributions of the end values of the segment [xl, xr], integrated over [a, b].

    Between table points sigma_0 is linear in E, i.e. in x^2:
        sigma_0(x) = sigma_l h_l(x) + sigma_r h_r(x),
        h_l = (xr^2 - x^2) / (xr^2 - xl^2),   h_r = (x^2 - xl^2) / (xr^2 - xl^2).
    Returns Int x^2 h_l K dx and Int x^2 h_r K dx (factored, no cancellation).
    """
    d = (xr - xl) * (xr + xl)
    n = max(1, math.ceil((b - a) / PANEL_WIDTH))
    h = (b - a) / n
    left = right = 0.0
    for p in range(n):
        mid = a + (p + 0.5) * h
        for t, w in zip(*gauss_rule(h), strict=True):
            x = mid + 0.5 * h * t
            f = w * x * x * kernel(x, y) / d
            left += f * (xr - x) * (xr + x)
            right += f * (x - xl) * (x + xl)
    return 0.5 * h * left, 0.5 * h * right


def broadening_weights(x: Sequence[float], y: float) -> dict[int, float]:
    """Weights w_i with sigma(y, T) = sum_i w_i sigma_0[i].

    x are the table energies in reduced speed (x_i = sqrt(alpha E_i)),
    increasing; y is the reduced speed of the output energy.
    """
    lo, hi = max(0.0, y - Z_CUTOFF), y + Z_CUTOFF
    w: dict[int, float] = defaultdict(float)

    # below the table: sigma_0 = sigma_0[0] * x_0 / x  (1/v), so x^2 sigma_0 = sigma_0[0] x_0 x
    if lo < x[0]:
        x0 = x[0]
        w[0] += _integrate(lambda t: x0 * t, lo, min(x0, hi), y)

    # inside: sigma_0 linear in E, i.e. in x^2, between neighbouring points
    i = max(bisect.bisect_right(x, lo) - 1, 0)
    while i < len(x) - 1 and x[i] < hi:
        a, b = max(x[i], lo), min(x[i + 1], hi)
        if b > a:
            left, right = _segment_weights(x[i], x[i + 1], a, b, y)
            w[i] += left
            w[i + 1] += right
        i += 1

    # above the table: constant
    if hi > x[-1]:
        w[len(x) - 1] += _integrate(lambda t: t * t, max(x[-1], lo), hi, y)

    norm = 1.0 / (y * y * math.sqrt(math.pi))
    return {k: v * norm for k, v in w.items()}


def broadened_at(
    table: PointwiseTable, awr: float, temperature: float, energies: Sequence[float]
) -> dict[int, list[float]]:
    """Broadened cross sections of every reaction in `table` at `energies`."""
    a = alpha(awr, temperature)
    x = [math.sqrt(a * e) for e in table.energies]
    out: dict[int, list[float]] = {mt: [] for mt in table.xs}
    for e in energies:
        w = broadening_weights(x, math.sqrt(a * e))
        for mt, values in table.xs.items():
            out[mt].append(math.fsum(wi * values[i] for i, wi in w.items()))
    return out


def broaden(
    table: PointwiseTable, awr: float, temperature: float, emax: float | None = None
) -> PointwiseTable:
    """The table at `temperature` (K), on the same energy grid.

    Points above `emax` keep their 0 K values (broadening is negligible at
    high energy, where cross sections are smooth on the Doppler scale).
    Broadening is linear, so MT1 stays the sum of the partials (to rounding).
    """
    if temperature == 0.0:
        return PointwiseTable(
            list(table.energies), {k: list(v) for k, v in table.xs.items()}
        )
    n = (
        len(table.energies)
        if emax is None
        else bisect.bisect_right(table.energies, emax)
    )
    hot = broadened_at(table, awr, temperature, table.energies[:n])
    xs = {mt: hot[mt] + list(values[n:]) for mt, values in table.xs.items()}
    return PointwiseTable(list(table.energies), xs)


# ------------------------------------------------ tables at temperature
def thin(
    energies: Sequence[float],
    values: Sequence[Sequence[float]],
    tol: float,
    atol: float,
) -> list[int]:
    """Indices of a subset of points that still represents the table.

    Greedy: from the last kept point i, extend the chord i -> j as far as
    every skipped point k (i < k < j) stays within tol |sigma_k| + atol of
    it, for every reaction; then keep j. Broadened tables are smoother than
    the 0 K ones they are seeded with, so many points become redundant.
    """
    n = len(energies)
    keep = [0]
    i = 0
    while i < n - 1:
        j = i + 1
        while j + 1 < n and _chord_holds(energies, values, i, j + 1, tol, atol):
            j += 1
        keep.append(j)
        i = j
    return keep


def _chord_holds(
    e: Sequence[float],
    v: Sequence[Sequence[float]],
    i: int,
    j: int,
    tol: float,
    atol: float,
) -> bool:
    """Do all points strictly between i and j lie within tolerance of the chord i -> j?"""
    for k in range(i + 1, j):
        t = (e[k] - e[i]) / (e[j] - e[i])
        for vi, vk, vj in zip(v[i], v[k], v[j], strict=True):
            if abs(vk - (vi + t * (vj - vi))) > tol * abs(vk) + atol:
                return False
    return True


SEED_WIDTHS = (1.0, 2.0, 3.0)  # extra seeds at E_i +- k Gamma_D in sparse regions


def doppler_seeds(
    energies: Sequence[float], awr: float, temperature: float, n: int
) -> list[float]:
    """Starting grid for a broadened table: the first n 0 K points, plus
    points at 1, 2 and 3 Doppler widths on either side of each of them
    wherever the 0 K grid is sparse.

    Broadening rounds every corner of the 0 K curve over a few Gamma_D (the
    kernel is ~exp(-(dE / Gamma_D)^2)). Between two corners far apart the
    broadened curve is straight, so a midpoint test can be fooled: a concave
    corner at one end and a convex one at the other shift the chord by
    opposite amounts that cancel at the midpoint. Seeds inside the rounded
    regions make the bisection look there.
    """
    seeds = set(energies[:n])
    for i in range(n):
        e = energies[i]
        gd = doppler_width(e, awr, temperature)
        for k in SEED_WIDTHS:
            if i + 1 < n and e + k * gd < energies[i + 1] - 0.5 * gd:
                seeds.add(e + k * gd)
            if i > 0 and e - k * gd > energies[i - 1] + 0.5 * gd:
                seeds.add(e - k * gd)
    return sorted(seeds)


def broaden_table(
    table: PointwiseTable,
    awr: float,
    temperature: float,
    *,
    tol: float = 1e-3,
    atol: float = 1e-5,
    emax: float | None = None,
    min_rel_width: float = 1e-9,
    use_rust: bool = True,
) -> PointwiseTable:
    """The table at `temperature` (K) on its own linear-linear grid.

    1. Evaluate the broadened cross sections at the 0 K points (up to emax),
       plus a few Doppler widths around them where they are sparse
       (`doppler_seeds`).
    2. Refine by the midpoint test (as in reconstruction): needed where the
       broadened curve bends more than the 0 K one, e.g. the low-energy
       free-gas rise of scattering cross sections.
    3. Thin: drop points that linear interpolation reproduces within tol.
    Points above emax are copied unbroadened.
    """
    mts = table.mts
    n = (
        len(table.energies)
        if emax is None
        else bisect.bisect_right(table.energies, emax)
    )
    seeds = doppler_seeds(table.energies, awr, temperature, n)
    if use_rust:
        from fission import _rust

        grid, cols = _rust.doppler_broaden_table(
            table.energies,
            [table.xs[mt] for mt in mts],
            awr,
            temperature,
            seeds,
            tol,
            atol,
            min_rel_width,
        )
        energies: list[float] = grid.tolist()
        xs = {mt: c.tolist() for mt, c in zip(mts, cols, strict=True)}
    else:
        a = alpha(awr, temperature)
        x = [math.sqrt(a * e) for e in table.energies]
        columns = [table.xs[mt] for mt in mts]

        def f(e: float) -> tuple[float, ...]:
            w = broadening_weights(x, math.sqrt(a * e))
            return tuple(
                math.fsum(wi * col[i] for i, wi in w.items()) for col in columns
            )

        grid_py, values = adaptive_grid(f, seeds, tol, atol, min_rel_width)
        keep = thin(grid_py, values, tol, atol)
        energies = [grid_py[k] for k in keep]
        xs = {mt: [values[k][m] for k in keep] for m, mt in enumerate(mts)}
    for mt in mts:
        xs[mt] += table.xs[mt][n:]
    return PointwiseTable(energies + list(table.energies[n:]), xs)


# ------------------------------------------------- closed form (reference)
def _gauss_tails(a: float) -> list[float]:
    """F_k(a) = Int_a^inf z^k exp(-z^2) dz for k = 0..4."""
    if a == math.inf:
        return [0.0] * 5
    e = math.exp(-a * a)
    f = [0.5 * math.sqrt(math.pi) * math.erfc(a), 0.5 * e]
    for k in range(2, 5):  # integration by parts
        f.append(0.5 * (k - 1) * f[k - 2] + 0.5 * a ** (k - 1) * e)
    return f


def _poly_gauss(c: Sequence[float], shift: float, lo: float, hi: float) -> float:
    """Int_lo^hi p(x) exp(-(x - shift)^2) dx for p(x) = sum_n c[n] x^n (n <= 4)."""
    # x = z + shift:  x^n = sum_k C(n, k) z^k shift^(n-k)
    d = [0.0] * 5
    for n, cn in enumerate(c):
        for k in range(n + 1):
            d[k] += cn * math.comb(n, k) * shift ** (n - k)
    fa, fb = _gauss_tails(lo - shift), _gauss_tails(hi - shift)
    return sum(dk * (fa[k] - fb[k]) for k, dk in enumerate(d))


def cullen_weisbin(x: Sequence[float], sigma: Sequence[float], y: float) -> float:
    """Broadened sigma at y by the closed-form (SIGMA1) integration of every segment.

    Same extrapolation as `broadening_weights`, but no cutoff. Accurate for
    moderate y and segments that are not extremely narrow; used in tests.
    """
    pieces: list[tuple[tuple[float, ...], float, float]] = [
        ((0.0, sigma[0] * x[0]), 0.0, x[0])  # 1/v below the table
    ]
    for i in range(len(x) - 1):
        b = (sigma[i + 1] - sigma[i]) / (x[i + 1] ** 2 - x[i] ** 2)  # sigma = A + B x^2
        pieces.append(((0.0, 0.0, sigma[i] - b * x[i] ** 2, 0.0, b), x[i], x[i + 1]))
    pieces.append(((0.0, 0.0, sigma[-1]), x[-1], math.inf))  # constant above
    total = math.fsum(
        _poly_gauss(c, y, lo, hi) - _poly_gauss(c, -y, lo, hi) for c, lo, hi in pieces
    )
    return total / (y * y * math.sqrt(math.pi))


# ------------------------------------------------------------------- CLI
def main(argv: list[str] | None = None) -> None:
    from fission.data import EndfLibrary
    from fission.endf import Evaluation
    from fission.endf.pointwise import pointwise

    parser = argparse.ArgumentParser(
        prog="python -m fission.processing.doppler",
        description="Doppler-broaden the cross sections around one resonance.",
    )
    parser.add_argument("nuclide", help="e.g. U238")
    parser.add_argument("energy", type=float, help="resonance energy (eV)")
    parser.add_argument(
        "--temps", type=float, nargs="+", default=[293.6, 600.0, 900.0, 1200.0]
    )
    parser.add_argument(
        "--mt", type=int, default=102, help="reaction to show (default capture)"
    )
    args = parser.parse_args(argv)

    ev = Evaluation.from_file(EndfLibrary().neutron(args.nuclide))
    awr, e0 = ev.header.awr, args.energy
    gd_max = doppler_width(e0, awr, max(args.temps))
    inner = (max(e0 - 4 * gd_max, 1e-5), e0 + 4 * gd_max)  # where values are reported
    margin = 2 * Z_CUTOFF * gd_max  # extra 0 K data needed around it
    cold = pointwise(
        ev, mts=[args.mt], lo=max(inner[0] - margin, 1e-5), hi=inner[1] + margin
    )
    grid = [e for e in cold.energies if inner[0] <= e <= inner[1]]
    print(
        f"{args.nuclide}  AWR = {awr:.4f}  MT{args.mt} around {e0} eV  ({len(grid)} points)"
    )
    print(
        f"{'T (K)':>8} {'Gamma_D (eV)':>13} {'peak (b)':>11} {'at E (eV)':>10} {'Int dE (b eV)':>14}"
    )
    for temp in [0.0, *args.temps]:
        if temp == 0.0:
            values = [cold(args.mt, e) for e in grid]
        else:
            values = broadened_at(cold, awr, temp, grid)[args.mt]
        peak = max(range(len(grid)), key=values.__getitem__)
        area = math.fsum(
            0.5 * (values[i] + values[i + 1]) * (grid[i + 1] - grid[i])
            for i in range(len(grid) - 1)
        )
        gd = doppler_width(e0, awr, temp) if temp else 0.0
        print(
            f"{temp:8.1f} {gd:13.4f} {values[peak]:11.1f} {grid[peak]:10.4f} {area:14.3f}"
        )


if __name__ == "__main__":
    main()
