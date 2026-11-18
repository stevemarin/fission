"""Legendre polynomials and the ENDF Legendre representation of distributions.

ENDF writes an angular distribution, normalized to 1 over mu in [-1, 1], as

    f(mu) = sum_l (2l + 1)/2 * a_l * P_l(mu),     a_0 = 1,

so the coefficients are the moments a_l = Int f(mu) P_l(mu) dmu. In
particular a_1 is the average cosine. MF4 stores a_1..a_NL (a_0 implied);
MF6 LAW=1 stores unnormalized f_0..f_NA, with f_0 the emission probability.
"""

from __future__ import annotations

from collections.abc import Sequence

# Gauss-Legendre rule, exact for polynomials of degree <= 2*8 - 1 = 15
_NODES = (
    -0.9602898564975362, -0.7966664774136267, -0.525532409916329, -0.18343464249564978,
    0.18343464249564978, 0.525532409916329, 0.7966664774136267, 0.9602898564975362,
)  # fmt: skip
_WEIGHTS = (
    0.10122853629037706, 0.22238103445337443, 0.3137066458778869, 0.36268378337836166,
    0.36268378337836166, 0.3137066458778869, 0.22238103445337443, 0.10122853629037706,
)  # fmt: skip


def legendre(lmax: int, mu: float) -> list[float]:
    """[P_0(mu), ..., P_lmax(mu)] by Bonnet's recurrence
    (l+1) P_{l+1} = (2l+1) mu P_l - l P_{l-1}."""
    p = [1.0, mu][: lmax + 1]
    for l in range(1, lmax):
        p.append(((2 * l + 1) * mu * p[l] - l * p[l - 1]) / (l + 1))
    return p


def legendre_series(coefficients: Sequence[float], mu: float) -> float:
    """sum_l (2l+1)/2 c_l P_l(mu) for c = [c_0, c_1, ...]."""
    p = legendre(len(coefficients) - 1, mu)
    return sum(
        (2 * l + 1) / 2 * c * pl
        for l, (c, pl) in enumerate(zip(coefficients, p, strict=True))
    )


def moments_of_linear(
    mu: Sequence[float], f: Sequence[float], lmax: int
) -> list[float]:
    """a_l = Int f P_l dmu, l = 0..lmax, for f linear between the points (mu_i, f_i).

    On each segment the integrand is a polynomial of degree l + 1, so an
    8-point Gauss rule is exact for l <= 14.
    """
    if lmax > 14:
        raise ValueError("lmax > 14 needs a larger Gauss rule")
    a = [0.0] * (lmax + 1)
    for i in range(len(mu) - 1):
        m1, m2, f1, f2 = mu[i], mu[i + 1], f[i], f[i + 1]
        h = m2 - m1
        if h <= 0.0:
            continue
        for t, w in zip(_NODES, _WEIGHTS, strict=True):
            x = m1 + 0.5 * h * (t + 1.0)
            fx = f1 + (f2 - f1) * (x - m1) / h
            for l, pl in enumerate(legendre(lmax, x)):
                a[l] += 0.5 * h * w * fx * pl
    return a
