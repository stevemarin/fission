"""Elastic scattering on moving targets: the exact free-gas kernel.

The slowing-down march of chapter 16 scatters on nuclei at rest: a neutron
of energy E' leaves with an energy uniform in [alpha E', E']. In hot fuel the
U-238 nuclei move with a Maxwellian velocity distribution, and near a
resonance this matters, because the cross section depends on the *relative*
speed of neutron and nucleus, and resonances are only a few Doppler widths
wide.

Kernel. Speeds are measured in reduced units x = sqrt(A E / kT) (the
variable of Doppler broadening, ``doppler.alpha``): x_in, x_out for the
neutron before and after, z for the relative speed (relative energy
E_r = z^2 kT / A). Scattering is isotropic in the centre of mass. Averaging
over the target's Maxwellian and over the scattering angle gives [Ouisloumen
& Sanchez 1991; derivation in ch. 19]

    sigma_s(E_in -> E_out) = (A + 1)^2 / (2 A^2 E_in sqrt(pi))
                             Int_0^inf dz z sigma_0(E_r(z)) J(z),

    J(z) = Int_{w_lo}^{w_hi} exp(x_in^2/A - z^2/(A+1) - w^2/a) dw,

with a = A / (A + 1), w the speed of the centre of mass, w_lo = max(|x_in -
a z|, |x_out - a z|) and w_hi = min(x_in, x_out) + a z. sigma_0 is the
*0 K* elastic cross section: the target motion is in the kernel, not in the
cross section. The exponent is minus the squared speed of the target, so it
is never positive, and it is at most -(x_in - z)^2: only relative speeds
within a few units of x_in contribute (``Z_CUTOFF``). Integrated over E_out
the kernel gives the Doppler-broadened sigma_s(E_in) of ``doppler.py``, and
as T -> 0 it becomes the box 1 / ((1 - alpha) E_in) of the target at rest.

Computation. With w_hi ~ 2 x_in the upper limit contributes exp(-4 x^2),
nothing for x > 10. Then J depends on x_out only through which of
|x_in - a z| and |x_out - a z| is larger, and that switches once, at
z* = (x_in + x_out) / (2 a). On each side J is a function of one node only:

    J(z) = D(x_in, z)                                     (x_in farther)
    J(z) = exp((x_in^2 - x_out^2) / A) D(x_out, z)         (x_out farther)

    D(x, z) = exp(x^2/A - z^2/(A+1) - (x - a z)^2/a) sqrt(pi a)/2 erfcx(|x - a z| / sqrt(a)).

The factor exp((x_in^2 - x_out^2)/A) = exp((E_in - E_out)/kT) is detailed
balance. On a uniform z grid the integral over z becomes a sum with the
weights W_j = Int hat_j(z) z sigma_0(z) dz (hat functions; exact for the
0 K table, which is linear between its points), and the two pieces are a
prefix and a suffix sum of W D for one node each. Both are kept as
logarithms (the detailed-balance factor reaches e^300 at 1 keV).

The kernel lives on a mesh uniform in x (spacing ``dx``, a fraction of the
Doppler width, which is 1 in x): incoming nodes from the bottom of the march
to ``e_top`` (and TAIL below e_low, for neutrons that scatter up into the
range), outgoing nodes wider by the tails. ``MovingKernel`` keeps
P(E_m -> E_n) = sigma_s(E_m -> E_n) / sigma_s(E_m) (normalized on the mesh,
so that the march conserves neutrons) by outgoing node, for the march.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from fission.processing.doppler import GAUSS_RULES
from fission.processing.doppler import alpha as reduced_scale

SQRT_PI = math.sqrt(math.pi)
Z_CUTOFF = 6.0  # |x_in - z| beyond this: exp(-36) ~ 2e-16
TAIL = (
    2.0 * Z_CUTOFF + 1.0
)  # outgoing x beyond [sqrt(alpha) x_in, x_in] by at most this
ERFCX_SWITCH = 2.0  # series below, continued fraction above (both ~1e-13)
SERIES_TERMS = 60
FRACTION_TERMS = 80


def erfcx(y: np.ndarray) -> np.ndarray:
    """Scaled complementary error function exp(y^2) erfc(y), y >= 0.

    Below 2: exp(y^2) - 2/sqrt(pi) sum_n y (2 y^2)^n / (2n + 1)!! (the series
    of erf with positive terms). Above: the continued fraction erfcx(y) =
    1 / (sqrt(pi) (y + (1/2) / (y + 1 / (y + (3/2) / (y + ...))))),
    evaluated backwards from a fixed depth. Fixed term counts, so that the
    Rust port does the same operations. Relative error < 1e-13.
    """
    y = np.asarray(y, dtype=float)
    out = np.empty_like(y)
    small = y < ERFCX_SWITCH
    ys = y[small]
    t2 = 2.0 * ys * ys
    term = ys.copy()
    total = ys.copy()
    for n in range(1, SERIES_TERMS):
        term = term * t2 / (2 * n + 1)
        total = total + term
    out[small] = np.exp(ys * ys) - 2.0 / SQRT_PI * total
    yl = y[~small]
    t = yl.copy()
    for k in range(FRACTION_TERMS, 0, -1):
        t = yl + 0.5 * k / t
    out[~small] = 1.0 / (SQRT_PI * t)
    return out


@dataclass(frozen=True)
class MovingKernel:
    """Normalized free-gas elastic kernel of one nuclide on a uniform x mesh.

    Node n sits at x_n = base + n dx, energy E_n = x_n^2 / scale. Incoming
    nodes go up to `top` (the node at the top of the moving-target range);
    outgoing nodes reach TAIL above. Entries are stored by outgoing node:
    starts[n]:starts[n+1] index rows (incoming node m) and values
    P(E_m -> E_n) in 1/eV. norm[m] = Int sigma_s(E_m -> E) dE on the mesh
    (barns; nan for nodes without a row): the broadened sigma_s(E_m).
    """

    nuclide: str
    awr: float
    temperature: float
    scale: float  # A / kT (1/eV): x = sqrt(scale E)
    base: float
    dx: float
    top: int
    nodes: int
    starts: np.ndarray
    rows: np.ndarray
    values: np.ndarray
    norm: np.ndarray

    @property
    def x(self) -> np.ndarray:
        return self.base + self.dx * np.arange(self.nodes)

    @property
    def energies(self) -> np.ndarray:
        return self.x**2 / self.scale

    @property
    def e_top(self) -> float:
        return float(self.energies[self.top])

    def column(self, n: int) -> tuple[np.ndarray, np.ndarray]:
        s = slice(self.starts[n], self.starts[n + 1])
        return self.rows[s], self.values[s]

    def dense(self) -> np.ndarray:
        """P as a (nodes, nodes) array [m, n] (for tests)."""
        p = np.zeros((self.nodes, self.nodes))
        for n in range(self.nodes):
            rows, values = self.column(n)
            p[rows, n] = values
        return p


def hat_weights(
    z0: float,
    dz: float,
    count: int,
    energies: np.ndarray,
    sigma: np.ndarray,
    scale: float,
) -> np.ndarray:
    """W_j = Int hat_j(z) z sigma_0(z^2 / scale) dz on the grid z_j = z0 + j dz.

    sigma_0 is linear in E between the table points, so on every piece
    between grid nodes and table points the integrand is a polynomial of
    degree 4 in z: the 4-point Gauss rule is exact.
    """
    z_hi = z0 + (count - 1) * dz
    z_table = np.sqrt(scale * np.asarray(energies))
    inside = z_table[(z_table > z0) & (z_table < z_hi)]
    points = np.union1d(z0 + dz * np.arange(count), inside)
    lo, hi = points[:-1], points[1:]
    cell = np.minimum(np.floor((0.5 * (lo + hi) - z0) / dz).astype(np.int64), count - 2)
    nodes, weights = GAUSS_RULES[4]
    w = np.zeros(count)
    for t, wt in zip(nodes, weights, strict=True):
        z = 0.5 * (lo + hi) + 0.5 * (hi - lo) * t
        f = 0.5 * (hi - lo) * wt * z * np.interp(z * z / scale, energies, sigma)
        right = (z - (z0 + cell * dz)) / dz  # hat of the upper node
        w += np.bincount(cell, f * (1.0 - right), minlength=count)
        w += np.bincount(cell + 1, f * right, minlength=count)
    return w


def free_gas_kernel(
    nuclide: str,
    energies: np.ndarray,
    sigma0: np.ndarray,
    awr: float,
    temperature: float,
    e_low: float,
    e_top: float,
    dx: float = 0.125,
    dz: float = 1.0 / 64.0,
    use_rust: bool = False,
) -> MovingKernel:
    """The kernel of one nuclide for incoming energies in [e_low, e_top].

    energies, sigma0: the 0 K elastic table (MT2). dx: mesh spacing in x
    (Doppler widths); dz: spacing of the relative-speed grid.
    """
    if use_rust:
        from fission import _rust

        return MovingKernel(nuclide, awr, temperature, *_rust.free_gas_kernel(
            np.ascontiguousarray(energies, dtype=float),
            np.ascontiguousarray(sigma0, dtype=float),
            awr, temperature, e_low, e_top, dx, dz,
        ))  # fmt: skip
    scale, base, top, nodes, z0, count = _mesh(awr, temperature, e_low, e_top, dx, dz)
    w = hat_weights(z0, dz, count, energies, sigma0, scale)
    starts, rows, values, norm = _assemble(awr, scale, base, dx, top, nodes, z0, dz,
                                           math.sqrt(scale * e_low), w)  # fmt: skip
    return MovingKernel(nuclide, awr, temperature, scale, base, dx, top, nodes,
                        starts, rows, values, norm)  # fmt: skip


def _mesh(awr, temperature, e_low, e_top, dx, dz):
    """Mesh in x (node `top` at e_top) and the z grid that covers it."""
    scale = reduced_scale(awr, temperature)
    sqrt_alpha = (awr - 1.0) / (awr + 1.0)
    x_top, x_low = math.sqrt(scale * e_top), math.sqrt(scale * e_low)
    # incoming rows reach TAIL below e_low (upscattering into the range)
    below = math.ceil((x_top - (sqrt_alpha * (x_low - TAIL) - TAIL)) / dx)
    above = math.ceil(TAIL / dx)
    base = x_top - below * dx
    if base < 10.0:
        raise ValueError(
            f"x = {base:.1f} at the bottom of the mesh: the upper limit of the "
            "centre-of-mass integral is not negligible (raise e_low)"
        )
    nodes = below + above + 1
    x_max = base + (nodes - 1) * dx
    z0 = base - TAIL - Z_CUTOFF
    z1 = x_max / sqrt_alpha + TAIL + Z_CUTOFF
    count = math.ceil((z1 - z0) / dz) + 1
    return scale, base, below, nodes, z0, count


def _log_add(a: np.ndarray) -> np.ndarray:
    return np.logaddexp.accumulate(a)


def _assemble(awr, scale, base, dx, top, nodes, z0, dz, x_low, w):
    """Rows of every incoming node, then sorted by outgoing node."""
    a = awr / (awr + 1.0)
    sqrt_alpha = (awr - 1.0) / (awr + 1.0)
    front = math.sqrt(math.pi * a) / 2.0
    log_w = np.full(len(w), -np.inf)
    log_w[w > 0] = np.log(w[w > 0])
    span: dict[
        int, tuple[int, np.ndarray, np.ndarray]
    ] = {}  # n -> (i_lo, prefix, suffix)

    def node(s: int) -> tuple[int, np.ndarray, np.ndarray]:
        if s not in span:
            xs = base + s * dx
            i_lo = max(math.floor((xs - TAIL - Z_CUTOFF - z0) / dz), 0)
            i_hi = min(
                math.ceil((xs / sqrt_alpha + TAIL + Z_CUTOFF - z0) / dz), len(w) - 1
            )
            z = z0 + dz * np.arange(i_lo, i_hi + 1)
            gap = np.abs(xs - a * z)
            log_d = (
                xs * xs / awr - z * z / (awr + 1.0) - gap * gap / a
                + np.log(front * erfcx(gap / math.sqrt(a)))
                + log_w[i_lo : i_hi + 1]
            )  # fmt: skip
            span[s] = (i_lo, _log_add(log_d), _log_add(log_d[::-1])[::-1])
        return span[s]

    def prefix(s: int, j: int) -> float:  # log sum_{i <= j} W D_s(z_i)
        i_lo, pre, _ = node(s)
        if j < i_lo:
            return -math.inf
        return float(pre[min(j - i_lo, len(pre) - 1)])

    def suffix(s: int, j: int) -> float:  # log sum_{i >= j} W D_s(z_i)
        i_lo, _, suf = node(s)
        if j - i_lo >= len(suf):
            return -math.inf
        return float(suf[max(j - i_lo, 0)])

    first = max(math.floor((x_low - TAIL - base) / dx), 0)
    all_m, all_n, all_k = [], [], []
    norm = np.full(nodes, np.nan)
    for m in range(first, top + 1):
        xm = base + m * dx
        lo = max(math.ceil((sqrt_alpha * xm - TAIL - base) / dx), 0)
        hi = min(math.floor((xm + TAIL - base) / dx), nodes - 1)
        k = np.empty(hi - lo + 1)
        for idx, n in enumerate(range(lo, hi + 1)):
            xn = base + n * dx
            if n == m:
                k[idx] = math.exp(prefix(m, len(w)))
                continue
            j = math.floor(((xm + xn) / (2.0 * a) - z0) / dz)  # last z_i <= z*
            balance = (xm * xm - xn * xn) / awr
            if n < m:  # x_in farther below z*, x_out above
                k[idx] = math.exp(prefix(m, j)) + math.exp(balance + suffix(n, j + 1))
            else:
                k[idx] = math.exp(suffix(m, j + 1)) + math.exp(balance + prefix(n, j))
        k *= (awr + 1.0) ** 2 / (2.0 * awr * awr * xm * xm / scale * SQRT_PI)
        xn = base + dx * np.arange(lo, hi + 1)
        total = float(
            np.sum(k * 2.0 * xn * dx / scale)
        )  # Int K dE: dE = 2 x dx / scale
        norm[m] = total
        all_m.append(np.full(hi - lo + 1, m))
        all_n.append(np.arange(lo, hi + 1))
        all_k.append(k / total)
        for s in [s for s in span if s < lo - 1]:
            del span[s]
    m_all = np.concatenate(all_m)
    n_all = np.concatenate(all_n)
    k_all = np.concatenate(all_k)
    order = np.lexsort((m_all, n_all))
    starts = np.zeros(nodes + 1, dtype=np.int64)
    np.add.at(starts, n_all + 1, 1)
    return np.cumsum(starts), m_all[order], k_all[order], norm
