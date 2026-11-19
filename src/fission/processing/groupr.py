"""Group averaging: pointwise data -> multigroup constants (like NJOY's GROUPR).

Every group constant is an integral over a group of the pointwise cross
section sigma(E), a weighting spectrum w(E) and a smooth kernel K(E):

    cross section      sigma_g          = Int_g sigma w      / Int_g w        K = 1
    neutron production (nu sigma_f)_g   = Int_g sigma_f w nu / Int_g w        K = nu(E)
    fission spectrum   chi_g'           ~ Int sigma_f w nu chi(E -> g')       K = nu chi(E -> g')
    scattering matrix  sigma_l,g->g'    = Int_g sigma_s w K_l(E -> g') / Int_g w

sigma(E) has up to a million points (resonances); the kernels are smooth.
So each kernel is tabulated on its own coarse grid z_j (refined by the
midpoint test, as in chapter 7) and taken as linear in between:
K(E) = sum_j K_j h_j(E) with hat functions h_j. Then

    Int_g sigma w K = sum_j K_j A[g, j],     A[g, j] = Int_g sigma w h_j dE,

and only the A integrals see the fine sigma grid (`hat_integrals`, ported to
Rust). On each piece of the union grid sigma and h_j are linear, so four
moments of w per piece serve every reaction.

Conventions: group 0 is the highest energy; S[l, g, g'] is the transfer from
group g to g' with P_l moments of the LAB scattering cosine. Elastic
scattering is computed for a target at rest (no upscattering); thermal
scattering kernels come in a later step.
"""

from __future__ import annotations

import itertools
import math
from bisect import bisect_right
from collections.abc import Callable, Sequence
from dataclasses import dataclass, field
from functools import cache

import numpy as np

from fission.endf import Evaluation
from fission.endf.distributions.angular import (
    CM,
    AngularDistribution,
    LegendrePanel,
    angular_distribution,
)
from fission.endf.distributions.energy import TabulatedLaw, energy_distribution
from fission.endf.distributions.incident import (
    interpolate_panels,
    interpolate_parameter,
)
from fission.endf.distributions.kinematics import cm_energy, separation_energy
from fission.endf.distributions.mf6 import KALBACH, Continuum, TwoBody, energy_angle
from fission.endf.distributions.nu import delayed_groups, nu_bar
from fission.endf.interpolation import Tabulated1D
from fission.endf.pointwise import PointwiseTable
from fission.endf.reconstruct import adaptive_grid

from .doppler import GAUSS_RULES
from .groups import GroupStructure
from .weights import Weight

MAX_RATIO = 1.2  # integration pieces are split until E_hi / E_lo <= 1.2


# ------------------------------------------------------------ hat integrals
def _pieces(points: Sequence[float]):
    """Consecutive pieces of a sorted grid, split to a ratio of at most MAX_RATIO."""
    for a, b in itertools.pairwise(points):
        if b <= a:
            continue
        n = max(1, math.ceil(math.log(b / a) / math.log(MAX_RATIO))) if a > 0 else 1
        for k in range(n):
            lo = a * (b / a) ** (k / n) if a > 0 else a + (b - a) * k / n
            hi = a * (b / a) ** ((k + 1) / n) if a > 0 else a + (b - a) * (k + 1) / n
            yield lo, (b if k == n - 1 else hi)


def hat_integrals(
    x: Sequence[float],
    columns: Sequence[Sequence[float]],
    weight: Weight,
    groups: GroupStructure,
    z: Sequence[float],
) -> dict[tuple[int, int], list[float]]:
    """A[(g, j)][c] = Int_g sigma_c(E) w(E) h_j(E) dE.

    sigma_c is linear between the points x (0 outside), h_j are the hat
    functions of the increasing grid z (0 outside [z_0, z_last]).
    """
    lo = max(x[0], z[0], groups.boundaries[-1])
    hi = min(x[-1], z[-1], groups.boundaries[0])
    cuts = {
        e for e in (*x, *z, *groups.boundaries, *weight.breakpoints) if lo <= e <= hi
    }
    points = sorted(cuts)
    nodes, weights = GAUSS_RULES[8]
    acc: dict[tuple[int, int], list[float]] = {}
    for a, b in _pieces(points):
        mid = 0.5 * (a + b)
        i = min(bisect_right(x, mid) - 1, len(x) - 2)
        j = min(bisect_right(z, mid) - 1, len(z) - 2)
        g = groups.find(mid)
        if g is None:
            continue
        x0, x1, z0, z1 = x[i], x[i + 1], z[j], z[j + 1]
        m00 = m01 = m10 = m11 = 0.0  # Int w (1-u)(1-v), (1-u) v, u (1-v), u v
        for t, wt in zip(nodes, weights, strict=True):
            e = 0.5 * (a + b) + 0.5 * (b - a) * t
            f = wt * weight(e)
            u = (e - x0) / (x1 - x0)
            v = (e - z0) / (z1 - z0)
            m00 += f * (1 - u) * (1 - v)
            m01 += f * (1 - u) * v
            m10 += f * u * (1 - v)
            m11 += f * u * v
        h = 0.5 * (b - a)
        left = acc.setdefault((g, j), [0.0] * len(columns))
        right = acc.setdefault((g, j + 1), [0.0] * len(columns))
        for c, col in enumerate(columns):
            s0, s1 = col[i], col[i + 1]
            left[c] += h * (s0 * m00 + s1 * m10)
            right[c] += h * (s0 * m01 + s1 * m11)
    return acc


def _hat_integrals(x, columns, weight, groups, z, use_rust: bool):
    if use_rust:
        from fission import _rust

        gs, js, values = _rust.group_hat_integrals(
            list(x), [list(c) for c in columns], _weight_params(weight),
            list(groups.boundaries), list(z), MAX_RATIO,
        )  # fmt: skip
        return {
            (int(g), int(j)): list(v) for g, j, v in zip(gs, js, values, strict=True)
        }
    return hat_integrals(x, columns, weight, groups, z)


def _weight_params(weight: Weight) -> list[float]:
    """Weight function as numbers for Rust: [kind, ...parameters]."""
    from .weights import Constant, InverseE

    if isinstance(weight, Constant):
        return [0.0]
    if isinstance(weight, InverseE):
        return [1.0]
    return [2.0, weight.eb, weight.kt, weight.ab, weight.ec, weight.tf, weight.ac]


@dataclass(frozen=True)
class Kernel:
    """A smooth function of incident energy, vector valued, linear between z_j."""

    z: tuple[float, ...]
    values: np.ndarray  # shape (len(z), M)

    @classmethod
    def constant(cls, lo: float, hi: float, value: float = 1.0) -> Kernel:
        return cls((lo, hi), np.full((2, 1), value))


def group_integrals(
    x: Sequence[float],
    columns: Sequence[Sequence[float]],
    weight: Weight,
    groups: GroupStructure,
    kernel: Kernel,
    use_rust: bool = True,
) -> np.ndarray:
    """out[g, c, m] = Int_g sigma_c w K_m dE."""
    a = _hat_integrals(x, columns, weight, groups, kernel.z, use_rust)
    out = np.zeros((groups.size, len(columns), kernel.values.shape[1]))
    for (g, j), v in a.items():
        out[g] += np.outer(v, kernel.values[j])
    return out


def build_kernel(
    f: Callable[[float], Sequence[float] | np.ndarray],
    seeds: Sequence[float],
    tol: float = 1e-4,
    atol: float = 1e-8,
    adaptive: bool = True,
) -> Kernel:
    """Tabulate a smooth vector function of E: by midpoint refinement, or
    (adaptive=False) just at the seeds, for data that are themselves given at
    a grid of incident energies and interpolated linearly between them."""
    points = sorted(set(seeds))
    if not adaptive:
        return Kernel(
            tuple(points), np.array([list(f(e)) for e in points], dtype=float)
        )
    z, values = adaptive_grid(lambda e: tuple(f(e)), points, tol, atol, 1e-6)
    return Kernel(tuple(z), np.array(values, dtype=float))


def _seeds(
    groups: GroupStructure, lo: float, hi: float, per_group: int = 2
) -> list[float]:
    """Group boundaries in [lo, hi] plus a few points inside each group."""
    out = [lo, hi]
    for g in range(groups.size):
        a, b = groups.bounds(g)
        for k in range(per_group + 1):
            e = a * (b / a) ** (k / per_group)
            if lo <= e <= hi:
                out.append(e)
    return sorted(set(out))


# ------------------------------------------------------------ two-body kinematics
_MU, _MU_W = (np.asarray(v) for v in np.polynomial.legendre.leggauss(16))


def _legendre_matrix(lmax: int, mu: np.ndarray) -> np.ndarray:
    p = np.zeros((lmax + 1, mu.size))
    p[0] = 1.0
    if lmax:
        p[1] = mu
    for l in range(1, lmax):
        p[l + 1] = ((2 * l + 1) * mu * p[l] - l * p[l - 1]) / (l + 1)
    return p


def _panel_values(panel, mu: np.ndarray) -> np.ndarray:
    if isinstance(panel, LegendrePanel):
        a = np.array([1.0, *panel.coefficients])
        return np.polynomial.legendre.legval(mu, a * (2 * np.arange(a.size) + 1) / 2)
    t = panel.table
    return np.interp(mu, t.x, t.y)


def _angular_values(dist, energy: float, mu: np.ndarray) -> np.ndarray:
    """f(mu, E) at many cosines (MF4 or MF6 LAW=2), vectorized."""
    if isinstance(dist, AngularDistribution) and dist.grid is None:
        return np.full(mu.shape, 0.5)
    b = dist.grid.bracket(energy)
    f1 = _panel_values(dist.panels[b.i], mu)
    if b.t == 0.0:
        return f1
    f2 = _panel_values(dist.panels[min(b.i + 1, len(dist.panels) - 1)], mu)
    return (1 - b.t) * f1 + b.t * f2


@cache
def _gauss(n: int) -> tuple[np.ndarray, np.ndarray]:
    x, w = np.polynomial.legendre.leggauss(n)
    return np.asarray(x), np.asarray(w)


def _cosine_pieces(
    dist, energy: float, m1: float, m2: float
) -> tuple[list[float], int]:
    """Breaks of [m1, m2] at the tabulated cosines of the distribution, and a
    Gauss order exact for its Legendre degree (times the P_l(mu_lab) factor)."""
    breaks = [m1, m2]
    order = 16
    if not (isinstance(dist, AngularDistribution) and dist.grid is None):
        b = dist.grid.bracket(energy)
        for panel in {
            id(q): q
            for q in (dist.panels[b.i], dist.panels[min(b.i + 1, len(dist.panels) - 1)])
        }.values():
            if isinstance(panel, LegendrePanel):
                order = max(order, (len(panel.coefficients) + 1) // 2 + 4)
            else:
                breaks += [m for m in panel.table.x if m1 < m < m2]
    return sorted(set(breaks)), order


def two_body_kernel(
    dist, awr: float, q: float, groups: GroupStructure, lmax: int, energy: float
) -> np.ndarray:
    """K_l(E -> g'), shape (lmax+1, G), for a two-body reaction with a CM
    angular distribution. E'_lab is linear in mu_cm, so the cosine range
    feeding each group g' is found exactly; it is integrated by Gauss-Legendre,
    split at the cosines of a tabulated distribution."""
    k = np.zeros((lmax + 1, groups.size))
    e_cm = cm_energy(energy, awr, q)
    if e_cm <= 0.0:
        return k
    a1 = awr + 1.0
    c0 = e_cm + energy / a1**2  # E'_lab = c0 + c1 mu_cm
    c1 = 2.0 * math.sqrt(energy * e_cm) / a1
    top = groups.find(min(c0 + c1, groups.boundaries[0]))
    bottom = groups.find(max(c0 - c1, groups.boundaries[-1]))
    if top is None or bottom is None:
        return k
    for gp in range(top, bottom + 1):
        lo, hi = groups.bounds(gp)
        m1, m2 = max(-1.0, (lo - c0) / c1), min(1.0, (hi - c0) / c1)
        if m2 <= m1:
            continue
        breaks, order = _cosine_pieces(dist, energy, m1, m2)
        x, w = _gauss(order)
        for a, b in itertools.pairwise(breaks):
            mu = 0.5 * (a + b) + 0.5 * (b - a) * x
            e_lab = c0 + c1 * mu
            mu_lab = np.sqrt(e_cm / e_lab) * mu + np.sqrt(energy / e_lab) / a1
            f = _angular_values(dist, energy, mu) * 0.5 * (b - a) * w
            k[:, gp] += _legendre_matrix(lmax, np.clip(mu_lab, -1, 1)) @ f
    return k


# ------------------------------------------------------------ continuum (MF6 LAW=1)
def _continuum_panel(cont: Continuum, energy: float):
    """The distribution at an incident energy, on a grid of E' points:
    (e_out, components [n, NA+1], histogram?). Components are densities
    (f_0, Legendre f_l) or, for Kalbach-Mann, (f_0, r)."""
    b = cont.grid.bracket(energy)
    p1 = cont.panels[b.i]
    p2 = cont.panels[min(b.i + 1, len(cont.panels) - 1)]
    na = max(len(r) for r in (*p1.b, *p2.b))
    if b.unit_base:
        (a1, z1), (a2, z2) = p1.range, p2.range
        xs = {(e - a1) / (z1 - a1) for e in p1.e_out} | {
            (e - a2) / (z2 - a2) for e in p2.e_out
        }
        lo, hi = a1 + b.t * (a2 - a1), z1 + b.t * (z2 - z1)
        grid = sorted(lo + x * (hi - lo) for x in xs)
    else:
        grid = sorted(set(p1.e_out) | set(p2.e_out))
    comps = np.zeros((len(grid), na))
    for j in range(na):
        c1, c2 = p1.component(j), p2.component(j)
        interp = (
            interpolate_parameter
            if (cont.lang == KALBACH and j == 1)
            else interpolate_panels
        )
        for n, e in enumerate(grid):
            comps[n, j] = interp(c1, p1.range, c2, p2.range, b.t, b.unit_base, e)
    return np.array(grid), comps, p1.lep == 1


def _kalbach_slope_np(
    energy: float, e_out_cm: np.ndarray, za: int, awr: float
) -> np.ndarray:
    """Kalbach's slope (see kinematics.kalbach_slope) for many outgoing energies."""
    s_n = _separation_energy(za)
    e_a = energy * awr / (awr + 1.0) / 1e6 + s_n
    e_b = e_out_cm * (awr + 1.0) / awr / 1e6 + s_n
    x1 = min(e_a, 130.0) * e_b / e_a
    x3 = min(e_a, 41.0) * e_b / e_a
    return 0.04 * x1 + 1.8e-6 * x1**3 + 6.7e-7 * 1.0 * 0.5 * x3**4


@cache
def _separation_energy(za: int) -> float:
    return separation_energy(za + 1, za)


def continuum_kernel(
    cont: Continuum,
    frame: int,
    awr: float,
    groups: GroupStructure,
    lmax: int,
    energy: float,
) -> np.ndarray:
    """K_l(E -> g') for continuum emission, shape (lmax+1, G).

    In the LAB frame the Legendre components are integrated over each group
    directly. In the CM frame (E'_cm, mu_cm) is sampled with Gauss points,
    transformed to (E'_lab, mu_lab) and binned into groups; the binning keeps
    the total (P0 conservation) exact and is accurate to the resolution of
    the quadrature at group boundaries.
    """
    k = np.zeros((lmax + 1, groups.size))
    if energy < cont.grid.energies[0]:
        return k
    grid, comps, histogram = _continuum_panel(cont, energy)
    if grid[-1] <= grid[0]:
        return k
    gn, gw = GAUSS_RULES[4]
    nodes, wts = [], []
    for a, b in itertools.pairwise(grid):
        if b > a:
            nodes.extend(0.5 * (a + b) + 0.5 * (b - a) * np.asarray(gn))
            wts.extend(0.5 * (b - a) * np.asarray(gw))
    e_out, w_out = np.array(nodes), np.array(wts)
    if histogram:
        idx = np.clip(np.searchsorted(grid, e_out, "right") - 1, 0, len(grid) - 2)
        c = comps[idx]
    else:
        c = np.column_stack(
            [np.interp(e_out, grid, comps[:, j]) for j in range(comps.shape[1])]
        )
    f0 = c[:, 0]
    if frame == 1:  # LAB: Legendre components integrate directly
        if cont.lang == KALBACH:
            raise NotImplementedError("Kalbach-Mann data must be in the CM frame")
        lab_e = e_out
        moments = np.zeros((lmax + 1, e_out.size))
        for l in range(min(lmax + 1, c.shape[1])):
            moments[l] = c[:, l] if l else f0
        contrib = moments * w_out
    else:  # CM: sample cosines, transform each (E'_cm, mu_cm) to the LAB
        mu_x, mu_w = _MU, _MU_W
        mu = mu_x[None, :]
        if cont.lang == KALBACH:
            a = _kalbach_slope_np(energy, e_out, cont.za, awr)[:, None]
            r = c[:, 1][:, None]
            shape = np.where(
                a < 1e-8,
                0.5 * (1 + r * a * mu),
                a
                / (2 * np.sinh(np.maximum(a, 1e-8)))
                * (np.cosh(a * mu) + r * np.sinh(a * mu)),
            )
            f = f0[:, None] * shape
        else:
            # exact for the Legendre degree of the data (up to 64 in ENDF/B-VIII.1)
            mu_x, mu_w = _gauss(max(16, c.shape[1] // 2 + 4))
            mu = mu_x[None, :]
            p = _legendre_matrix(c.shape[1] - 1, mu_x)  # (L, n)
            coef = c * (2 * np.arange(c.shape[1]) + 1) / 2.0
            f = coef @ p
        a1 = awr + 1.0
        e_cm = e_out[:, None]
        lab_e = e_cm + energy / a1**2 + 2 * np.sqrt(energy * e_cm) * mu / a1
        lab_mu = np.sqrt(e_cm / lab_e) * mu + np.sqrt(energy / lab_e) / a1
        weight = f * w_out[:, None] * mu_w[None, :]
        pl = _legendre_matrix(lmax, np.clip(lab_mu, -1, 1).ravel()).reshape(
            lmax + 1, *lab_mu.shape
        )
        contrib = (pl * weight[None]).reshape(lmax + 1, -1)
        lab_e = lab_e.ravel()
    asc = np.asarray(groups.boundaries[::-1])
    idx = np.searchsorted(asc, lab_e, "right") - 1
    g_idx = groups.size - 1 - np.clip(idx, 0, groups.size - 1)
    inside = (lab_e >= asc[0]) & (lab_e <= asc[-1])
    for l in range(lmax + 1):
        k[l] = np.bincount(g_idx[inside], contrib[l][inside], groups.size)
    return k


# ------------------------------------------------------------ library
@dataclass
class GroupLibrary:
    """Infinitely dilute multigroup data of one nuclide at one temperature."""

    nuclide: str
    temperature: float
    groups: GroupStructure
    weight: str
    lmax: int
    flux: np.ndarray  # Int_g w dE
    sigma: dict[int, np.ndarray]  # MT -> sigma_g (b)
    scatter_elastic: np.ndarray  # (lmax+1, G, G)
    scatter_inelastic: np.ndarray  # (lmax+1, G, G), includes (n,xn) multiplicities
    nu_fission: np.ndarray | None = None  # (nu_t sigma_f)_g
    nu_prompt: np.ndarray | None = None  # nu_p per group (fission-rate weighted)
    nu_delayed: np.ndarray | None = None
    chi_prompt: np.ndarray | None = None
    chi_delayed: np.ndarray | None = None  # (families, G)
    chi: np.ndarray | None = None  # total, (nu_p chi_p + nu_d chi_d)/nu_t
    decay_constants: np.ndarray | None = None
    notes: list[str] = field(default_factory=list)
    thermal: str = "target at rest"  # elastic kernel below thermal_cutoff
    thermal_cutoff: float = 0.0

    @property
    def scatter(self) -> np.ndarray:
        return self.scatter_elastic + self.scatter_inelastic


def shape_table(
    table: PointwiseTable, energies: Sequence[float], shape: Sequence[float]
) -> PointwiseTable:
    """sigma(E) F(E) for every reaction, on the union of both grids (within
    the table's range): the input of `group_library` with a flux shape."""
    x = np.asarray(table.energies)
    sx, sf = np.asarray(energies), np.asarray(shape)
    u = np.union1d(x, sx[(sx >= x[0]) & (sx <= x[-1])])
    f = np.interp(u, sx, sf)
    return PointwiseTable(
        u.tolist(),
        {
            mt: (np.interp(u, x, np.asarray(v)) * f).tolist()
            for mt, v in table.xs.items()
        },
    )


def _fold(
    table: Tabulated1D | None, yield_: Tabulated1D | None, k: np.ndarray, e: float
):
    s = table(e) if table is not None else 1.0
    y = yield_(e) if yield_ is not None else 1.0
    return (s * y * k).ravel()


def _transfer(
    f, groups, lmax, seeds, ones, weight, use_rust, x=None, column=None, adaptive=True
):
    """Int_g [column] w K(E) for a kernel K(E) -> (lmax+1, G) matrix rows."""
    return _integrate_kernel(
        build_kernel(f, seeds, tol=1e-3, atol=1e-5, adaptive=adaptive),
        groups, lmax, ones, weight, use_rust, x, column,
    )  # fmt: skip


def _integrate_kernel(kern, groups, lmax, ones, weight, use_rust, x=None, column=None):
    if x is None:
        x, column = ones
    out = group_integrals(x, [column], weight, groups, kern, use_rust)[:, 0, :]
    g_count = groups.size
    return out.reshape(g_count, lmax + 1, g_count).transpose(1, 0, 2)


def _angular_for_rust(dist) -> tuple[list[float], list[int], list[int], list]:
    """An angular distribution as plain data for the Rust kernels."""
    if isinstance(dist, AngularDistribution) and dist.grid is None:
        return [], [], [], []
    g = dist.grid
    panels: list[tuple[int, list[float], list[float]]] = []
    for panel in dist.panels:
        if isinstance(panel, LegendrePanel):
            panels.append((0, list(panel.coefficients), []))
        else:
            panels.append((1, list(panel.table.x), list(panel.table.y)))
    return list(g.energies), list(g.breaks), list(g.laws), panels


def two_body_kernel_table(
    dist,
    awr: float,
    q: float,
    groups: GroupStructure,
    lmax: int,
    seeds: Sequence[float],
    sigma: Tabulated1D | None = None,
    adaptive: bool = True,
    use_rust: bool = True,
) -> Kernel:
    """The two-body kernel tabulated over incident energy (times sigma(E), if
    given): refined by the midpoint test, or at the seeds only."""
    points = sorted(set(seeds))
    if use_rust:
        from fission import _rust

        z, rows = _rust.two_body_table(
            *_angular_for_rust(dist), awr, q, list(groups.boundaries), lmax, points,
            sigma, adaptive, 1e-3, 1e-5,
        )  # fmt: skip
        return Kernel(tuple(z.tolist()), np.array(rows))
    return build_kernel(
        lambda e: _fold(sigma, None, two_body_kernel(dist, awr, q, groups, lmax, e), e),
        points, tol=1e-3, atol=1e-5, adaptive=adaptive,
    )  # fmt: skip


def group_library(
    ev: Evaluation,
    table: PointwiseTable,
    groups: GroupStructure,
    weight: Weight,
    *,
    nuclide: str = "",
    temperature: float = 0.0,
    lmax: int = 3,
    use_rust: bool = True,
    flux_shape: tuple[Sequence[float], Sequence[float]] | None = None,
) -> GroupLibrary:
    """Multigroup constants from a (Doppler-broadened) pointwise table.

    flux_shape = (energies, F): the flux inside each group is w(E) F(E)
    instead of w(E), e.g. the narrow-resonance flux F = 1/Sigma_t(E) of a
    mixture. The table must then already be multiplied by F (`shape_table`);
    F itself weights the flux and the smooth (folded) reactions.
    """
    lo, hi = groups.boundaries[-1], groups.boundaries[0]
    flux_k = Kernel.constant(lo, hi)
    ones = (
        ([lo, hi], [1.0, 1.0])
        if flux_shape is None
        else (list(flux_shape[0]), list(flux_shape[1]))
    )
    flux = group_integrals(ones[0], [ones[1]], weight, groups, flux_k, use_rust)[
        :, 0, 0
    ]
    mts = table.mts
    xs = group_integrals(
        table.energies, [table.xs[mt] for mt in mts], weight, groups, flux_k, use_rust
    )
    sigma = {mt: xs[:, c, 0] / flux for c, mt in enumerate(mts)}
    awr = ev.header.awr
    g_count = groups.size
    notes: list[str] = []

    # elastic: pointwise sigma (resonances) times the two-body kernel
    ad = angular_distribution(ev, 2)
    seeds = _seeds(groups, table.energies[0], min(table.energies[-1], hi))
    alpha = ((awr - 1) / (awr + 1)) ** 2
    seeds += [
        b / alpha for b in groups.boundaries if alpha > 0 and lo <= b / alpha <= hi
    ]
    s_el = _integrate_kernel(
        two_body_kernel_table(ad, awr, 0.0, groups, lmax, seeds, use_rust=use_rust),
        groups, lmax, ones, weight, use_rust, table.energies, table.xs[2],
    ) / flux[None, :, None]  # fmt: skip

    # inelastic and (n,xn): smooth MF3 cross sections folded into the kernels
    s_in = np.zeros((lmax + 1, g_count, g_count))
    for mt in sorted(m for m in ev.reactions if 51 <= m <= 90):
        cs = ev.cross_section(mt)
        dist: AngularDistribution | Continuum | TwoBody | None = None
        if (4, mt) in ev.sections:
            dist = angular_distribution(ev, mt)
            if dist.frame != CM and not dist.isotropic:
                raise NotImplementedError(f"MT{mt}: two-body data in the LAB frame")
        elif (6, mt) in ev.sections:
            products = energy_angle(ev, mt).products
            dist = products[0].distribution if products else None
        if dist is None or not isinstance(dist, AngularDistribution | TwoBody):
            dist = AngularDistribution(mt, CM, None, ())
        thr = max(cs.threshold, lo)
        level_seeds = _seeds(groups, thr, hi, 4)
        if isinstance(dist, AngularDistribution | TwoBody) and dist.grid is not None:
            level_seeds += [e for e in dist.grid.energies if thr <= e <= hi]
        kern = two_body_kernel_table(
            dist,
            awr,
            cs.qi,
            groups,
            lmax,
            level_seeds,
            cs.table,
            adaptive=False,
            use_rust=use_rust,
        )
        s_in += _integrate_kernel(kern, groups, lmax, ones, weight, use_rust)
    for mf, mt in sorted(ev.sections):
        if mf != 6 or mt == 18 or 51 <= mt <= 90 or mt not in ev.reactions:
            continue
        ea = energy_angle(ev, mt)
        cs = ev.cross_section(mt)
        for p in ea.products:
            if not isinstance(p.distribution, Continuum):
                notes.append(f"MT{mt}: neutron product with LAW={p.law} skipped")
                continue
            cont, frame = p.distribution, 1 if ea.frame == 1 else CM
            thr = max(cont.grid.energies[0], cs.threshold, lo)
            seeds = _seeds(groups, thr, hi, 2) + [
                e for e in cont.grid.energies if thr <= e <= hi
            ]
            s_in += _transfer(
                lambda e, c=cont, fr=frame, t=cs.table, y=p.yield_: _fold(t, y, continuum_kernel(c, fr, awr, groups, lmax, e), e),
                groups, lmax, seeds, ones, weight, use_rust, adaptive=False,
            )  # fmt: skip
    s_in /= flux[None, :, None]

    lib = GroupLibrary(
        nuclide,
        temperature,
        groups,
        weight.describe(),
        lmax,
        flux,
        sigma,
        s_el,
        s_in,
        notes=notes,
    )
    if 18 in table.xs and (1, 452) in ev.sections:
        _fission(ev, table, groups, weight, lib, use_rust)
    return lib


def _fission(ev, table, groups, weight, lib: GroupLibrary, use_rust: bool) -> None:
    """nu, nu sigma_f and the fission spectra, weighted by the fission source."""
    x, sf = table.energies, table.xs[18]
    lo, hi = max(x[0], groups.boundaries[-1]), min(x[-1], groups.boundaries[0])
    nt, npr, nd = (nu_bar(ev, mt) for mt in (452, 456, 455))
    chi_p = energy_distribution(ev, 18)
    dg = delayed_groups(ev)
    g_count = groups.size

    # prompt spectrum per tabulated incident energy, integrated over each group:
    # chi(E -> g') is linear in E between these energies (INT = 2)
    (part,) = chi_p.partials
    law = part.law
    if not isinstance(law, TabulatedLaw) or set(law.grid.laws) != {2}:
        raise NotImplementedError("prompt spectrum must be tabulated, lin-lin in E")
    panels = [
        np.array([_panel_integral(pnl, *groups.bounds(gp)) for gp in range(g_count)])
        for pnl in law.panels
    ]
    families = [
        np.array(_family_in_groups(dg, k, groups))
        for k in range(len(dg.decay_constants))
    ]

    def f(e: float) -> np.ndarray:
        b = law.grid.bracket(e)
        chi_e = (1 - b.t) * panels[b.i] + b.t * panels[min(b.i + 1, len(panels) - 1)]
        vp, vd = npr(e), nd(e)
        delayed = [
            vd * fam_part.probability(e) * fam
            for fam_part, fam in zip(dg.spectra.partials, families, strict=True)
        ]
        return np.concatenate([[nt(e), vp, vd], vp * chi_e, *delayed])

    seeds = _seeds(groups, lo, hi, 2) + [e for e in law.grid.energies if lo <= e <= hi]
    for nu in (nt, npr, nd):
        if nu.table is not None:
            seeds += [e for e in nu.table.x if lo <= e <= hi]
    kern = build_kernel(f, seeds, adaptive=False)
    out = group_integrals(x, [sf], weight, groups, kern, use_rust)[:, 0, :]
    fission_rate = group_integrals(
        x, [sf], weight, groups, Kernel.constant(lo, hi), use_rust
    )[:, 0, 0]
    nu_t, nu_p, nu_d = out[:, 0], out[:, 1], out[:, 2]
    with np.errstate(invalid="ignore", divide="ignore"):
        lib.nu_fission = nu_t / lib.flux
        lib.nu_prompt = np.where(fission_rate > 0, nu_p / fission_rate, 0.0)
        lib.nu_delayed = np.where(fission_rate > 0, nu_d / fission_rate, 0.0)
    prod_p = out[:, 3 : 3 + g_count].sum(axis=0)  # prompt neutrons born in g'
    lib.chi_prompt = prod_p / prod_p.sum()
    nfam = len(dg.decay_constants)
    fam = out[:, 3 + g_count :].reshape(g_count, nfam, g_count).sum(axis=0)
    lib.chi_delayed = fam / fam.sum(axis=1, keepdims=True)
    lib.chi = (prod_p + fam.sum(axis=0)) / (prod_p.sum() + fam.sum())
    lib.decay_constants = np.array(dg.decay_constants)


def _integral(dist, energy: float, a: float, b: float) -> float:
    """Int_a^b f(E -> E') dE' for an MF5 distribution (all its partials)."""
    total = 0.0
    for part in dist.partials:
        p = part.probability(energy)
        if p == 0.0:
            continue
        law = part.law
        if isinstance(law, TabulatedLaw):
            bk = law.grid.bracket(energy)
            p1, p2 = law.panels[bk.i], law.panels[min(bk.i + 1, len(law.panels) - 1)]
            if bk.unit_base:
                raise NotImplementedError("unit-base fission spectra")
            total += p * (
                (1 - bk.t) * _panel_integral(p1, a, b)
                + bk.t * _panel_integral(p2, a, b)
            )
        else:
            raise NotImplementedError("fission spectrum law")
    return total


def _panel_integral(panel, a: float, b: float) -> float:
    """Int_a^b of a lin-lin (or histogram) panel."""
    xs, ys = panel.x, panel.y
    total = 0.0
    for i in range(len(xs) - 1):
        x0, x1 = max(xs[i], a), min(xs[i + 1], b)
        if x1 <= x0:
            continue
        if panel.law == 1:
            total += ys[i] * (x1 - x0)
            continue
        s = (ys[i + 1] - ys[i]) / (xs[i + 1] - xs[i])
        y0, y1 = ys[i] + s * (x0 - xs[i]), ys[i] + s * (x1 - xs[i])
        total += 0.5 * (y0 + y1) * (x1 - x0)
    return total


def _family_in_groups(dg, k: int, groups: GroupStructure) -> list[float]:
    """Delayed family k: its spectrum integrated over each group. The LF=5
    data used here have a constant theta, so this does not depend on E."""
    law = dg.spectra.partials[k].law
    theta = law.theta(1e6)
    if any(abs(law.theta(e) - theta) > 1e-9 * theta for e in law.theta.x):
        raise NotImplementedError("energy-dependent delayed spectra")
    norm = law.g.integral()
    return [
        _panel_integral(law.g, a / theta, b / theta) / norm
        for a, b in (groups.bounds(g) for g in range(groups.size))
    ]
