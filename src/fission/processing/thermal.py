"""Thermal scattering: group transfer matrices with upscattering.

Below a few eV the target nuclei are not at rest. A neutron can gain energy
in a collision (upscattering), and in water it exchanges energy with
molecular rotations and vibrations. The double-differential cross section
per atom, integrated over the azimuth, is (ENDF-102 eq. 7.6)

    sigma(E -> E', mu) = sigma_b / (2 kT) sqrt(E'/E) S_asym(alpha, beta),
    alpha = (E + E' - 2 mu sqrt(E E')) / (A kT),   beta = (E' - E) / kT,

with S_asym = exp(-beta/2) S(alpha, beta). Two models:

* FreeGas: atoms moving freely with a Maxwellian distribution,
  S_asym = exp(-(alpha + beta)^2 / (4 alpha)) / sqrt(4 pi alpha). Its total
  is the Doppler-broadened elastic cross section of chapter 8; as in NJOY,
  the kernel is normalized to the pointwise (broadened) elastic table.
* Tabulated: S(alpha, beta) from MF7 (H in H2O), absolute cross sections.

The cosine integral is done in s = sqrt(alpha), where the integrand is
smooth even at alpha -> 0:

    Int dmu sigma(E -> E', mu) P_l(mu) = sigma_b A / (2E) Int ds s S_asym(s^2, beta) P_l(mu(s)).

Detailed balance holds exactly for a symmetric S:
E exp(-E/kT) sigma(E -> E') = E' exp(-E'/kT) sigma(E' -> E).
"""

from __future__ import annotations

import math
from dataclasses import dataclass, replace

import numpy as np

from fission.endf.distributions.sab import KT0, ScatteringLaw

from .doppler import BOLTZMANN_EV_PER_K, GAUSS_RULES
from .groupr import GroupLibrary, Kernel, _legendre_matrix, group_integrals
from .groups import GroupStructure
from .weights import Weight

SQRT_4PI = math.sqrt(4.0 * math.pi)
THERMAL_CUTOFF = 4.0  # eV: an XMAS-172 boundary; MF7 data for water reach 10 eV
# Nuclides lighter than a neutron (H-1) scatter on moving targets above the
# thermal range as well, as in OpenMC (targets at rest only if awr > 1): the
# epithermal flux is higher by about 2.7 kT/E (chapter 19). A free gas up to
# 1 keV leaves less than 1e-4 of that.
LIGHT_FREE_GAS_TOP = 1000.0  # eV


@dataclass(frozen=True)
class FreeGas:
    """Free-gas scattering law for an atom of mass ratio `awr`, per unit
    free-atom cross section (sigma_free = 1 b)."""

    awr: float
    temperature: float

    @property
    def kt(self) -> float:
        return BOLTZMANN_EV_PER_K * self.temperature

    @property
    def mass(self) -> float:
        return self.awr

    @property
    def sigma_bound(self) -> float:
        return ((self.awr + 1.0) / self.awr) ** 2

    def s_asym(self, alpha: np.ndarray, beta: np.ndarray) -> np.ndarray:
        a = np.maximum(alpha, 1e-300)
        return np.exp(-((a + beta) ** 2) / (4.0 * a)) / (SQRT_4PI * np.sqrt(a))

    def describe(self) -> str:
        return f"free gas (A = {self.awr:.4f}, {self.temperature:g} K)"


@dataclass(frozen=True)
class Tabulated:
    """S(alpha, beta) from an MF7 file (per principal atom)."""

    law: ScatteringLaw

    @property
    def temperature(self) -> float:
        return self.law.temperature

    @property
    def kt(self) -> float:
        return BOLTZMANN_EV_PER_K * self.law.temperature

    @property
    def mass(self) -> float:
        return self.law.a0

    @property
    def sigma_bound(self) -> float:
        return self.law.sigma_bound

    def s_asym(self, alpha: np.ndarray, beta: np.ndarray) -> np.ndarray:
        law = self.law
        if law.lasym:
            raise NotImplementedError("asymmetric S(alpha, beta)")
        scale = self.kt / KT0 if law.lat == 1 else 1.0
        a, b = alpha * scale, np.abs(beta) * scale  # coordinates of the table
        ga, gb = law.alpha, law.beta
        inside = (a <= ga[-1]) & (b <= gb[-1])
        # ln S is linear in alpha and in beta between table points (INT = 4)
        ac = np.clip(a, ga[0], ga[-1])
        ia = np.clip(np.searchsorted(ga, ac, "right") - 1, 0, ga.size - 2)
        ib = np.clip(np.searchsorted(gb, b, "right") - 1, 0, gb.size - 2)
        ta = (ac - ga[ia]) / (ga[ia + 1] - ga[ia])
        tb = np.clip((b - gb[ib]) / (gb[ib + 1] - gb[ib]), 0.0, 1.0)
        s = law.ln_s
        ln_s = (1 - tb) * ((1 - ta) * s[ib, ia] + ta * s[ib, ia + 1]) + tb * (
            (1 - ta) * s[ib + 1, ia] + ta * s[ib + 1, ia + 1]
        )
        # below the alpha grid: the diffusive (liquid) extrapolation of NJOY
        c_liq = (s[0, 0] - s[1, 0]) * ga[0] / gb[1] ** 2 if s[0, 0] > s[0, 1] else 0.0
        small = (a < ga[0]) & (b <= 0.2) & (c_liq > 0)
        ln_s = np.where(
            small,
            s[0, 0]
            + 0.5 * np.log(ga[0] / np.maximum(a, 1e-300))
            - c_liq * b**2 / np.maximum(a, 1e-300),
            ln_s,
        )
        tab = np.exp(np.maximum(ln_s, -700.0) - beta / 2.0)
        tab = np.where(ln_s <= -225.0, 0.0, tab)
        # NJOY's rule: far out (alpha A or beta >= 30) a "negligible" table
        # value at any corner means the table has run out, not that S is 0
        corners = np.minimum(
            np.minimum(s[ib, ia], s[ib, ia + 1]),
            np.minimum(s[ib + 1, ia], s[ib + 1, ia + 1]),
        )
        exhausted = ((a * law.a0 >= 30.0) | (b >= 30.0)) & (corners <= -225.0) & ~small
        inside = inside & ~exhausted
        # outside the table: short-collision-time approximation (eq. 7.8)
        t_ratio = law.t_effective / law.temperature
        am = np.maximum(alpha, 1e-300)
        arg = (am - np.abs(beta)) ** 2 / (4.0 * am * t_ratio) + (
            np.abs(beta) + beta
        ) / 2.0
        sct = np.exp(-arg) / (SQRT_4PI * np.sqrt(am * t_ratio))
        return np.where(inside, tab, sct)

    def describe(self) -> str:
        return f"S(alpha,beta) ZA={self.law.za} at {self.law.temperature:g} K"


type Model = FreeGas | Tabulated


def ddxs(model: Model, energy: float, e_out: np.ndarray, mu: np.ndarray) -> np.ndarray:
    """sigma(E -> E', mu) per eV and unit cosine (b), for sigma_free of the model."""
    kt, a = model.kt, model.mass
    alpha = (energy + e_out - 2.0 * mu * np.sqrt(energy * e_out)) / (a * kt)
    beta = (e_out - energy) / kt
    return (
        model.sigma_bound
        / (2.0 * kt)
        * np.sqrt(e_out / energy)
        * model.s_asym(alpha, beta)
    )


def _outgoing_nodes(model: Model, energy: float, lo: float, groups: GroupStructure):
    """Gauss nodes and weights in E' covering all outgoing energies that matter."""
    kt = model.kt
    hi = energy + 40.0 * kt
    width = 2.0 * math.sqrt(energy * kt / model.mass)  # free-gas (Doppler) width
    points = {lo, hi, energy}
    points |= {b for b in groups.boundaries if lo < b < hi}
    for scale in (width, *(kt * 10.0**k for k in range(-4, 2))):
        for c in (0.1, 0.3, 1.0, 3.0, 10.0):
            for e in (energy - c * scale, energy + c * scale):
                if lo < e < hi:
                    points.add(e)
    decades = math.log10(hi / lo)
    points |= {
        lo * (hi / lo) ** (k / (60 * decades)) for k in range(int(60 * decades) + 1)
    }
    edges = np.array(sorted(points))
    t, w = (np.asarray(v) for v in GAUSS_RULES[4])
    a, b = edges[:-1], edges[1:]
    nodes = (0.5 * (a + b))[:, None] + (0.5 * (b - a))[:, None] * t[None, :]
    weights = (0.5 * (b - a))[:, None] * w[None, :]
    return nodes.ravel(), weights.ravel()


def thermal_kernel(
    model: Model, groups: GroupStructure, lmax: int, energy: float
) -> tuple[np.ndarray, float]:
    """K_l(E -> g') (shape (lmax+1, G)) and the total scattering cross section
    sigma_s(E), both per atom, for the model's cross-section scale."""
    kt, a = model.kt, model.mass
    lo = min(groups.boundaries[-1], energy) * 1e-3
    e_out, w_out = _outgoing_nodes(model, energy, lo, groups)
    beta = (e_out - energy) / kt
    s_min = np.abs(np.sqrt(energy) - np.sqrt(e_out)) / math.sqrt(a * kt)
    s_max = (np.sqrt(energy) + np.sqrt(e_out)) / math.sqrt(a * kt)
    # pieces of width <= 1 in s (the scale of the free-gas peak), 8 Gauss points each
    n_pieces = int(np.ceil(np.max(s_max - s_min)))
    t, w = (np.asarray(v) for v in GAUSS_RULES[8])
    piece = (s_max - s_min)[:, None] / n_pieces
    starts = s_min[:, None] + piece * np.arange(n_pieces)[None, :]
    s = (
        starts[:, :, None] + 0.5 * piece[:, :, None] * (t + 1.0)[None, None, :]
    ).reshape(e_out.size, -1)
    ws = np.broadcast_to(
        0.5 * piece[:, :, None] * w[None, None, :], (e_out.size, n_pieces, w.size)
    ).reshape(e_out.size, -1)
    alpha = s * s
    mu = (energy + e_out[:, None] - a * kt * alpha) / (
        2.0 * np.sqrt(energy * e_out[:, None])
    )
    f = s * model.s_asym(alpha, beta[:, None]) * ws  # (n_out, n_s)
    pref = model.sigma_bound * a / (2.0 * energy)
    pl = _legendre_matrix(lmax, np.clip(mu, -1.0, 1.0).ravel()).reshape(
        lmax + 1, *mu.shape
    )
    per_node = pref * np.einsum("lij,ij->li", pl, f) * w_out[None, :]  # (lmax+1, n_out)
    total = float(per_node[0].sum())
    asc = np.asarray(groups.boundaries[::-1])
    idx = np.searchsorted(asc, e_out, "right") - 1
    inside = (e_out >= asc[0]) & (e_out <= asc[-1])
    g_idx = groups.size - 1 - np.clip(idx, 0, groups.size - 1)
    k = np.zeros((lmax + 1, groups.size))
    for l in range(lmax + 1):
        k[l] = np.bincount(g_idx[inside], per_node[l][inside], groups.size)
    return k, total


def _model_for_rust(
    model: Model,
) -> tuple[int, list[float], list[float], list[float], list[float]]:
    """A thermal model as plain data for the Rust kernels."""
    if isinstance(model, FreeGas):
        return 0, [model.awr, model.kt, model.sigma_bound], [], [], []
    law = model.law
    if law.lasym:
        raise NotImplementedError("asymmetric S(alpha, beta)")
    scale = model.kt / KT0 if law.lat == 1 else 1.0
    params = [
        model.kt,
        law.a0,
        law.sigma_bound,
        scale,
        law.t_effective / law.temperature,
    ]
    return 1, params, law.alpha.tolist(), law.beta.tolist(), law.ln_s.ravel().tolist()


def thermal_kernels(
    model: Model, groups: GroupStructure, lmax: int, energies, use_rust: bool = True
) -> list[tuple[np.ndarray, float]]:
    """`thermal_kernel` at many incident energies (in Rust, on all cores)."""
    if not use_rust:
        return [thermal_kernel(model, groups, lmax, e) for e in energies]
    from fission import _rust

    rows = _rust.thermal_table(
        *_model_for_rust(model), list(groups.boundaries), lmax, list(energies)
    )
    g = groups.size
    return [(r[:-1].reshape(lmax + 1, g), float(r[-1])) for r in rows]


def free_gas_total(awr: float, temperature: float, energy: float) -> float:
    """Closed form for sigma_free = 1 b: (1 + 1/2y^2) erf(y) + exp(-y^2)/(y sqrt(pi)),
    y^2 = A E / kT (chapter 8)."""
    y = math.sqrt(awr * energy / (BOLTZMANN_EV_PER_K * temperature))
    return (1 + 0.5 / y**2) * math.erf(y) + math.exp(-(y**2)) / (y * math.sqrt(math.pi))


def thermal_groups(groups: GroupStructure, cutoff: float = THERMAL_CUTOFF) -> list[int]:
    return [
        g for g in range(groups.size) if groups.bounds(g)[1] <= cutoff * (1 + 1e-12)
    ]


def with_thermal(
    lib: GroupLibrary,
    model: Model,
    weight: Weight,
    elastic: tuple[list[float], list[float]] | None = None,
    cutoff: float = THERMAL_CUTOFF,
    use_rust: bool = True,
    du: float = 0.02,
    flux_shape: tuple[list[float], list[float]] | None = None,
) -> GroupLibrary:
    """The library with its thermal elastic rows replaced by a thermal kernel.

    For a free gas, `elastic` = (energies, sigma_el) of the broadened pointwise
    table normalizes the kernel; for S(alpha, beta) the cross sections are
    absolute and `elastic` is not used. With a flux shape (see
    `groupr.group_library`) the elastic table must already include it.
    The P0 rows then include upscattering,
    and sigma[2] (and the total) in the thermal groups become the thermal
    scattering cross section.
    """
    groups, lmax = lib.groups, lib.lmax
    rows = thermal_groups(groups, cutoff)
    lo, top = groups.boundaries[-1], cutoff
    seeds: set[float] = {lo, top}
    for g in rows:
        a, b = groups.bounds(g)
        n = max(
            3, math.ceil(math.log(b / a) / du)
        )  # kernel points at most du apart in lethargy
        seeds |= {a * (b / a) ** (k / n) for k in range(n + 1)}
    z = sorted(seeds)
    values = []
    for k, total in thermal_kernels(model, groups, lmax, z, use_rust):
        if isinstance(model, FreeGas):
            k = k / total  # shape only; the broadened elastic table supplies the size
            total = 1.0
        values.append(np.concatenate([k.ravel(), [total]]))
    kern = Kernel(tuple(z), np.array(values))
    if isinstance(model, FreeGas):
        if elastic is None:
            raise ValueError(
                "a free-gas kernel needs the elastic table to normalize it"
            )
        x, col = elastic
    else:
        x, col = flux_shape if flux_shape is not None else ([lo, top], [1.0, 1.0])
    out = (
        group_integrals(x, [col], weight, groups, kern, use_rust)[:, 0, :]
        / lib.flux[:, None]
    )
    g_count = groups.size
    matrix = out[:, :-1].reshape(g_count, lmax + 1, g_count).transpose(1, 0, 2)
    sigma_th = out[:, -1]

    elastic_m = lib.scatter_elastic.copy()
    sigma = {mt: v.copy() for mt, v in lib.sigma.items()}
    for g in rows:
        elastic_m[:, g, :] = matrix[:, g, :]
        if not isinstance(model, FreeGas):
            sigma[1][g] += sigma_th[g] - sigma[2][g]
            sigma[2][g] = sigma_th[g]
    notes = [*lib.notes, f"thermal: {model.describe()} below {cutoff:g} eV"]
    return replace(
        lib, sigma=sigma, scatter_elastic=elastic_m, notes=notes,
        thermal=model.describe(), thermal_cutoff=cutoff,
    )  # fmt: skip
