"""Resonance self-shielding of a heterogeneous lattice by equivalence theory.

Inside a fuel pin the flux dips at every large resonance, so the group cross
sections of chapter 10 (computed with a smooth flux) overestimate the
reaction rates. Equivalence theory [Bell & Glasstone, sections 2.8c and
8.4b] reduces the heterogeneous problem to the homogeneous one of chapter 12:

* The probability that a neutron born in the fuel makes its next collision
  outside is written in the **rational form**

      P_F->M(E) = Sigma_e / (Sigma_t,F(E) + Sigma_e),     Sigma_e = (1 - C) / l,

  with l = 4 V / S the mean chord of the fuel region and C the **Dancoff
  correction**: the fraction by which the other pins shadow a pin, i.e. the
  reduction of the current entering a black pin compared with the same pin
  alone in the moderator [Bell & Glasstone eq. 2.122].

* With that form, the narrow-resonance flux in the fuel is that of a
  homogeneous mixture whose background per absorber atom is

      sigma_0,i = ( sum_{j != i} N_j sigma_p,j + Sigma_e ) / N_i

  [Bell & Glasstone eqs. 8.87-8.88; the "equivalence principle"]: the
  other nuclides of the fuel act through their potential scattering
  sigma_p = 4 pi a^2 (a: the ENDF scattering radius), and the moderator
  outside acts through the escape cross section Sigma_e. The group
  constants of nuclide i are then averaged with phi(E) = w(E) / (sigma_t,i(E)
  + sigma_0,i), the flux shape of chapter 12.

The Dancoff correction is computed exactly for the actual lattice with a
one-group fixed-source MOC problem: the resonant region is made black, every
other region gets its potential-scattering cross section Sigma_p and the
source Sigma_p / (4 pi), so that without the black region the flux would be
1 everywhere. A pin alone in such a medium absorbs S / 4 per unit height
(the inward partial current of an isotropic unit flux); in the lattice it
absorbs (1 - C) S / 4. Both currents are computed with the same tracks, the
unshadowed one with an optically thick medium around the black pin, so that
the tracking error cancels in their ratio.

Approximations: the rational form for P_F->M, narrow resonances, flat
sources in the fuel, and no interference between the resonances of
different nuclides (each sees the others' potential scattering only).
"""

from __future__ import annotations

import math
from dataclasses import dataclass, replace
from functools import cache

import numpy as np

from fission.pincell.moving import MovingKernel
from fission.processing.cache import cached_doppler
from fission.processing.groupr import GroupLibrary
from fission.processing.groups import GroupStructure
from fission.processing.library import cached_groups
from fission.transport.geometry import Geometry
from fission.transport.moc import FOUR_PI, Problem, fixed_source

BLACK = 1.0e3  # 1/cm: a 0.4 cm pellet is 400 mean free paths thick
THICK = 1.0e3  # 1/cm: a medium that hides every other pin
SLOWING_DOWN_RANGE = (
    0.996,
    24787.5,
)  # eV, XMAS-172 boundaries around U-238's resolved range
MOVING_TOP = 1000.0  # eV: moving targets below (OpenMC's resonance-scattering default)


@cache
def potential_scattering(nuclide: str) -> float:
    """sigma_p = 4 pi a^2 (barns) from the first positive scattering radius in MF2."""
    from fission.data import EndfLibrary
    from fission.endf import Evaluation

    ev = Evaluation.from_file(EndfLibrary().neutron(nuclide))
    for isotope in ev.resonances:
        for rng in isotope.ranges:
            if rng.ap > 0:
                return 4 * math.pi * rng.ap**2  # (1e-12 cm)^2 = 1e-24 cm^2 = barn
    raise ValueError(f"{nuclide}: no scattering radius in MF2")


def surface(geometry: Geometry, material: str) -> float:
    """Total length (cm, per unit height) of the circles bounding a material."""
    total = 0.0
    for cell in geometry.cells:
        for i, m in enumerate(cell.materials[:-1]):  # rings inside the square
            if m != material:
                continue
            if i > 0 and cell.materials[i - 1] != material:
                total += 2 * math.pi * cell.radii[i - 1]
            if cell.materials[i + 1] != material:
                total += 2 * math.pi * cell.radii[i]
        if cell.materials[-1] == material:
            raise ValueError(f"{material!r} fills the cell: no closed surface")
    return total


@dataclass(frozen=True)
class Escape:
    """Equivalence data of one resonant material."""

    material: str
    dancoff: float  # C
    chord: float  # mean chord 4 V / S (cm)
    sigma_e: float  # escape cross section (1 - C) / chord (1/cm)


def dancoff(
    problem: Problem, sigma_p: dict[str, float], material: str, use_rust: bool = True
) -> Escape:
    """Dancoff correction of `material` by a black-region MOC calculation.

    problem: a tracking (``moc.tracking``) of the lattice; sigma_p: the
    potential-scattering macroscopic cross section (1/cm) of every material.
    """
    g = problem.geometry
    materials = np.array(g.materials())
    black = materials == material

    def absorbed(sigma_outside: np.ndarray) -> float:
        sigma = np.where(black, BLACK, sigma_outside)
        source = np.where(
            black, 0.0, sigma / FOUR_PI
        )  # unit flux without the black pin
        flux, _ = fixed_source(
            problem, sigma[:, None], source[:, None], use_rust=use_rust
        )
        return float((problem.areas * sigma * flux[:, 0])[black].sum())

    sigma = np.array([sigma_p[m] for m in materials])
    sigma = np.maximum(sigma, 1e-12)  # a near-void gap still needs Sigma > 0
    in_lattice = absorbed(sigma)
    # The unshadowed current S / 4, measured with the same tracks: a black pin
    # in an optically thick medium (everything arriving has psi = 1 / 4 pi).
    # Using it instead of the exact S / 4 cancels most of the tracking error.
    alone = absorbed(np.full_like(sigma, THICK))
    s = surface(g, material)
    volume = float(problem.areas[black].sum())
    c = 1.0 - in_lattice / alone
    chord = 4.0 * volume / s
    return Escape(material, c, chord, (1.0 - c) / chord)


def sigma0(densities: dict[str, float], sigma_e: float) -> dict[str, float]:
    """Background cross section (barns) of every nuclide of a resonant material."""
    potential = {n: d * potential_scattering(n) for n, d in densities.items()}
    total = sum(potential.values())
    return {n: (total - potential[n] + sigma_e) / d for n, d in densities.items()}


def shielded(
    nuclide: str, temperature: float, background: float, thermal: str
) -> GroupLibrary:
    """Group constants averaged with the NR flux w(E) / (sigma_t(E) + sigma_0)."""
    table = cached_doppler(nuclide, temperature)
    energies = np.asarray(table.energies)
    shape = 1.0 / (np.asarray(table.xs[1]) + background)
    return cached_groups(
        nuclide,
        temperature,
        thermal=thermal,
        flux_shape=(energies.tolist(), shape.tolist()),
        shape_name=f"s0-{background:.6g}",
    )


# ---------------------------------------------------------------- beyond the rational form
def escape_probability(
    problem: Problem,
    sigma_p: dict[str, float],
    material: str,
    sigmas: np.ndarray | None = None,
    use_rust: bool = True,
) -> tuple[np.ndarray, np.ndarray]:
    """P(Sigma): probability that a neutron born uniformly in `material` makes
    its first collision outside it, when the material's total cross section
    is Sigma and everything else has its potential scattering (purely
    absorbing one-group problems, flat source; Bell & Glasstone sec. 8.4a)."""
    if sigmas is None:
        sigmas = np.logspace(-3, 4, 43)
    g = problem.geometry
    materials = np.array(g.materials())
    inside = materials == material
    outside = np.maximum(np.array([sigma_p[m] for m in materials]), 1e-12)
    volume = float(problem.areas[inside].sum())
    escape = []
    for sig in sigmas:
        sigma = np.where(inside, sig, outside)
        source = np.where(
            inside, 1.0 / FOUR_PI, 0.0
        )  # one neutron per cm^3 in the region
        flux, _ = fixed_source(
            problem, sigma[:, None], source[:, None], use_rust=use_rust
        )
        collided = float((problem.areas * sig * flux[:, 0])[inside].sum())
        escape.append(1.0 - collided / volume)
    return np.asarray(sigmas), np.asarray(escape)


def ring_escape(
    problem: Problem,
    sigma_p: dict[str, float],
    material: str,
    sigmas: np.ndarray | None = None,
    use_rust: bool = True,
) -> tuple[np.ndarray, np.ndarray]:
    """E_k(Sigma): probability that a neutron born in ring k of `material`
    (labels from ``Geometry.split_rings``) makes its first collision outside
    the material, for every ring at once; shape (rings, len(sigmas)).

    One sweep per Sigma, by reciprocity: for first-flight probabilities
    V_j Sigma_j P_j->k = V_k Sigma_k P_k->j (both sides are the same double
    integral of Sigma_j Sigma_k exp(-tau) / (4 pi r^2) over the two regions,
    and the optical path tau is symmetric), so a
    source Sigma_j per unit volume in every outside region j gives in ring k
    the collision rate sum_j V_j Sigma_j P_j->k = V_k Sigma E_k. The
    volume average of E_k over the rings is ``escape_probability``.
    """
    if sigmas is None:
        sigmas = np.logspace(-3, 4, 43)
    g = problem.geometry
    labels = np.array(g.split_rings(material).materials())
    materials = np.array(g.materials())
    inside = materials == material
    rings = sorted({str(x) for x in labels[inside]}, key=lambda s: int(s.split("/")[1]))
    outside = np.maximum(np.array([sigma_p[m] for m in materials]), 1e-12)
    out = np.zeros((len(rings), len(sigmas)))
    for i, sig in enumerate(sigmas):
        sigma = np.where(inside, sig, outside)
        source = np.where(inside, 0.0, sigma / FOUR_PI)
        flux, _ = fixed_source(
            problem, sigma[:, None], source[:, None], use_rust=use_rust
        )
        collisions = problem.areas * sig * flux[:, 0]
        for k, ring in enumerate(rings):
            mask = labels == ring
            out[k, i] = collisions[mask].sum() / (sig * problem.areas[mask].sum())
    return np.asarray(sigmas), out


def collision_matrix(
    problem: Problem,
    sigma_p: dict[str, float],
    material: str,
    sigmas: np.ndarray | None = None,
    use_rust: bool = True,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """P_j->i(Sigma) between the rings of `material`: the probability that a
    neutron born uniformly in ring j makes its first collision in ring i,
    when the whole material has total cross section Sigma (everything else:
    pure absorbers Sigma_p). Shape (len(sigmas), R, R), plus the ring volumes.
    1 - sum_i P_j->i is ring j's escape E_j (``ring_escape``).
    One fixed-source sweep per (Sigma, source ring).
    """
    if sigmas is None:
        sigmas = np.logspace(-3, 4, 43)
    g = problem.geometry
    labels = np.array(g.split_rings(material).materials())
    inside = np.array(g.materials()) == material
    rings = sorted({str(x) for x in labels[inside]}, key=lambda x: int(x.split("/")[1]))
    masks = [labels == r for r in rings]
    volumes = np.array([problem.areas[m].sum() for m in masks])
    outside = np.maximum(np.array([sigma_p[m] for m in g.materials()]), 1e-12)
    out = np.zeros((len(sigmas), len(rings), len(rings)))
    for s_index, sig in enumerate(sigmas):
        sigma = np.where(inside, sig, outside)
        source = np.array([np.where(m, 1.0 / FOUR_PI, 0.0) for m in masks]).T
        flux, _ = fixed_source(
            problem,
            np.repeat(sigma[:, None], len(rings), axis=1),
            source,
            use_rust=use_rust,
        )  # fmt: skip  (one "group" per source ring)
        for j in range(len(rings)):
            collisions = problem.areas * sig * flux[:, j]
            for i, m in enumerate(masks):
                out[s_index, j, i] = collisions[m].sum() / volumes[j]
    return np.asarray(sigmas), out, volumes


@dataclass
class CoupledTable:
    """First-collision probabilities between the pellet rings and the clad,
    tabulated against both materials' total cross sections."""

    labels: list[str]  # resonant regions: "fuel/0", ..., "clad"
    volumes: np.ndarray  # (R,)
    sigma_fuel: np.ndarray  # (nf,) 1/cm
    sigma_clad: np.ndarray  # (nc,)
    collide: np.ndarray  # (nf, nc, R, R): P_j->i at [., ., j, i]

    def at(self, fuel: np.ndarray, clad: np.ndarray) -> np.ndarray:
        """P_j->i at pointwise (Sigma_fuel, Sigma_clad), bilinear in the logs:
        shape (n, R, R)."""
        x, y = np.log(self.sigma_fuel), np.log(self.sigma_clad)
        lx = np.clip(np.log(fuel), x[0], x[-1])
        ly = np.clip(np.log(clad), y[0], y[-1])
        i = np.clip(np.searchsorted(x, lx) - 1, 0, len(x) - 2)
        j = np.clip(np.searchsorted(y, ly) - 1, 0, len(y) - 2)
        tx = ((lx - x[i]) / (x[i + 1] - x[i]))[:, None, None]
        ty = ((ly - y[j]) / (y[j + 1] - y[j]))[:, None, None]
        c = self.collide
        return ((1 - tx) * (1 - ty) * c[i, j] + tx * (1 - ty) * c[i + 1, j]
                + (1 - tx) * ty * c[i, j + 1] + tx * ty * c[i + 1, j + 1])  # fmt: skip


def coupled_collisions(
    problem: Problem,
    sigma_p: dict[str, float],
    sigma_fuel: np.ndarray | None = None,
    sigma_clad: np.ndarray | None = None,
    use_rust: bool = True,
) -> CoupledTable:
    """P_j->i between the fuel rings and the clad, for every pair (Sigma_fuel,
    Sigma_clad) of a grid: one fixed-source sweep per pair, with one source
    "group" per resonant region (everything else: pure absorbers Sigma_p).

    The clad is heavy (zirconium): a neutron scattered there stays close to
    its energy, so in a fuel resonance the clad's flux dips too. Treating the
    clad as a region of the slowing-down calculation, instead of as part of
    the narrow-resonance outside, needs these probabilities at both cross
    sections, since both vary with energy."""
    if sigma_fuel is None:
        sigma_fuel = np.logspace(-3, 4, 43)
    if sigma_clad is None:
        sigma_clad = np.logspace(-2, 3, 21)
    g = problem.geometry
    labels_fsr = np.array(g.split_rings("fuel").materials())
    rings = sorted({str(x) for x in labels_fsr if str(x).startswith("fuel/")},
                   key=lambda x: int(x.split("/")[1]))  # fmt: skip
    labels = [*rings, "clad"]
    masks = [labels_fsr == x for x in labels]
    volumes = np.array([problem.areas[m].sum() for m in masks])
    base = np.maximum(np.array([sigma_p[m] for m in g.materials()]), 1e-12)
    source = np.array([np.where(m, 1.0 / FOUR_PI, 0.0) for m in masks]).T
    fuel = np.any(np.array(masks[:-1]), axis=0)
    clad = masks[-1]
    out = np.zeros((len(sigma_fuel), len(sigma_clad), len(labels), len(labels)))
    for a, sf in enumerate(sigma_fuel):
        for b, sc in enumerate(sigma_clad):
            sigma = np.where(fuel, sf, np.where(clad, sc, base))
            flux, _ = fixed_source(
                problem, np.repeat(sigma[:, None], len(labels), axis=1), source,
                use_rust=use_rust,
            )  # fmt: skip
            collisions = problem.areas[:, None] * sigma[:, None] * flux  # (fsr, j)
            for i, m in enumerate(masks):
                out[a, b, :, i] = collisions[m].sum(axis=0) / volumes
    return CoupledTable(
        labels, volumes, np.asarray(sigma_fuel), np.asarray(sigma_clad), out
    )


@dataclass
class CoupledFlux:
    """The pointwise solution of the fuel-rings + clad slowing down."""

    energy: np.ndarray  # ascending grid (eV)
    level: np.ndarray  # psi_M
    flux: dict[str, np.ndarray]  # label -> phi on the grid
    source: dict[str, np.ndarray]  # label -> collision source S


def coupled_flux(
    fuel: tuple[dict[str, float], float],
    clad: tuple[dict[str, float], float],
    table: CoupledTable,
    cell: list[tuple[float, dict[str, float], float]],
    resonant: tuple[int, int],
    energy_range: tuple[float, float] = SLOWING_DOWN_RANGE,
    use_rust: bool = True,
    moving: tuple[str, ...] = (),
    moving_top: float = MOVING_TOP,
) -> CoupledFlux:
    """Flux of every pellet ring and of the clad, slowing down together.

    fuel, clad: (number densities, temperature); cell: (volume, densities,
    temperature) of every material of the cell, `resonant` the indices of the
    fuel and the clad in it (the others deplete psi_M by their absorption).
    Narrow resonances outside `energy_range`, as in ``region_flux``: there
    every region has phi_r = E_r + (1 - E_r) Sigma_p,r / Sigma_t,r with its
    own escape E_r = 1 - sum_i P_r->i.

    moving: fuel nuclides that scatter on moving targets below `moving_top`
    (``moving_kernels``; ``march``).
    """
    mats = {"fuel": fuel, "clad": clad}
    tables = {
        (m, n): cached_doppler(n, t) for m, (dens, t) in mats.items() for n in dens
    }
    grid = np.unique(np.concatenate([np.asarray(t.energies) for t in tables.values()]))
    e_low, e_top = energy_range
    alphas = {
        key: ((mass_ratio(key[1]) - 1) / (mass_ratio(key[1]) + 1)) ** 2
        for key in tables
    }
    reach = e_top * max(1 / a for a in alphas.values()) * 1.001
    grid = np.union1d(grid, [e_low, e_top, reach])
    kernels = moving_kernels(moving, fuel[1], e_low, moving_top, use_rust)
    grid = _with_mesh(grid, kernels.values(), e_low, reach)

    def column(key, mt):
        t = tables[key]
        return mats[key[0]][0][key[1]] * np.interp(grid, np.asarray(t.energies),
                                                   np.asarray(t.xs[mt]))  # fmt: skip

    total = {m: sum(column(k, 1) for k in tables if k[0] == m) for m in mats}
    sigma_p = {
        m: sum(d * potential_scattering(n) for n, d in mats[m][0].items()) for m in mats
    }
    is_clad = np.array([x == "clad" for x in table.labels])
    collide = table.at(total["fuel"], total["clad"])  # (grid, R, R)
    escape = 1.0 - collide.sum(axis=2)
    totals = np.where(is_clad[None, :], total["clad"][:, None], total["fuel"][:, None])
    potentials = np.where(is_clad, sigma_p["clad"], sigma_p["fuel"])
    flux = escape + (1.0 - escape) * potentials[None, :] / totals  # NR everywhere
    sources = np.broadcast_to(potentials, flux.shape).copy()
    level = np.ones_like(grid)
    inside = (grid >= e_low) & (grid <= reach)
    order = np.nonzero(inside)[0][::-1]
    start = int(np.searchsorted(-grid[order], -e_top, side="right"))
    scatter = [np.where(is_clad[None, :] == (k[0] == "clad"), column(k, 2)[order][:, None], 0.0)
               for k in tables]  # fmt: skip
    power = sum(
        v * sum(d * slowing_down_power(n) for n, d in dens.items())
        for v, dens, _ in cell
    )
    share = sum(cell[i][0] for i in resonant) / power
    outside = np.zeros(len(order))
    for r, (v, dens, temp) in enumerate(cell):
        if r in resonant:
            continue
        for n, d in dens.items():
            t = cached_doppler(n, temp)
            e_t = np.asarray(t.energies)
            absorption = np.interp(grid[order], e_t, np.asarray(t.xs[1])) - np.interp(
                grid[order], e_t, np.asarray(t.xs[2]))  # fmt: skip
            outside += v * d * np.maximum(absorption, 0.0) / power
    levels = np.ones(len(order))
    psi = march(grid[order], totals[order], scatter, list(alphas.values()), escape[order],
                flux[order][:start], start, share, outside, levels, collide[order],
                table.volumes, use_rust,
                moving=[kernels.get(k[1]) if k[0] == "fuel" else None for k in tables])  # fmt: skip
    flux[order] = psi
    level[order] = levels
    # S_j from Sigma_t,i V_i (psi_i - E_i psi_M) = sum_j V_j P_ji S_j
    vol = table.volumes
    lhs = np.einsum("nji,j->nij", collide[order], vol)
    rhs = totals[order] * (psi - escape[order] * levels[:, None]) * vol[None, :]
    sources[order] = np.linalg.solve(lhs, rhs[..., None])[..., 0]
    return CoupledFlux(
        grid, level,
        {x: flux[:, r] for r, x in enumerate(table.labels)},
        {x: sources[:, r] for r, x in enumerate(table.labels)},
    )  # fmt: skip


def moving_kernels(
    nuclides: tuple[str, ...],
    temperature: float,
    e_low: float,
    e_top: float,
    use_rust: bool = True,
) -> dict[str, MovingKernel]:
    """Free-gas elastic kernels (``moving.free_gas_kernel``) of `nuclides`
    at `temperature`, for incoming energies in [e_low, e_top], from their
    0 K elastic tables."""
    from fission.endf.cache import cached_pointwise
    from fission.pincell.moving import free_gas_kernel
    from fission.processing.cache import COLD

    out = {}
    for n in nuclides:
        cold = cached_pointwise(n, **COLD)
        out[n] = free_gas_kernel(n, np.asarray(cold.energies), np.asarray(cold.xs[2]),
                                 mass_ratio(n), temperature, e_low, e_top,
                                 use_rust=use_rust)  # fmt: skip
    return out


def _with_mesh(grid, kernels, lo: float, hi: float) -> np.ndarray:
    """The grid plus every kernel's mesh energies in [lo, hi] (``march``
    needs them as grid points)."""
    for k in kernels:
        e = k.energies
        grid = np.union1d(grid, e[(e >= lo) & (e <= hi)])
    return grid


def slowing_down_power(nuclide: str) -> float:
    """xi sigma_p (barns): mean lethargy gain per collision times sigma_p."""
    a = mass_ratio(nuclide)
    alpha = ((a - 1) / (a + 1)) ** 2
    xi = 1.0 if alpha == 0 else 1.0 + alpha * math.log(alpha) / (1.0 - alpha)
    return xi * potential_scattering(nuclide)


@cache
def mass_ratio(nuclide: str) -> float:
    """AWR: the nuclide's mass in neutron masses (ENDF MF1)."""
    from fission.data import EndfLibrary
    from fission.endf import Evaluation

    return Evaluation.from_file(EndfLibrary().neutron(nuclide)).header.awr


def region_flux(
    densities: dict[str, float],
    temperature: float,
    table: tuple[np.ndarray, np.ndarray],
    slowing_down: bool = False,
    energy_range: tuple[float, float] = SLOWING_DOWN_RANGE,
    cell: list[tuple[float, dict[str, float], float]] | None = None,
    region: int = 0,
    rings: np.ndarray | None = None,
    details: dict[str, np.ndarray] | None = None,
    collisions: tuple[np.ndarray, np.ndarray] | None = None,
    use_rust: bool = True,
    moving: tuple[str, ...] = (),
    moving_top: float = MOVING_TOP,
) -> tuple[list[float], list[float]] | list[tuple[list[float], list[float]]]:
    """Flux shape of a resonant region from the slowing-down balance.

    With flat sources and narrow resonances (moderator flux 1), collisions in
    region F at energy E balance as Sigma_t phi_F = (1 - P) Sigma_p + P Sigma_t
    [Bell & Glasstone eqs. 8.83-8.85], so

        phi_F(E) = P(Sigma_t(E)) + (1 - P(Sigma_t(E))) Sigma_p / Sigma_t(E),

    with Sigma_t(E) the pointwise total of the whole mixture (all its
    resonances, so overlapping resonances of different nuclides interact)
    and P the escape probability tabulated by ``escape_probability``. The
    rational form P = Sigma_e / (Sigma_t + Sigma_e) turns this into the
    equivalence flux (Sigma_p + Sigma_e) / (Sigma_t + Sigma_e).

    slowing_down: replace the NR flux in `energy_range` by the numerical
    solution of ``march``. cell: (volume, densities, temperature) of every
    region of the cell, `region` the index of this one; it sets how fast
    absorption depletes the moderator flux (without it the moderator flux
    stays 1).

    rings: E_k(Sigma) of every ring of the region (``ring_escape``, on the
    same Sigma grid as `table`); then one shape per ring is returned. The
    region's own collision source S(E) is taken flat over its rings, and the
    reciprocity argument behind the formula above gives in ring k

        phi_k(E) = E_k psi_M + (1 - E_k) S / Sigma_t,

    with S = Sigma_p and psi_M = 1 under NR, and S recovered from the
    marched balance Sigma_t psi = (1 - P) S + P Sigma_t psi_M otherwise.
    Rings near the surface see more of the moderator (larger E_k), so their
    flux dips less in a resonance: the rim absorbs more per atom.

    collisions: (P_j->i table on the Sigma grid of `table`, ring volumes)
    from ``collision_matrix``. With slowing_down and rings, every ring then
    keeps its own scattering source (``march`` with several regions)
    instead of sharing the region's flat one: in a resonance with a large
    scattering width the source itself dips towards the centre.

    details: if given, filled with the pointwise pieces on the grid
    ("energy", "total", "source" = S / Sigma_t, "level" = psi_M, and
    "flux" or "flux/k" per ring, "source/k" with ring sources), for
    ``sph_factors``.

    moving (with slowing_down): nuclides that scatter on moving targets
    below `moving_top` (``moving_kernels``; ``march``).
    """
    tables = {n: cached_doppler(n, temperature) for n in densities}
    grid = np.unique(np.concatenate([np.asarray(t.energies) for t in tables.values()]))
    e_low, e_top = energy_range
    alphas = {n: ((mass_ratio(n) - 1) / (mass_ratio(n) + 1)) ** 2 for n in densities}
    reach = (
        e_top * max(1 / a for a in alphas.values()) * 1.001
    )  # top of the widest window
    kernels: dict[str, MovingKernel] = {}
    if slowing_down:
        grid = np.union1d(grid, [e_low, e_top, reach])
        kernels = moving_kernels(moving, temperature, e_low, moving_top, use_rust)
        grid = _with_mesh(grid, kernels.values(), e_low, reach)
    total = np.zeros_like(grid)
    for n, t in tables.items():
        total += densities[n] * np.interp(
            grid, np.asarray(t.energies), np.asarray(t.xs[1])
        )
    sigma_p = sum(d * potential_scattering(n) for n, d in densities.items())
    sigmas, escape = table
    p = np.interp(np.log(total), np.log(sigmas), escape)
    flux = p + (1.0 - p) * sigma_p / total
    source = sigma_p / total  # S / Sigma_t
    level = np.ones_like(grid)  # psi_M
    ring_flux: dict[int, np.ndarray] = {}  # marched ring fluxes, if any
    ring_source: dict[int, np.ndarray] = {}  # their S_k / Sigma_t

    def shapes() -> (
        tuple[list[float], list[float]] | list[tuple[list[float], list[float]]]
    ):
        if details is not None:
            details.update(energy=grid, total=total, source=source, level=level)
        if rings is None:
            if details is not None:
                details["flux"] = flux
            return grid.tolist(), flux.tolist()
        out = []
        for k, e_k in enumerate(rings):
            if k in ring_flux:
                phi = ring_flux[k]
            else:
                e = np.interp(np.log(total), np.log(sigmas), e_k)
                phi = e * level + (1.0 - e) * source
            if details is not None:
                details[f"flux/{k}"] = phi
                if k in ring_source:
                    details[f"source/{k}"] = ring_source[k]
            out.append((grid.tolist(), phi.tolist()))
        return out

    if not slowing_down:
        return shapes()
    # ultra-fine slowing down in [e_low, e_top]; NR outside, and as the
    # history that feeds the first windows below e_top
    inside = (grid >= e_low) & (grid <= reach)
    order = np.nonzero(inside)[0][::-1]  # descending energies
    energies = grid[order]
    scatter = [
        densities[n] * np.interp(energies, np.asarray(t.energies), np.asarray(t.xs[2]))
        for n, t in tables.items()
    ]
    start = int(np.searchsorted(-energies, -e_top, side="right"))
    share, outside = 0.0, None
    if cell is not None:
        power = sum(
            v * sum(d * slowing_down_power(n) for n, d in dens.items())
            for v, dens, _ in cell
        )
        share = cell[region][0] / power
        outside = np.zeros_like(energies)
        for r, (v, dens, temp) in enumerate(cell):
            if r == region:
                continue
            for n, d in dens.items():
                t = cached_doppler(n, temp)
                e_t = np.asarray(t.energies)
                absorption = np.interp(energies, e_t, np.asarray(t.xs[1])) - np.interp(
                    energies, e_t, np.asarray(t.xs[2])
                )
                outside += v * d * np.maximum(absorption, 0.0) / power
    levels = np.ones_like(energies)
    if rings is not None and collisions is not None:
        # every ring with its own source: NR history, then the coupled march
        matrix, volumes = collisions
        log_t = np.log(total[order])
        escape_k = np.array([np.interp(log_t, np.log(sigmas), e_k) for e_k in rings]).T
        collide = np.empty((len(order), len(rings), len(rings)))
        for j in range(len(rings)):
            for i in range(len(rings)):
                collide[:, j, i] = np.interp(log_t, np.log(sigmas), matrix[:, j, i])
        nr = escape_k + (1.0 - escape_k) * (sigma_p / total[order])[:, None]
        psi_rings = march(
            energies, total[order], scatter, [alphas[n] for n in tables], escape_k,
            nr[:start], start, share, outside, levels, collide, volumes, use_rust,
            moving=[kernels.get(n) for n in tables],
        )  # fmt: skip
        # S_j / Sigma_t from Sigma_t V_i (psi_i - E_i psi_M) = sum_j V_j P_ji S_j
        lhs = np.einsum("nji,j->nij", collide, volumes)
        rhs = (psi_rings - escape_k * levels[:, None]) * volumes[None, :]
        sources = np.linalg.solve(lhs, rhs[..., None])[..., 0]
        level[order] = levels
        for k in range(len(rings)):
            e = np.interp(np.log(total), np.log(sigmas), rings[k])
            ring_flux[k] = e * level + (1.0 - e) * source
            ring_flux[k][order] = psi_rings[:, k]
            ring_source[k] = source.copy()
            ring_source[k][order] = sources[:, k]
        flux = flux.copy()
        flux[order] = psi_rings @ (volumes / volumes.sum())
        return shapes()
    psi = march(
        energies, total[order], scatter, [alphas[n] for n in tables], p[order],
        flux[order][:start], start, share, outside, levels, use_rust=use_rust,
        moving=[kernels.get(n) for n in tables],
    )  # fmt: skip
    flux = flux.copy()
    flux[order] = psi
    level[order] = levels
    pk = p[order]
    source[order] = np.where(
        pk < 1.0 - 1e-9, (psi - pk * levels) / np.maximum(1.0 - pk, 1e-9), psi
    )  # S / Sigma_t from Sigma_t psi = (1 - P) S + P Sigma_t psi_M
    return shapes()


# ---------------------------------------------------------------- SPH


def sph_factors(
    problem: Problem,
    sigma_p: dict[str, float],
    names: np.ndarray,
    energy: np.ndarray,
    level: np.ndarray,
    references: dict[str, tuple[np.ndarray, np.ndarray]],
    sigma_t: dict[str, np.ndarray],
    groups: GroupStructure,
    energy_range: tuple[float, float] = SLOWING_DOWN_RANGE,
    tol: float = 1e-8,
    max_iterations: int = 100,
    use_rust: bool = True,
) -> dict[str, np.ndarray]:
    """Superhomogenization factors mu_g of the regions of a resonant material.

    The group constant Sigma_g = <Sigma phi> / <phi> preserves the reaction
    rate only if the multigroup calculation reproduces the group flux <phi>.
    It does not: in the self-shielding model the region's flux

        phi(E) = P(Sigma_t(E)) psi_M + (1 - P(Sigma_t(E))) S / Sigma_t(E)

    depends non-linearly on Sigma_t, and the same model with one constant
    cross section per group gives another flux. SPH (superhomogenization)
    multiplies every cross section of region k in group g by mu_k,g, chosen
    so that the multigroup version of the model reproduces the reference:

        mu_k,g phi^MG_k,g(mu) = <phi_k>_g.

    The multigroup flux is then <phi_k>_g / mu_k,g and every rate mu Sigma_g
    <phi>_g / mu is preserved [Hebert, Applied Reactor Physics, SPH method].

    The multigroup model is the one-group fixed-source problem behind P
    (``escape_probability``, ``ring_escape``) for every group at once: the
    regions of `material` get mu_k,g Sigma_k,g and the group-averaged flat
    source <S>_g; everything else is a pure absorber Sigma_p with source
    Sigma_p <psi_M>_g. Then phi^MG_k is its flux in region k, and mu is
    iterated as mu <- <phi_k> / phi^MG_k(mu) (the usual SPH iteration). With
    a single region this is phi^MG = P(mu Sigma) <psi_M> + (1 - P) <S> / (mu
    Sigma); with rings the rings shadow each other, which is why they must be
    solved together. Averages are in lethargy (the weighting function is 1/E
    throughout `energy_range`); groups outside it get mu = 1.

    names: the region label of every FSR (``Geometry.split_rings`` labels
    for rings); energy, level: the pointwise grid and psi_M; references:
    label -> (pointwise flux phi_k, collision source S_k) of every resonant
    region; sigma_t: label -> group total (1/cm, without SPH). Regions not
    in `references` are the pure absorbers Sigma_p of the model material.
    """
    g_all = problem.geometry
    labels = sorted(references)
    e = energy
    lo_range, hi_range = energy_range
    chosen = [
        g for g in range(groups.size)
        if groups.bounds(g)[0] >= lo_range * 0.999
        and groups.bounds(g)[1] <= hi_range * 1.001
    ]  # fmt: skip
    if not chosen:
        return {label: np.ones(groups.size) for label in labels}

    def average(values: np.ndarray) -> np.ndarray:
        out = []
        for g in chosen:
            lo, hi = groups.bounds(g)
            # the group's grid points plus its edges (linear interpolation,
            # the grid's own law), so that narrow groups are covered too
            points = np.concatenate([[lo], e[(e > lo) & (e < hi)], [hi]])
            du = np.diff(np.log(points))
            v = np.interp(points, e, values)
            out.append(np.sum(0.5 * du * (v[1:] + v[:-1])) / du.sum())
        return np.array(out)

    level = average(level)
    sources = {label: average(references[label][1]) for label in labels}
    reference = {label: average(references[label][0]) for label in labels}
    masks = {label: names == label for label in labels}
    base = np.maximum(np.array([sigma_p.get(m, 0.0) for m in g_all.materials()]), 1e-12)
    sig = np.repeat(base[:, None], len(chosen), axis=1)
    q = base[:, None] * level[None, :]  # outside: Sigma_p psi_M
    for label in labels:
        q[masks[label]] = sources[label]
    mu = {label: np.ones(len(chosen)) for label in labels}
    # boundary fluxes carried from one iteration to the next: mu changes
    # little, so each solve starts close to its answer
    psi = np.zeros((len(problem.tracks), 2, len(problem.polar.sin), len(chosen)))
    for _ in range(max_iterations):
        for label in labels:
            sig[masks[label]] = mu[label] * sigma_t[label][chosen]
        flux, _ = fixed_source(problem, sig, q / FOUR_PI, use_rust=use_rust,
                               boundary=psi)  # fmt: skip
        change = 0.0
        for label in labels:
            m = masks[label]
            a = problem.areas[m]
            phi = (a[:, None] * flux[m]).sum(axis=0) / a.sum()
            new = reference[label] / phi
            change = max(change, float(np.max(np.abs(new / mu[label] - 1))))
            mu[label] = new
        if change < tol:
            break
    out = {}
    for label in labels:
        full = np.ones(groups.size)
        full[chosen] = mu[label]
        out[label] = full
    return out


def with_sph(lib: GroupLibrary, mu: np.ndarray) -> GroupLibrary:
    """A library with every cross section of group g multiplied by mu_g
    (rows of the transfer matrices: scattering out of g)."""
    rows = mu[None, :, None]
    return replace(
        lib,
        sigma={mt: v * mu for mt, v in lib.sigma.items()},
        scatter_elastic=lib.scatter_elastic * rows,
        scatter_inelastic=lib.scatter_inelastic * rows,
        nu_fission=None if lib.nu_fission is None else lib.nu_fission * mu,
        notes=[*lib.notes, "SPH factors applied"],
    )


# ---------------------------------------------------------------- slowing down numerically


def march(
    energies: np.ndarray,
    total: np.ndarray,
    scatter: list[np.ndarray],
    alphas: list[float],
    escape: np.ndarray,
    history: np.ndarray,
    start: int,
    share: float = 0.0,
    outside_absorption: np.ndarray | None = None,
    levels: np.ndarray | None = None,
    collide: np.ndarray | None = None,
    volumes: np.ndarray | None = None,
    use_rust: bool = False,
    moving: list[MovingKernel | None] | None = None,
    tol: float = 1e-10,
    max_passes: int = 200,
) -> np.ndarray:
    """Solve the slowing-down balance of a region on an ultra-fine grid.

    In lethargy u (flux per unit lethargy psi, asymptotically 1), with the
    moderator outside treated as narrow-resonance (no resonances, its flux
    psi_M smooth):

        Sigma_t(u) psi(u) = (1 - P(u)) S(u) + P(u) Sigma_t(u) psi_M(u),

        S(u) = sum_i 1 / (1 - alpha_i) Int_{u - Delta_i}^{u} Sigma_s,i(u') psi(u') e^{-(u - u')} du',

    the exact elastic source of each nuclide i of the region (isotropic in the
    centre of mass, target at rest; Delta_i = ln(1 / alpha_i)) [Bell &
    Glasstone sec. 8.4c]. Unlike NR or WR, scattering inside a resonance is
    followed: a neutron scattered by U-238 stays near the resonance.

    energies: descending; total, scatter[i] (= N_i sigma_s,i), escape (P) on
    that grid. psi is given (``history``) for the first `start` points
    (energies above the range) and marched below. Each nuclide keeps the
    cumulative integral C_i(u) = Int f_i e^{u'} du' (trapezoid rule), so that
    its window integral is C_i(u) - C_i(u - Delta_i); the unknown psi(u)
    enters linearly through the trapezoid's last half-step. Pure scattering
    conserves neutrons, so nothing damps a bias: both ends of every window
    must be integrated with the same rule.

    psi_M is 1 at the top and falls as absorption removes neutrons from the
    slowing-down density q of the cell: with q = Q psi_M and Q = sum_r V_r
    xi_r Sigma_p,r the cell's slowing-down power,

        d psi_M / du = -( share Sigma_a(u) psi(u) + outside_absorption(u) psi_M ),

    share = V_region / Q and outside_absorption = sum_{r outside} V_r
    Sigma_a,r(u) / Q (absorption of the other regions, e.g. boron in the
    water). With share = 0 and no outside absorption psi_M stays 1, Bell &
    Glasstone's asymptotic normalization. levels: if given (length n), filled
    with psi_M at every point (1 on the history).

    Several regions (pellet rings, or rings and the clad; volumes V_i): escape
    has shape (n, R) (E_i, ``ring_escape``) and collide (n, R, R) holds
    P_j->i, the probability that a neutron born in ring j first collides in
    ring i (``collision_matrix``). Every ring keeps its own scattering
    source S_i, and the balance becomes a small linear system per point,

        Sigma_t,i V_i psi_i = sum_j V_j P_j->i S_j + Sigma_t,i V_i E_i psi_M,

    which for one region is the equation above (P_1->1 = 1 - P). history
    and the result then have shape (n, R); absorption depletes psi_M with
    the volume-averaged absorption. Regions of different compositions (fuel
    and clad) give total and scatter[i] shape (n, R): each nuclide's
    macroscopic scattering in each region (zero where it is absent); a 1D
    array means the same in every region.

    moving: per nuclide, None (target at rest) or a ``MovingKernel``: below
    its top energy the nuclide scatters on moving targets. Its source is then

        S(u) = E Int Sigma_s(E') phi(E') P(E' -> E) dE'
             = E sum_m P(E_m -> E) Phi_m,     Phi_m = Int hat_m Sigma_s psi du',

    with P on the kernel's mesh (hat functions in x = sqrt(A E / kT), the
    mesh energies must be grid points) and S interpolated linearly in x
    between the nodes. A moving target can also *raise* the neutron's
    energy, so the source at u depends on psi below u, which the march has
    not reached yet. The march is therefore repeated (Gauss-Seidel): the
    first pass treats every target at rest, and every later pass takes Phi
    from this pass where the march has passed and from the previous pass
    elsewhere, until psi changes by less than `tol` (relative). The at-rest
    window of that nuclide then covers only scattering from above the top.

    use_rust: the same march in Rust (``slowing.rs``, ~1e-12 parity).
    """
    if use_rust:
        return _rust_march(energies, total, scatter, alphas, escape, history, start,
                           share, outside_absorption, levels, collide, volumes,
                           moving, tol, max_passes)  # fmt: skip
    single = collide is None
    if single:
        collide = (1.0 - escape)[:, None, None]
        escape = escape[:, None]
        volumes = np.ones(1)
    regions = escape.shape[1]
    totals = _per_region(total, regions)
    scatters = [_per_region(sc, regions) for sc in scatter]
    hist = np.asarray(history, dtype=float).reshape(start, regions)
    kernels = [None] * len(alphas) if moving is None else list(moving)
    targets = [_Targets(k, energies, start) if k is not None else None for k in kernels]
    passes = [None] if not any(targets) else [None, *range(max_passes)]
    phi = [None if t is None else np.zeros((t.kernel.nodes, regions)) for t in targets]
    psi = None
    for p in passes:
        rest = p is None  # the first pass treats every target at rest
        new, phi_new = _march_pass(energies, totals, scatters, alphas, escape, hist,
                                   start, share, outside_absorption, levels, collide,
                                   volumes, targets, phi, rest, single)  # fmt: skip
        if psi is not None:
            change = float(np.max(np.abs(new - psi) / np.maximum(np.abs(new), 1e-300)))
            psi, phi = new, phi_new
            if change < tol:
                break
        else:
            psi, phi = new, phi_new
    else:
        if any(targets):
            raise RuntimeError(
                f"moving-target iterations: no convergence in {max_passes}"
            )
    assert psi is not None
    return psi[:, 0] if single else psi


class _Targets:
    """Where every march point sits on the mesh of a ``MovingKernel``."""

    def __init__(self, kernel: MovingKernel, energies: np.ndarray, start: int):
        self.kernel = kernel
        # node coordinate of every point (x - base) / dx
        self.p = (np.sqrt(kernel.scale * energies) - kernel.base) / kernel.dx
        # segments [k-1, k] below the top of the range: their scattering goes
        # through the kernel (Phi), the ones above through the at-rest window
        top = kernel.e_top * (1.0 + 1e-12)
        self.on = (energies >= kernel.e_top * (1.0 - 1e-12)).astype(float)
        self.below = np.zeros(len(energies), dtype=bool)
        self.below[1:] = energies[:-1] <= top
        self.below[:start] = False
        n = kernel.nodes
        self.cell = np.clip(np.floor(self.p + 1e-9).astype(np.int64), 0, n - 2)
        for k in np.nonzero(self.below)[0]:
            if self.p[k - 1] - self.cell[k] > 1.0 + 1e-6:
                raise ValueError("the kernel's mesh energies must be grid points")


def _march_pass(energies, totals, scatters, alphas, escape, hist, start, share,
                outside_absorption, levels, collide, volumes, targets, phi_old,
                rest, single):  # fmt: skip
    """One march down the grid (``march``); returns psi (n, R) and the hat
    moments Phi of every moving nuclide accumulated on the way."""
    u = np.log(energies[0] / energies)
    n = len(energies)
    regions = escape.shape[1]
    v = np.asarray(volumes, dtype=float)
    weights = v / v.sum()
    psi = np.zeros((n, regions))
    psi[:start] = hist
    deltas = [float(np.log(1.0 / a)) for a in alphas]
    if u[start - 1] < max(deltas):
        raise ValueError(
            "the history must span the widest scattering window (ln 1/alpha)"
        )
    factors = [1.0 / (1.0 - a) for a in alphas]
    # on[i][k]: segment [k-1, k] scatters at rest (1) or through the kernel (0)
    on = [np.ones(n) if rest or t is None else t.on for t in targets]
    cumulative = [np.zeros((n, regions)) for _ in alphas]
    f = [sc * psi * np.exp(u)[:, None] for sc in scatters]  # valid on the history
    for c, fi in zip(cumulative, f, strict=True):
        c[1:start] = np.cumsum(
            0.5 * np.diff(u[:start])[:, None] * (fi[1:start] + fi[: start - 1]), axis=0
        )
    phi_new = [
        None if t is None else np.zeros_like(p)
        for t, p in zip(targets, phi_old, strict=True)
    ]
    sources = [None if t is None else np.full((t.kernel.nodes, regions), np.nan)
               for t in targets]  # fmt: skip
    pointer = [0] * len(alphas)  # first grid index with u >= u_k - Delta_i
    level = 1.0  # psi_M
    if levels is not None:
        levels[:start] = 1.0
    for k in range(start, n):
        h = u[k] - u[k - 1]
        eu = np.exp(u[k])
        known = np.zeros(regions)  # S_r without the unknown psi_r(u_k)
        coefficient = np.zeros(regions)  # d S_r / d psi_r(u_k)
        for i, delta in enumerate(deltas):
            lower = u[k] - delta
            j = pointer[i]
            while u[j + 1] < lower:
                j += 1
            pointer[i] = j
            c = cumulative[i]
            o = on[i]
            if j + 1 < k:
                # C at the window's lower edge with the same trapezoid rule as
                # C itself (the integrand interpolated linearly), so that both
                # ends of the window are treated alike and their errors cancel
                t = (lower - u[j]) / (u[j + 1] - u[j])
                g_lower = f[i][j] + t * (f[i][j + 1] - f[i][j])
                c_lower = c[j] + o[j + 1] * 0.5 * (lower - u[j]) * (f[i][j] + g_lower)
                window = c[k - 1] + o[k] * 0.5 * h * f[i][k - 1] - c_lower
                weight_k = o[k] * 0.5 * h  # of f_k in the window integral
            else:
                # the whole window lies inside the last step (a very heavy
                # nuclide): the trapezoid from the lower edge to u_k, with the
                # integrand at the edge interpolated between f_{k-1} and f_k
                t = (lower - u[k - 1]) / h
                width = u[k] - lower
                window = o[k] * 0.5 * width * (1.0 - t) * f[i][k - 1]
                weight_k = o[k] * 0.5 * width * (1.0 + t)
            known += factors[i] * window / eu
            coefficient += factors[i] * weight_k * scatters[i][k]
            tg = targets[i]
            if tg is not None and not rest:
                known += energies[k] * _moving_source(tg, k, phi_new[i], phi_old[i],
                                                      sources[i])  # fmt: skip
        # the moderator level, depleted by absorption up to the previous point
        if share > 0.0 or outside_absorption is not None:
            absorption = totals[k - 1] - sum(sc[k - 1] for sc in scatters)
            absorbed = share * float(weights @ (absorption * psi[k - 1]))
            if outside_absorption is not None:
                absorbed += outside_absorption[k - 1] * level
            level -= h * absorbed
        if levels is not None:
            levels[k] = level
        if single:
            p = escape[k, 0]
            psi[k, 0] = ((1.0 - p) * known[0] + p * totals[k, 0] * level) / (
                totals[k, 0] - (1.0 - p) * coefficient[0]
            )
        else:
            # Sigma_t,i V_i psi_i - sum_j V_j P_ji coef_j psi_j
            #     = sum_j V_j P_ji known_j + Sigma_t,i V_i E_i psi_M
            transfer = collide[k].T * v[None, :]  # [i, j] = V_j P_j->i
            a = np.diag(totals[k] * v) - transfer * coefficient[None, :]
            b = transfer @ known + totals[k] * v * escape[k] * level
            psi[k] = np.linalg.solve(a, b)
        for i in range(len(deltas)):
            f[i][k] = scatters[i][k] * psi[k] * eu
            cumulative[i][k] = cumulative[i][k - 1] + on[i][k] * 0.5 * h * (
                f[i][k - 1] + f[i][k]
            )
            tg = targets[i]
            if tg is not None and tg.below[k]:
                # hat moments Phi_m = Int hat_m Sigma_s psi du of this segment
                j, p0, p1 = tg.cell[k], tg.p[k - 1], tg.p[k]
                g0, g1 = scatters[i][k - 1] * psi[k - 1], scatters[i][k] * psi[k]
                t0, t1 = p0 - j, p1 - j
                phi_new[i][j] += 0.5 * h * ((1.0 - t0) * g0 + (1.0 - t1) * g1)
                phi_new[i][j + 1] += 0.5 * h * (t0 * g0 + t1 * g1)
    for i, tg in enumerate(targets):
        if tg is not None:
            _below_the_range(tg, energies[0], scatters[i][-1] * psi[-1], phi_new[i])
    return psi, phi_new


def _below_the_range(tg: _Targets, e0: float, g: np.ndarray, phi: np.ndarray) -> None:
    """Hat moments below the last grid point, with Sigma_s psi continued flat.

    Neutrons below the march (e.g. under 1 eV) also scatter up into it;
    without them the flux would dip at the bottom of the range. The segments
    run from the last point down through the mesh nodes, with the same
    trapezoid rule as the march's own segments.
    """
    kernel = tg.kernel
    top = tg.p[-1]
    node = math.floor(top + 1e-9)
    if abs(top - round(top)) < 1e-9:
        node = round(top) - 1
    while node >= 0:
        e_hi = (kernel.base + top * kernel.dx) ** 2 / kernel.scale
        e_lo = (kernel.base + node * kernel.dx) ** 2 / kernel.scale
        h = math.log(e0 / e_lo) - math.log(e0 / e_hi)
        t0 = top - node
        phi[node] += 0.5 * h * ((1.0 - t0) * g + g)
        phi[node + 1] += 0.5 * h * (t0 * g)
        top, node = float(node), node - 1


def _moving_source(tg: _Targets, k: int, phi_new: np.ndarray, phi_old: np.ndarray,
                   cache: np.ndarray) -> np.ndarray:  # fmt: skip
    """s(E_k) = sum_m P(E_m -> E_k) Phi_m per unit energy, interpolated
    linearly in x between the mesh nodes around E_k (zero off the mesh).

    Gauss-Seidel: hats the march has already passed (x_{m-1} >= x_{k-1}) use
    this pass's Phi, the others the previous pass's. s_n is computed the
    first time the march needs it and kept for the pass.
    """
    kernel = tg.kernel
    p = tg.p[k]
    if p < 0.0 or p > kernel.nodes - 1:
        return np.zeros(phi_new.shape[1])
    j = int(tg.cell[k])
    theta = max(p - j, 0.0)
    above = tg.p[k - 1]
    q = math.floor(above + 1e-9)
    done = q + 1 if abs(above - round(above)) < 1e-9 else q + 2  # first complete hat
    out = np.zeros(phi_new.shape[1])
    for node, weight in ((j, 1.0 - theta), (j + 1, theta)):
        if np.isnan(cache[node, 0]):
            rows, values = kernel.column(node)
            current = np.where((rows >= done)[:, None], phi_new[rows], phi_old[rows])
            cache[node] = values @ current
        out += weight * cache[node]
    return out


def _per_region(a: np.ndarray, regions: int) -> np.ndarray:
    """(n,) -> (n, R) by repetition; (n, R) unchanged."""
    a = np.asarray(a, dtype=float)
    return np.repeat(a[:, None], regions, axis=1) if a.ndim == 1 else a


def _rust_march(energies, total, scatter, alphas, escape, history, start, share,
                outside_absorption, levels, collide, volumes, moving=None, tol=1e-10,
                max_passes=200):  # fmt: skip
    from fission import _rust

    single = collide is None
    regions = 1 if single else escape.shape[1]

    def flat(a) -> np.ndarray:
        return np.ascontiguousarray(a, dtype=float).ravel()

    totals = _per_region(total, regions)
    scatters = np.array([_per_region(sc, regions) for sc in scatter])
    psi, level = _rust.slowing_down_march(
        flat(energies), flat(totals), flat(scatters), flat(alphas),
        flat(escape), flat([] if single else collide),
        flat(np.ones(1) if single else volumes), flat(history), start, share,
        flat([] if outside_absorption is None else outside_absorption),
        [] if moving is None else list(moving), tol, max_passes,
    )  # fmt: skip
    if levels is not None:
        levels[:] = level
    psi = psi.reshape(len(energies), regions)
    return psi[:, 0] if single else psi
