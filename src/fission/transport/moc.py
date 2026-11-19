"""The method of characteristics: multigroup k-eigenvalue transport in 2D.

**Equations.** Flat-source regions i with isotropic sources (per steradian)

    Q_i,g = [ sum_g' Sigma_s,g'->g phi_i,g' + chi_i,g / k * sum_g' nu Sigma_f,i,g' phi_i,g' ] / (4 pi).

Along a segment of 2D length l in region i, a neutron moving at polar
angle theta (from the z axis) travels l / sin(theta). With
tau = Sigma_t l / sin(theta), the angular flux changes by

    delta psi = (psi_in - Q / Sigma_t) (1 - exp(-tau)),    psi_out = psi_in - delta psi

(1 - exp(-tau) is computed as -expm1(-tau), accurate for thin segments).

**Scalar flux.** Integrating psi over a region and over all directions
(derivation in docs ch. 15) gives

    phi_i,g = 4 pi Q_i,g / Sigma_t,i,g
              + 1 / (Sigma_t,i,g A_i) * sum over segments in i of W delta psi,

with direction-and-strip weight W = 2 pi w_a omega_p delta_a sin(theta_p):
w_a the azimuthal weight (sum 1 over [0, pi)), omega_p the polar weight
(sum 1 over (0, pi/2)), delta_a the track spacing. Each track is swept both
ways, and the polar angles above and below the plane are folded together,
so the weights of all directions add up to 2 pi * 2 = 4 pi.

**Boundaries.** The flux leaving a sweep becomes, in the next iteration,
the incoming flux of the sweep it links to (reflective sides); at vacuum
sides it is lost: leakage = sum W psi_out.

**Areas.** Segment lengths are rescaled region by region so that the
weighted track areas equal the exact areas (chapter 14).

**Iteration.** Power iteration: from phi and k form Q, sweep all tracks
once, update phi and then k from the ratio of total fission production. Each
sweep updates scattering and fission sources together (one inner iteration
per outer).
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

from fission.model import MacroXS, Model

from . import cmfd
from .geometry import Geometry
from .trace import Crossings, Segments, crossings, trace, track_areas
from .tracks import TrackSet

FOUR_PI = 4.0 * math.pi


@dataclass(frozen=True)
class Polar:
    """Polar quadrature over (0, pi/2): sin(theta_p) and weights (sum 1)."""

    sin: np.ndarray
    weight: np.ndarray

    @classmethod
    def gauss(cls, n: int) -> Polar:
        """Gauss-Legendre in mu = cos(theta) on (0, 1)."""
        x, w = np.polynomial.legendre.leggauss(n)
        mu = (x + 1) / 2
        return cls(np.sqrt(1 - mu * mu), w / 2)


@dataclass
class Problem:
    """Everything the sweep needs, as arrays."""

    geometry: Geometry
    tracks: TrackSet
    segments: Segments  # lengths already area-corrected
    polar: Polar
    areas: np.ndarray  # exact FSR areas (cm^2)
    material: np.ndarray  # material index of each FSR
    xs: list[MacroXS]  # by material index
    crossings: Crossings | None = None  # cell-face crossings (for CMFD currents)

    @property
    def groups(self) -> int:
        return self.xs[0].groups

    def per_region(self, name: str) -> np.ndarray:
        """A cross section by region: shape (regions, groups) or (regions, G, G)."""
        table = np.array([getattr(x, name) for x in self.xs])
        return table[self.material]

    def weights(self) -> np.ndarray:
        """W for every (track, polar): 2 pi w_a omega_p delta_a sin(theta_p)."""
        a = self.tracks.angle
        wa = np.array([g.weight * g.spacing for g in self.tracks.angles])[a]
        return 2 * math.pi * wa[:, None] * (self.polar.weight * self.polar.sin)[None, :]


def build(
    model: Model,
    geometry: Geometry,
    tracks: TrackSet,
    n_polar: int = 3,
    use_rust: bool = True,
    faces: bool = True,
) -> Problem:
    """Trace, correct the segment lengths, and map regions to cross sections.

    faces: also find where tracks cross the cell faces, so that sweeps can
    tally net currents (needed by CMFD).
    """
    names = sorted(set(geometry.materials()))
    xs = []
    for name in names:
        mat = model.materials[name]
        if mat.macroscopic is None:
            raise ValueError(f"material {name!r} has no macroscopic cross sections")
        xs.append(model.multigroup[mat.macroscopic])
    index = {n: i for i, n in enumerate(names)}
    material = np.array([index[m] for m in geometry.materials()])
    segments = trace(geometry, tracks, use_rust)
    areas = geometry.areas()
    estimated = tracks.weights @ track_areas(geometry, tracks, segments)
    segments.length = segments.length * (areas / estimated)[segments.region]
    faces_crossed = crossings(geometry, tracks, segments) if faces else None
    return Problem(
        geometry, tracks, segments, Polar.gauss(n_polar), areas, material, xs,
        faces_crossed,
    )  # fmt: skip


# ---------------------------------------------------------------- one sweep
def sweep(
    problem: Problem,
    sigma_t: np.ndarray,
    source: np.ndarray,
    psi_in: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, float, np.ndarray | None]:
    """Sweep every track both ways (Python reference).

    sigma_t, source: (regions, groups); source is Q, per steradian.
    psi_in: (tracks, 2, polar, groups) incoming flux of each sweep
    (direction 0 forward, 1 backward).
    Returns (tally of W delta psi per region, next psi_in, leakage, net
    current through every cell face in the +x / +y direction, or None if
    the problem has no face crossings). A current is sum W psi over the
    sweeps crossing the face, the same flow measure as the leakage.
    """
    tr, seg, cross = problem.tracks, problem.segments, problem.crossings
    weights = problem.weights()
    inv_sin = 1.0 / problem.polar.sin
    ratio = source / sigma_t
    tally = np.zeros_like(source)
    psi_next = np.zeros_like(psi_in)
    currents = None if cross is None else np.zeros((cross.n_faces, source.shape[1]))
    leakage = 0.0
    for t in range(len(tr)):
        regions, lengths = seg.of(t)
        at: dict[int, list[tuple[int, float]]] = {}  # position -> (face, sign)
        if cross is not None:
            for c in range(cross.start[t], cross.start[t + 1]):
                at.setdefault(int(cross.position[c]), []).append(
                    (int(cross.face[c]), float(cross.sign[c]))
                )

        n = len(regions)
        for d, forward in enumerate((True, False)):
            psi = psi_in[t, d].copy()  # (polar, groups)
            direction = 1.0 if forward else -1.0
            _tally_faces(currents, at, weights[t], 0 if forward else n, psi, direction)
            order = range(n) if forward else range(n - 1, -1, -1)
            for k in order:
                r = regions[k]
                tau = sigma_t[r][None, :] * (lengths[k] * inv_sin)[:, None]
                delta = (psi - ratio[r][None, :]) * -np.expm1(-tau)
                tally[r] += weights[t] @ delta
                psi -= delta
                _tally_faces(currents, at, weights[t], k + 1 if forward else k, psi,
                             direction)  # fmt: skip
            nxt = tr.next_fwd[t] if forward else tr.next_bwd[t]
            nfw = tr.fwd_forward[t] if forward else tr.bwd_forward[t]
            if nxt < 0:
                leakage += float((weights[t] @ psi).sum())
            else:
                psi_next[nxt, 0 if nfw else 1] = psi
    return tally, psi_next, leakage, currents


def _tally_faces(
    currents: np.ndarray | None,
    at: dict[int, list[tuple[int, float]]],
    weights: np.ndarray,
    position: int,
    psi: np.ndarray,
    direction: float,
) -> None:
    """Add the signed flow sum_p W psi to the faces crossed at this position."""
    if currents is None:
        return
    for face, sign in at.get(position, ()):
        currents[face] += direction * sign * (weights @ psi)


def _rust_sweep(problem: Problem, sigma_t, source, psi_in):
    from fission import _rust

    tr, cross = problem.tracks, problem.crossings
    faces = None
    if cross is not None:
        faces = (cross.start, cross.position, cross.face, cross.sign, cross.n_faces)
    tally, psi_next, leakage, currents = _rust.moc_sweep(
        problem.segments.region, problem.segments.length, problem.segments.start,
        tr.next_fwd, tr.fwd_forward.astype(np.int64), tr.next_bwd,
        tr.bwd_forward.astype(np.int64), np.ascontiguousarray(problem.weights().ravel()),
        problem.polar.sin, np.ascontiguousarray(sigma_t.ravel()),
        np.ascontiguousarray(source.ravel()), np.ascontiguousarray(psi_in.ravel()),
        problem.groups, faces,
    )  # fmt: skip
    return (
        tally.reshape(source.shape),
        psi_next.reshape(psi_in.shape),
        leakage,
        None if cross is None else currents.reshape(cross.n_faces, problem.groups),
    )


# ---------------------------------------------------------------- eigenvalue
@dataclass
class Solution:
    problem: Problem
    k: float
    flux: np.ndarray  # (regions, groups), scalar flux
    iterations: int
    leakage: float  # neutrons leaving through vacuum sides (per unit height)
    history: list[float] = field(default_factory=list)  # k per iteration
    coarse_iterations: list[int] = field(
        default_factory=list
    )  # CMFD outers per iteration

    def rates(self, name: str) -> np.ndarray:
        """Reaction rate per region (integrated over its area), summed over groups."""
        return self.problem.areas * (self.problem.per_region(name) * self.flux).sum(
            axis=1
        )

    def balance(self) -> float:
        """production / (absorption + leakage): equals k when converged.

        Absorption here is total - out-scattering, the removal the transport
        equation actually uses; tabulated absorption columns agree with it
        only to their printed digits.
        """
        production = self.rates("nu_fission").sum()
        removal = self.problem.per_region("total") - self.problem.per_region(
            "scatter"
        ).sum(axis=2)
        absorption = (self.problem.areas * (removal * self.flux).sum(axis=1)).sum()
        return production / (absorption + self.leakage)

    def cell_rates(self, name: str = "fission") -> np.ndarray:
        """A reaction rate summed per cell, shape (ny, nx), row 0 at the top."""
        g = self.problem.geometry
        per_cell = np.bincount(
            g.cell_of_regions(), self.rates(name), minlength=len(g.cells)
        )
        return per_cell.reshape(g.ny, g.nx)[::-1]


WINDOW = 10  # iterations over which the convergence ratio is measured


def remaining_error(changes: list[tuple[float, float]]) -> tuple[float, float]:
    """Estimated distance to the converged k and fission source.

    If successive changes shrink by a ratio rho per iteration, the remaining
    error after a change d is d (rho + rho^2 + ...) = d rho / (1 - rho).
    rho is measured over the last WINDOW iterations; unaccelerated MOC has
    rho close to 1, so a small change per iteration alone means little.
    """
    if len(changes) <= WINDOW:
        return math.inf, math.inf
    out = []
    for j in range(2):
        now, before = changes[-1][j], changes[-1 - WINDOW][j]
        if now == 0.0:
            out.append(0.0)
            continue
        rho = (now / before) ** (1 / WINDOW) if before > 0 else 1.0
        out.append(now * rho / (1 - rho) if rho < 1 else math.inf)
    return out[0], out[1]


def solve(
    problem: Problem,
    tol_k: float = 1e-6,
    tol_source: float = 1e-5,
    max_iterations: int = 10000,
    use_rust: bool = True,
    verbose: bool = False,
    cmfd: bool = True,
) -> Solution:
    """Power iteration with one transport sweep per iteration.

    cmfd: accelerate every iteration with a coarse-mesh (pin-cell) problem
    (needs the face crossings, ``build(faces=True)``); otherwise plain power
    iteration.

    Stops when the *estimated remaining* errors (``remaining_error``) of k and
    of the fission source are below tol_k and tol_source.
    """
    sigma_t = problem.per_region("total")
    if np.any(sigma_t <= 0):
        raise ValueError("every region needs a positive total cross section")
    scatter = problem.per_region("scatter")  # (regions, g, g')
    nu_fission = problem.per_region("nu_fission")
    chi = problem.per_region("chi")
    areas = problem.areas
    regions, groups = sigma_t.shape
    flux = np.ones((regions, groups))
    psi = np.zeros((len(problem.tracks), 2, len(problem.polar.sin), groups))
    k = 1.0
    production = (nu_fission * flux).sum(axis=1)
    history: list[float] = []
    changes: list[tuple[float, float]] = []  # (|dk|, source change) per iteration
    leakage = 0.0
    accelerate = cmfd and problem.crossings is not None
    if cmfd and not accelerate:
        raise ValueError("CMFD needs face crossings: build the problem with faces=True")
    coarse_iterations: list[int] = []
    cmfd_state: dict = {}
    for it in range(1, max_iterations + 1):
        in_scatter = np.einsum("rgh,rg->rh", scatter, flux)
        source = (in_scatter + chi * (production / k)[:, None]) / FOUR_PI
        if use_rust:
            tally, psi, leakage, currents = _rust_sweep(problem, sigma_t, source, psi)
        else:
            tally, psi, leakage, currents = sweep(problem, sigma_t, source, psi)
        flux = FOUR_PI * source / sigma_t + tally / (sigma_t * areas[:, None])
        if accelerate:
            assert currents is not None
            flux, psi, k_new, outer = cmfd_step(problem, flux, currents, psi, k, use_rust,
                                                cmfd_state)  # fmt: skip
            coarse_iterations.append(outer)
        else:
            k_new = k * (areas @ (nu_fission * flux).sum(axis=1)) / (areas @ production)
        new_production = (nu_fission * flux).sum(axis=1)
        # relative change of the fission source, where there is a real one
        fissile = production > 1e-6 * production.max()
        change = float(
            np.max(np.abs(new_production[fissile] / production[fissile] - 1))
        )
        changes.append((abs(k_new - k), change))
        k, production = k_new, new_production
        history.append(k)
        err_k, err_s = remaining_error(changes)
        converged = err_k < tol_k and err_s < tol_source
        if verbose and (it % (10 if accelerate else 100) == 0 or converged):
            print(f"  iteration {it:5d}: k = {k:.6f}, estimated error in k {err_k:.1e}, "
                  f"in the source {err_s:.1e}")  # fmt: skip
        if converged:
            break
    return Solution(problem, k, flux, it, leakage, history, coarse_iterations)


def cmfd_step(
    problem: Problem,
    flux: np.ndarray,
    currents: np.ndarray,
    psi: np.ndarray,
    k: float,
    use_rust: bool = True,
    state: dict | None = None,
    relax: float = cmfd.RELAX,
) -> tuple[np.ndarray, np.ndarray, float, int]:
    """Solve the coarse problem and rescale the MOC fluxes with it.

    state carries the previous correction D^ between calls (relaxation).
    Returns (flux, incoming angular fluxes, k, coarse outer iterations).
    """
    state = {} if state is None else state
    coarse = cmfd.homogenize(problem, flux)
    dtx, dhx, dty, dhy = cmfd.couplings(
        problem, coarse, currents, state.get("dh"), relax
    )
    state["dh"] = (dhx, dhy)
    op = cmfd.operator(coarse, dtx, dhx, dty, dhy)
    if use_rust:
        k_c, phi_c, outer = cmfd.solve_coarse_rust(op, coarse.flux, k)
    else:
        k_c, phi_c, outer = cmfd.solve_coarse(op, coarse.flux, k)
    groups = problem.groups
    ratio = (phi_c / coarse.flux).reshape(groups, -1)  # (groups, cells)
    cell = problem.geometry.cell_of_regions()
    flux = flux * ratio[:, cell].T
    # incoming angular fluxes: scaled by the ratio of the cell they enter
    seg_cell = cell[problem.segments.region]
    start = problem.segments.start
    entry = np.stack(
        [seg_cell[start[:-1]], seg_cell[start[1:] - 1]], axis=1
    )  # (tracks, 2)
    psi = psi * np.transpose(ratio[:, entry], (1, 2, 0))[:, :, None, :]
    return flux, psi, k_c, outer
