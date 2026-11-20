"""Type stubs for the Rust extension module (rust/src/lib.rs)."""

from collections.abc import Sequence
from typing import Any

import numpy as np
import numpy.typing as npt

__version__: str

type Array = npt.NDArray[np.float64]

class ReichMoore:
    def __init__(
        self, spin_groups: Sequence[object], potentials: Sequence[object], spi: float
    ) -> None: ...
    def __call__(self, energy: float) -> tuple[float, float, float]: ...
    def evaluate(self, energies: Array) -> tuple[Array, Array, Array]: ...

class Mlbw:
    def __init__(
        self, spin_groups: Sequence[object], potentials: Sequence[object], spi: float
    ) -> None: ...
    def __call__(self, energy: float) -> tuple[float, float, float]: ...
    def evaluate(self, energies: Array) -> tuple[Array, Array, Array]: ...

def reconstruct(
    model: ReichMoore | Mlbw,
    backgrounds: Sequence[object | None],
    seeds: Sequence[float],
    tol: float = ...,
    atol: float = ...,
    min_rel_width: float = ...,
) -> tuple[Array, Array, Array, Array]: ...
def tab1_eval(table: object, xs: Array) -> Array: ...
def k_constant() -> float: ...
def doppler_weights(x: Sequence[float], y: float) -> tuple[int, Array]: ...
def doppler_evaluate(
    energies: Sequence[float],
    columns: Sequence[Sequence[float]],
    awr: float,
    temperature: float,
    at: Sequence[float],
) -> list[Array]: ...
def doppler_thin(
    energies: Sequence[float],
    columns: Sequence[Sequence[float]],
    tol: float,
    atol: float,
) -> list[int]: ...
def doppler_broaden_table(
    energies: Sequence[float],
    columns: Sequence[Sequence[float]],
    awr: float,
    temperature: float,
    seeds: Sequence[float],
    tol: float = ...,
    atol: float = ...,
    min_rel_width: float = ...,
) -> tuple[Array, list[Array]]: ...
def pointwise_table(
    model: ReichMoore | Mlbw | None,
    el: float,
    eh: float,
    resonance_slots: Sequence[int | None],
    tables: Sequence[object | None],
    total: int,
    partials: Sequence[int],
    seeds: Sequence[float],
    unsplittable: Sequence[tuple[float, float]],
    tol: float = ...,
    atol: float = ...,
    min_rel_width: float = ...,
) -> tuple[Array, list[Array]]: ...
def fsum(values: Sequence[float]) -> float: ...
def group_hat_integrals(
    x: Sequence[float],
    columns: Sequence[Sequence[float]],
    weight: Sequence[float],
    boundaries: Sequence[float],
    z: Sequence[float],
    max_ratio: float,
) -> tuple[list[int], list[int], list[list[float]]]: ...
def two_body_kernel(
    energies: Sequence[float],
    breaks: Sequence[int],
    laws: Sequence[int],
    panels: Sequence[tuple[int, Sequence[float], Sequence[float]]],
    awr: float,
    q: float,
    boundaries: Sequence[float],
    lmax: int,
    energy: float,
) -> Array: ...
def two_body_table(
    energies: Sequence[float],
    breaks: Sequence[int],
    laws: Sequence[int],
    panels: Sequence[tuple[int, Sequence[float], Sequence[float]]],
    awr: float,
    q: float,
    boundaries: Sequence[float],
    lmax: int,
    seeds: Sequence[float],
    sigma: object | None = None,
    adaptive: bool = True,
    tol: float = ...,
    atol: float = ...,
) -> tuple[Array, list[Array]]: ...
def thermal_table(
    kind: int,
    params: Sequence[float],
    alpha: Sequence[float],
    beta: Sequence[float],
    ln_s: Sequence[float],
    boundaries: Sequence[float],
    lmax: int,
    energies: Sequence[float],
) -> list[Array]: ...
def infinite_medium(
    total: Sequence[float],
    scatter: Sequence[float],
    fission: Sequence[float],
    tol: float = ...,
) -> tuple[float, Array, int]: ...
def slowing_down_march(
    energies: Array,
    total: Array,
    scatter: Array,
    alphas: Array,
    escape: Array,
    collide: Array,
    volumes: Array,
    history: Array,
    start: int,
    share: float,
    outside: Array,
    kernels: list[Any] = ...,  # MovingKernel | None per nuclide
    tol: float = ...,
    max_passes: int = ...,
) -> tuple[Array, Array]: ...
def free_gas_kernel(
    energies: Array,
    sigma0: Array,
    awr: float,
    temperature: float,
    e_low: float,
    e_top: float,
    dx: float,
    dz: float,
) -> tuple[float, float, float, int, int, Any, Any, Array, Array]: ...
def gauss_legendre(n: int) -> tuple[Array, Array]: ...
def trace_tracks(
    pitch: float,
    nx: int,
    ny: int,
    radii: list[list[float]],
    sectors: list[int],
    offsets: list[int],
    x0: Array,
    y0: Array,
    x1: Array,
    y1: Array,
    gap: float = 0.0,
) -> tuple[npt.NDArray[np.int64], Array, npt.NDArray[np.int64]]: ...
def moc_sweep(
    region: npt.NDArray[np.int64],
    length: Array,
    start: npt.NDArray[np.int64],
    next_fwd: npt.NDArray[np.int64],
    fwd_forward: npt.NDArray[np.int64],
    next_bwd: npt.NDArray[np.int64],
    bwd_forward: npt.NDArray[np.int64],
    weights: Array,
    sin: Array,
    sigma_t: Array,
    source: Array,
    psi_in: Array,
    groups: int,
    faces: tuple[
        npt.NDArray[np.int64], npt.NDArray[np.int64], npt.NDArray[np.int64], Array, int
    ]
    | None = None,
) -> tuple[Array, Array, float, Array]: ...
def cmfd_solve(
    groups: int,
    ny: int,
    nx: int,
    diag: Array,
    east: Array,
    west: Array,
    north: Array,
    south: Array,
    inscatter: Array,
    fission: Array,
    chi: Array,
    phi: Array,
    k: float,
    tol: float = 1e-10,
    max_outer: int = 5000,
    inner: int = 4,
    omega: float = 1.5,
) -> tuple[float, Array, int]: ...
