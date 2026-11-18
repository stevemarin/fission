"""Type stubs for the Rust extension module (rust/src/lib.rs)."""

from collections.abc import Sequence

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
