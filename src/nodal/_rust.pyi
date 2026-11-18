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

def reconstruct(
    model: ReichMoore,
    backgrounds: Sequence[object | None],
    seeds: Sequence[float],
    tol: float = ...,
    atol: float = ...,
    min_rel_width: float = ...,
) -> tuple[Array, Array, Array, Array]: ...
def tab1_eval(table: object, xs: Array) -> Array: ...
def k_constant() -> float: ...
