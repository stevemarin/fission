"""Rust-accelerated versions of the resonance kernels (`nodal._rust`).

Each function takes and returns the same objects as its Python reference, so
the two are interchangeable; tests/test_rust_parity.py checks they agree.

    Python reference                      Rust port (here)
    reich_moore.reich_moore(rr, E)        reich_moore(rr)(E)
    reconstruct.reconstruct(ev, ...)      reconstruct(ev, ...)
"""

from __future__ import annotations

from .. import _rust
from .evaluation import Evaluation
from .reconstruct import BACKGROUND_MT, PointwiseXS, reconstruction_seeds
from .reich_moore import potential_terms, spin_groups
from .resonances import ResolvedRange


def reich_moore(rr: ResolvedRange) -> _rust.ReichMoore:
    """A Rust Reich-Moore evaluator for one resolved range.

    Built from the very same spin groups and potential terms as the Python
    reference; call it with an energy, or use .evaluate(numpy_array).
    """
    return _rust.ReichMoore(spin_groups(rr), potential_terms(rr), rr.spi)


def reconstruct(
    ev: Evaluation,
    lo: float | None = None,
    hi: float | None = None,
    tol: float = 1e-3,
    atol: float = 1e-5,
) -> PointwiseXS:
    """Pointwise resolved-range cross sections, computed in Rust."""
    (iso,) = ev.resonances
    rr = iso.resolved
    if rr is None:
        raise ValueError("no resolved resonance range")
    if rr.lrf != 3:
        raise NotImplementedError("the Rust port covers Reich-Moore (LRF=3) only")
    backgrounds = [
        ev.cross_section(mt).table if mt in ev.reactions else None
        for mt in BACKGROUND_MT.values()
    ]
    seeds = reconstruction_seeds(ev, lo, hi)
    energies, elastic, capture, fission = _rust.reconstruct(
        reich_moore(rr), backgrounds, seeds, tol, atol
    )
    return PointwiseXS(
        energies=energies.tolist(),
        elastic=elastic.tolist(),
        capture=capture.tolist(),
        fission=fission.tolist(),
    )
