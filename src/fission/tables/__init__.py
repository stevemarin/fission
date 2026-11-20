"""Cross-section tables for the nodal code (M3).

The state matrix (``fission.branches``) gives few-group constants at the
nodes of a grid in four variables. A nodal code needs them anywhere inside
it. Two ways are compared here:

**Multilinear interpolation.** In one variable, between nodes x_i and
x_i+1, Sigma(x) = (1 - t) Sigma_i + t Sigma_i+1 with t = (x - x_i) /
(x_i+1 - x_i). In four variables the same is done axis after axis: the
2^4 = 16 corners of the cell around x, each weighted with the product of
its four (1 - t) or t. Outside the grid the end segment is continued
(linear extrapolation). Interpolation is linear in the data, so every
linear relation between constants holds exactly at every point, e.g.
Sigma_t = Sigma_a + sum_g' Sigma_s,g->g'.

The interpolation variable need not be the state variable itself. Doppler
broadening widens a resonance by Gamma_D ~ sqrt(T) (chapter 8), and the
resonance integral grows roughly like sqrt(T): constants are closer to
linear in sqrt(T_fuel) than in T_fuel, so the fuel axis is interpolated in
sqrt(T) (``TRANSFORMS``).

**A fitted polynomial.** Sigma(u) = sum_j c_j m_j(u), with monomials m_j of
the normalized variables u (each axis mapped to [-1, 1]). The terms are
chosen greedily: starting from the constant, add the candidate that most
lowers the leave-one-out error over all the grid points. That error comes
without refitting: for least squares with hat matrix H = X (X^T X)^-1 X^T,
the residual of point i when it is left out is e_i / (1 - H_ii) (the PRESS
residual). One set of terms serves all constants, so the linear relations
hold here too. A fit needs a handful of coefficients per constant instead
of the whole grid, and is smooth between nodes; but it can miss what the
grid resolves (a strong curvature), which is what the check points test.
"""

from __future__ import annotations

import itertools
import json
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from fission.branches import Table
from fission.model import BRANCH_AXES

# the constants of a nodal node (two groups: [G], sigma_s [G, G])
KEYS = ("D", "total", "sigma_a", "sigma_s", "nu_sigma_f", "sigma_f",
        "kappa_sigma_f", "chi", "adf")  # fmt: skip
TRANSFORMS: dict[str, str] = {"fuel_temperature": "sqrt"}
FORMAT = 1


def coordinate(axis: str, value: float, transforms: dict[str, str]) -> float:
    """The interpolation variable of an axis: sqrt(T) or the value itself."""
    kind = transforms.get(axis, "linear")
    if kind == "sqrt":
        return float(np.sqrt(value))
    if kind == "linear":
        return float(value)
    raise ValueError(f"unknown transform {kind!r} for {axis}")


def multilinear(
    nodes: list[np.ndarray], values: np.ndarray, x: list[float]
) -> np.ndarray:
    """Interpolate `values` (shape = node counts + any trailing shape) at x.

    Axis by axis: find the segment [x_i, x_i+1] (the end segment outside the
    grid, which extrapolates), and replace the leading axis by
    (1 - t) v[i] + t v[i + 1]. An axis with a single node is constant."""
    v = np.asarray(values, dtype=float)
    for axis_nodes, xa in zip(nodes, x, strict=True):
        n = len(axis_nodes)
        if n == 1:
            v = v[0]
            continue
        i = int(np.clip(np.searchsorted(axis_nodes, xa) - 1, 0, n - 2))
        t = (xa - axis_nodes[i]) / (axis_nodes[i + 1] - axis_nodes[i])
        v = (1.0 - t) * v[i] + t * v[i + 1]
    return v


def k_infinity(c: dict[str, np.ndarray]) -> float:
    """Few-group k-infinity of a set of constants (no leakage)."""
    a = np.diag(c["total"]) - np.asarray(c["sigma_s"]).T
    return float(np.asarray(c["nu_sigma_f"]) @ np.linalg.solve(a, c["chi"]))


def _values(state) -> list[float]:
    if isinstance(state, dict):
        return [float(state[a]) for a in BRANCH_AXES]
    return [float(getattr(state, a)) for a in BRANCH_AXES]


@dataclass
class XSTable:
    """Few-group constants on a grid, interpolated multilinearly."""

    axes: dict[str, tuple[float, ...]]  # physical node values, every BRANCH_AXES
    data: dict[str, np.ndarray]  # key -> grid shape + the constant's shape
    transforms: dict[str, str] = field(default_factory=lambda: dict(TRANSFORMS))
    meta: dict = field(default_factory=dict)

    @classmethod
    def from_branches(
        cls, table: Table, transforms: dict[str, str] | None = None
    ) -> XSTable:
        keys = [
            k for k in (*KEYS, "k_infinity") if k in next(iter(table.points.values()))
        ]
        meta = {k: table.header[k] for k in ("model", "branches", "solver", "options",
                                             "spectrum", "boundaries_eV", "base")}  # fmt: skip
        return cls(
            {a: tuple(table.axes[a]) for a in BRANCH_AXES},
            {k: table.array(k) for k in keys},
            dict(TRANSFORMS if transforms is None else transforms),
            meta,
        )

    @property
    def shape(self) -> tuple[int, ...]:
        return tuple(len(self.axes[a]) for a in BRANCH_AXES)

    def nodes(self) -> list[np.ndarray]:
        return [np.array([coordinate(a, v, self.transforms) for v in self.axes[a]])
                for a in BRANCH_AXES]  # fmt: skip

    def coordinates(self, state) -> list[float]:
        return [coordinate(a, v, self.transforms)
                for a, v in zip(BRANCH_AXES, _values(state), strict=True)]  # fmt: skip

    def __call__(self, key: str, state) -> np.ndarray:
        """One constant at a state (a State, or a dict of the four variables)."""
        return multilinear(self.nodes(), self.data[key], self.coordinates(state))

    def constants(self, state) -> dict[str, np.ndarray]:
        return {k: self(k, state) for k in KEYS if k in self.data}

    def without(self, axis: str, index: int) -> XSTable:
        """The table with one node of one axis removed (leave one out)."""
        n = BRANCH_AXES.index(axis)
        values = list(self.axes[axis])
        del values[index]
        return XSTable(
            {**self.axes, axis: tuple(values)},
            {k: np.delete(v, index, axis=n) for k, v in self.data.items()},
            dict(self.transforms), dict(self.meta),
        )  # fmt: skip

    def write(self, path: Path) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        data = {
            "format": FORMAT,
            "interpolation": "multilinear",
            "axes": {a: list(v) for a, v in self.axes.items()},
            "axis_order": list(BRANCH_AXES),
            "transforms": self.transforms,
            "meta": self.meta,
            "data": {k: v.tolist() for k, v in self.data.items()},
        }
        path.write_text(json.dumps(data, indent=1))

    @classmethod
    def read(cls, path: Path) -> XSTable:
        d = json.loads(path.read_text())
        if d.get("format") != FORMAT:
            raise ValueError(
                f"{path}: table format {d.get('format')!r}, expected {FORMAT}"
            )
        return cls(
            {a: tuple(v) for a, v in d["axes"].items()},
            {k: np.asarray(v, dtype=float) for k, v in d["data"].items()},
            d["transforms"], d["meta"],
        )  # fmt: skip


# ---------------------------------------------------------------- polynomial fit
Term = tuple[int, ...]  # exponent of each axis


def candidate_terms(table: XSTable) -> list[Term]:
    """Monomials a grid can determine: linear in every axis with two or more
    nodes, squares of axes with three or more, cubes of axes with four or
    more, and the products of two different axes."""
    counts = table.shape
    d = len(counts)
    out: list[Term] = []
    for a in range(d):
        for power in (1, 2, 3):
            if counts[a] > power:
                out.append(tuple(power if j == a else 0 for j in range(d)))
    for a, b in itertools.combinations(range(d), 2):
        if counts[a] > 1 and counts[b] > 1:
            out.append(tuple(1 if j in (a, b) else 0 for j in range(d)))
    return out


@dataclass
class Fit:
    """A polynomial in the normalized variables u in [-1, 1]."""

    terms: list[Term]
    center: np.ndarray  # per axis, in interpolation coordinates
    half: np.ndarray  # half range per axis (1 for a single node)
    transforms: dict[str, str]
    coefficients: dict[str, np.ndarray]  # key -> (terms, *shape)

    def normalized(self, state) -> np.ndarray:
        x = np.array([coordinate(a, v, self.transforms)
                      for a, v in zip(BRANCH_AXES, _values(state), strict=True)])  # fmt: skip
        return (x - self.center) / self.half

    def __call__(self, key: str, state) -> np.ndarray:
        row = _design(self.normalized(state)[None, :], self.terms)[0]
        return np.tensordot(row, self.coefficients[key], axes=1)

    def constants(self, state) -> dict[str, np.ndarray]:
        return {k: self(k, state) for k in KEYS if k in self.coefficients}

    @property
    def size(self) -> int:
        return len(self.terms)


def _design(u: np.ndarray, terms: list[Term]) -> np.ndarray:
    """X[i, j] = prod_a u[i, a] ** terms[j][a]."""
    return np.array([[float(np.prod(row ** np.array(t))) for t in terms] for row in u])


def _grid_points(table: XSTable) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """(normalized coordinates of every node [N, axes], center, half range)."""
    nodes = table.nodes()
    center = np.array([(n.max() + n.min()) / 2 for n in nodes])
    half = np.array([(n.max() - n.min()) / 2 if len(n) > 1 else 1.0 for n in nodes])
    grid = np.array(list(itertools.product(*nodes)))
    return (grid - center) / half, center, half


def _targets(table: XSTable, keys: tuple[str, ...]) -> np.ndarray:
    """Every constant of every node as columns of one matrix [N, M]."""
    n = int(np.prod(table.shape))
    return np.hstack([table.data[k].reshape(n, -1) for k in keys])


def press(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Leave-one-out residuals of the least-squares fit y ~ x (all columns):
    e_i / (1 - H_ii), with H the hat matrix."""
    coef, *_ = np.linalg.lstsq(x, y, rcond=None)
    residual = y - x @ coef
    q, _ = np.linalg.qr(x)
    h = np.sum(q * q, axis=1)  # diagonal of H = Q Q^T
    return residual / (1.0 - h)[:, None]


def _score(x: np.ndarray, y: np.ndarray, scale: np.ndarray) -> float:
    """RMS of the relative leave-one-out errors over every constant."""
    e = press(x, y) / scale
    return float(np.sqrt(np.mean(e * e)))


def select(table: XSTable, keys: tuple[str, ...] = KEYS, max_terms: int = 16,
           ) -> tuple[Fit, list[tuple[Term, float]]]:  # fmt: skip
    """Greedy forward selection of terms by the leave-one-out error.

    Returns the fit with the fewest terms whose score is within 1% of the
    best, and the path [(term added, score)]."""
    keys = tuple(k for k in keys if k in table.data)
    u, _, _ = _grid_points(table)
    y = _targets(table, keys)
    scale = np.maximum(np.abs(y).mean(axis=0), 1e-30)
    keep = np.ptp(y, axis=0) > 1e-12 * scale  # constant columns (adf = 1) fit exactly
    y, scale = y[:, keep], scale[keep]
    d = u.shape[1]
    chosen: list[Term] = [tuple([0] * d)]
    path = [(chosen[0], _score(_design(u, chosen), y, scale))]
    pool = candidate_terms(table)
    while pool and len(chosen) < max_terms and len(chosen) + 1 < len(u):
        scores = [(_score(_design(u, [*chosen, t]), y, scale), t) for t in pool]
        best, term = min(scores)
        chosen.append(term)
        pool.remove(term)
        path.append((term, best))
    # the fewest terms within 1% of the best score: more terms that buy
    # (almost) nothing only fit noise
    history = np.array([s for _, s in path])
    n_best = int(np.argmax(history <= 1.01 * history.min() + 1e-12)) + 1
    return fit(table, chosen[:n_best], keys), path


def fit(table: XSTable, terms: list[Term], keys: tuple[str, ...] = KEYS) -> Fit:
    """Least-squares coefficients of the given terms for every constant."""
    keys = tuple(k for k in keys if k in table.data)
    u, center, half = _grid_points(table)
    x = _design(u, terms)
    n = len(u)
    coefficients = {}
    for k in keys:
        y = table.data[k].reshape(n, -1)
        coef, *_ = np.linalg.lstsq(x, y, rcond=None)
        coefficients[k] = coef.reshape((len(terms),) + table.data[k].shape[4:])
    return Fit(list(terms), center, half, dict(table.transforms), coefficients)


def describe(term: Term) -> str:
    names = ("Tf", "Tm", "rho", "B")
    parts = [n if p == 1 else f"{n}^{p}" for n, p in zip(names, term, strict=True) if p]
    return " ".join(parts) or "1"
