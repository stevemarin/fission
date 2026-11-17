r"""Nodal Expansion Method (NEM) kernels, vectorized over many nodes/faces.

Along each direction u, the transverse-integrated flux in a node of width h is
expanded on xi = x/h in [-1/2, 1/2]:

    phi(xi) = phibar + a1 f1 + a2 f2 + a3 f3 + a4 f4

    f1 = xi                         f3 = xi (xi - 1/2)(xi + 1/2)
    f2 = 3 xi^2 - 1/4               f4 = (xi^2 - 1/20)(xi^2 - 1/4)

All f_i have zero mean; f3, f4 vanish at the faces, so
    phi(+-1/2) = phibar +- a1/2 + a2/2
    J(+1/2)    = -(D/h) (a1 + 3 a2 + a3/2 + a4/5)
    J(-1/2)    = -(D/h) (a1 - 3 a2 + a3/2 - a4/5)

The 1D equation  -(D/h^2) phi'' + Sigma_r phi = Q(xi) - L(xi)  is weighted with
f1 and f2 (weighted residual), giving the two moment equations

    (Sr/12) a1 - (D/(2h^2) + Sr/120) a3 = S1
    (Sr/20) a2 - (D/(5h^2) + Sr/700) a4 = S2

where S_k = int f_k (Q - L) dxi. Q (fission + in-scatter) is built from the flux
moments  m1 = int f1 phi = a1/12 - a3/120,  m2 = int f2 phi = a2/20 - a4/700  of
the other groups. L is a quadratic transverse leakage fitted to the node and its
two neighbours. The zeroth moment (node balance) is enforced by CMFD.

Nodes are coupled through partial currents built on the heterogeneous surface
flux (f = ADF):  J+- = f phi_s / 4 +- J / 2.  Given the two incoming partial
currents, each node is an independent 4x4 problem (response-matrix form), which
is far better conditioned than imposing surface fluxes on optically thick nodes.
"""

import numpy as np


def _moment_rows(A, b, row, col, D, sr, h, S1, S2):
    dh2 = D / h**2
    A[:, row, col + 0] = sr / 12.0
    A[:, row, col + 2] = -0.5 * dh2 - sr / 120.0
    A[:, row + 1, col + 1] = sr / 20.0
    A[:, row + 1, col + 3] = -0.2 * dh2 - sr / 700.0
    b[:, row] = S1
    b[:, row + 1] = S2


def flux_moments(a):
    """First and second flux moments from expansion coefficients a (M, 4)."""
    return a[:, 0] / 12.0 - a[:, 2] / 120.0, a[:, 1] / 20.0 - a[:, 3] / 700.0


def one_node(D, h, sr, f, phi, S1, S2, jin_m, jin_p):
    """One-node NEM problem with incoming partial currents (all arrays (M,)).

    jin_m enters through the - face (J+ there), jin_p through the + face (J- there).
    Returns a (M, 4) and the outgoing partial currents jout_m, jout_p.
    """
    M = D.size
    A = np.zeros((M, 4, 4))
    b = np.zeros((M, 4))
    _moment_rows(A, b, 0, 0, D, sr, h, S1, S2)
    dh = D / h
    # J+(-1/2) = f phi(-1/2)/4 + J(-1/2)/2
    A[:, 2] = np.stack([-f / 8 - dh / 2, f / 8 + 1.5 * dh, -dh / 4, dh / 10], 1)
    b[:, 2] = jin_m - f * phi / 4
    # J-(+1/2) = f phi(+1/2)/4 - J(+1/2)/2
    A[:, 3] = np.stack([f / 8 + dh / 2, f / 8 + 1.5 * dh, dh / 4, dh / 10], 1)
    b[:, 3] = jin_p - f * phi / 4

    a = np.linalg.solve(A, b[..., None])[..., 0]
    ps_m = f * (phi - 0.5 * a[:, 0] + 0.5 * a[:, 1])
    ps_p = f * (phi + 0.5 * a[:, 0] + 0.5 * a[:, 1])
    return a, 0.5 * ps_m - jin_m, 0.5 * ps_p - jin_p


def leakage_moments(L0, La, Lb, a_lo, a_hi, b_lo, b_hi):
    """Moments int f1 L, int f2 L of the quadratic through three node averages.

    L0 is the node's own transverse-leakage average; La, Lb are the averages of
    two other nodes occupying [a_lo, a_hi] and [b_lo, b_hi] in the node's xi
    coordinate (usually its two neighbours, or two nodes on one side at a
    boundary). All arguments broadcast. The quadratic is L0 + l1 f1 + l2 f2,
    which already preserves L0.
    """
    def f_avg(lo, hi):  # averages of f1, f2 over [lo, hi]
        return 0.5 * (lo + hi), lo * lo + lo * hi + hi * hi - 0.25

    f1a, f2a = f_avg(a_lo, a_hi)
    f1b, f2b = f_avg(b_lo, b_hi)
    ra, rb = La - L0, Lb - L0
    det = f1a * f2b - f2a * f1b
    l1 = (ra * f2b - f2a * rb) / det
    l2 = (f1a * rb - ra * f1b) / det
    return l1 / 12.0, l2 / 20.0
