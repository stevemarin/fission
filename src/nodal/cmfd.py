r"""Coarse-mesh finite-difference (CMFD) global solver.

Partial-current CMFD (pCMFD). Across an interface between nodes L (at -) and
R (at +) the two partial currents are corrected separately:

    J+ = -(Dtilde/2)(phi_R - phi_L) + Dhat+ phi_L
    J- = +(Dtilde/2)(phi_R - phi_L) + Dhat- phi_R
    J  = J+ - J- = (Dtilde + Dhat+) phi_L - (Dtilde + Dhat-) phi_R

Dtilde is the ordinary finite-difference coupling; Dhat+- are the nonlinear
corrections from the nodal solve that make CMFD reproduce the NEM partial
currents. The classic net-current form, J = -Dtilde dphi - Dhat (phi_L + phi_R),
is the special case Dhat+ = -Dhat, Dhat- = +Dhat; it goes unstable once nodes
are several diffusion lengths thick (already the case for the thermal group at
20 cm in a PWR), while pCMFD stays stable. At a boundary face,
J_out = (Dtilde_b + Dhat_b) phi.

The eigenproblem  M phi = (1/k) F phi  is solved by Wielandt-shifted power
iteration; each outer solves  (M - F/k_s) x = b  by ILU-preconditioned BiCGSTAB.
"""

import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import LinearOperator, bicgstab, spilu, splu

from .geometry import Core


class CMFD:
    def __init__(self, core: Core, shift: float = 0.05):
        self.core = c = core
        self.shift = shift
        self.ks = None
        self._ilu, self._ilu_ks, self._ilu_stale, self._ilu_its = None, None, False, None
        self.Dt_i, self.Dp_i, self.Dm_i, self.area_i = [], [], [], []
        self.Dt_b, self.Dh_b, self.area_b = [], [], []
        for d in range(3):
            L, R = c.interfaces[d].L, c.interfaces[d].R
            DL, DR = c.D[L], c.D[R]
            hL, hR = c.h[L, d][:, None], c.h[R, d][:, None]
            self.Dt_i.append(2 * DL * DR / (DL * hR + DR * hL))
            self.Dp_i.append(np.zeros_like(DL))
            self.Dm_i.append(np.zeros_like(DL))
            self.area_i.append(c.vol[L] / c.h[L, d])

            b = c.boundaries[d]
            n = b.node
            kap = 2 * c.D[n] / c.h[n, d][:, None]
            cJ, cphi = b.cJ[:, None], b.cphi[:, None]
            self.Dt_b.append(kap * cphi / (cphi + cJ * kap))
            self.Dh_b.append(np.zeros_like(kap))
            self.area_b.append(c.vol[n] / c.h[n, d])

    # ------------------------------------------------------------------ matrices
    def matrices(self):
        """Loss matrix M and production matrix F (unknown index n*G + g)."""
        c, G, N = self.core, self.core.G, self.core.N
        rows, cols, vals = [], [], []

        def add(r, cc, v):
            rows.append(r.ravel())
            cols.append(cc.ravel())
            vals.append(v.ravel())

        g = np.arange(G)
        node_g = lambda n: n[:, None] * G + g[None, :]

        add(node_g(np.arange(N)), node_g(np.arange(N)), c.sigma_r * c.vol[:, None])
        for gf in range(G):
            for gt in range(G):
                if gf != gt:
                    v = -c.sigma_s[:, gf, gt] * c.vol
                    nz = v != 0
                    n = np.nonzero(nz)[0]
                    add(n * G + gt, n * G + gf, v[nz])
        for d in range(3):
            L, R = c.interfaces[d].L, c.interfaces[d].R
            a = self.area_i[d][:, None]
            cL, cR = self.Dt_i[d] + self.Dp_i[d], self.Dt_i[d] + self.Dm_i[d]
            iL, iR = node_g(L), node_g(R)
            add(iL, iL, a * cL)
            add(iL, iR, -a * cR)
            add(iR, iR, a * cR)
            add(iR, iL, -a * cL)
            n = c.boundaries[d].node
            add(node_g(n), node_g(n), self.area_b[d][:, None] * (self.Dt_b[d] + self.Dh_b[d]))
        M = sp.csc_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                          shape=(N * G, N * G))

        rows, cols, vals = [], [], []
        for gt in range(G):
            for gf in range(G):
                v = c.chi[:, gt] * c.nu_sigma_f[:, gf] * c.vol
                nz = v != 0
                n = np.nonzero(nz)[0]
                add(n * G + gt, n * G + gf, v[nz])
        F = sp.csc_matrix((np.concatenate(vals) if vals else [],
                           (np.concatenate(rows) if rows else [],
                            np.concatenate(cols) if cols else [])),
                          shape=(N * G, N * G))
        return M, F

    # ----------------------------------------------------------------- currents
    def partial_currents(self, phi, d):
        """Interface partial currents (J+, J-) along direction d, each (M, G)."""
        L, R = self.core.interfaces[d].L, self.core.interfaces[d].R
        half = 0.5 * self.Dt_i[d] * (phi[R] - phi[L])
        return -half + self.Dp_i[d] * phi[L], half + self.Dm_i[d] * phi[R]

    def set_partial_correction(self, d, phi, Jplus, Jminus):
        """Choose Dhat+- so CMFD reproduces the given interface partial currents."""
        L, R = self.core.interfaces[d].L, self.core.interfaces[d].R
        half = 0.5 * self.Dt_i[d] * (phi[R] - phi[L])
        self.Dp_i[d][:] = (Jplus + half) / phi[L]
        self.Dm_i[d][:] = (Jminus - half) / phi[R]

    def node_currents(self, phi):
        """Net currents (+u direction) on each node's - and + faces: Jm, Jp (N, 3, G)."""
        c = self.core
        Jm = np.zeros((c.N, 3, c.G))
        Jp = np.zeros((c.N, 3, c.G))
        for d in range(3):
            L, R = c.interfaces[d].L, c.interfaces[d].R
            Jplus, Jminus = self.partial_currents(phi, d)
            J = Jplus - Jminus
            Jp[L, d] = J
            Jm[R, d] = J
            b = c.boundaries[d]
            Jout = (self.Dt_b[d] + self.Dh_b[d]) * phi[b.node]
            plus = b.side > 0
            Jp[b.node[plus], d] = Jout[plus]
            Jm[b.node[~plus], d] = -Jout[~plus]
        return Jm, Jp

    def fd_surface_flux(self, phi):
        """Finite-difference estimate of the surface fluxes, psm/psp (N, 3, G)."""
        c = self.core
        psm = np.zeros((c.N, 3, c.G))
        psp = np.zeros((c.N, 3, c.G))
        for d in range(3):
            L, R = c.interfaces[d].L, c.interfaces[d].R
            wL = c.D[L] / c.h[L, d][:, None]
            wR = c.D[R] / c.h[R, d][:, None]
            ps = (wL * phi[L] + wR * phi[R]) / (wL + wR)
            psp[L, d] = ps
            psm[R, d] = ps
            b = c.boundaries[d]
            n = b.node
            kap = 2 * c.D[n] / c.h[n, d][:, None]
            cJ, cphi = b.cJ[:, None], b.cphi[:, None]
            ps = cJ * kap * phi[n] / (cphi + cJ * kap)
            plus = b.side > 0
            psp[n[plus], d] = ps[plus]
            psm[n[~plus], d] = ps[~plus]
        return psm, psp

    # --------------------------------------------------------------- eigensolve
    def _solver(self, A, ks, rtol):
        """Linear solve with A via ILU-preconditioned BiCGSTAB.

        The ILU factors are kept across calls (the Dhat terms change little from
        one nodal update to the next) and rebuilt only when the shift changes or
        the Krylov iteration count has doubled since the last rebuild.
        """
        if self._ilu is None or self._ilu_ks != ks or self._ilu_stale:
            self._ilu = spilu(A.tocsc(), drop_tol=1e-2, fill_factor=10)
            self._ilu_ks, self._ilu_stale, self._ilu_its = ks, False, None
        P = LinearOperator(A.shape, self._ilu.solve)

        def solve(b, x0):
            count = [0]
            x, info = bicgstab(A, b, x0=x0, M=P, rtol=rtol, atol=0.0, maxiter=500,
                               callback=lambda _: count.__setitem__(0, count[0] + 1))
            if info != 0:  # fall back to a direct solve
                x = splu(A.tocsc()).solve(b)
                self._ilu_stale = True
            elif self._ilu_its is None:
                self._ilu_its = max(count[0], 4)
            elif count[0] > 2 * self._ilu_its:
                self._ilu_stale = True
            return x
        return solve

    def eigensolve(self, phi, k, tol_k=1e-9, tol_s=1e-8, max_outer=2000, warm=True):
        """Wielandt-shifted power iteration. Returns (phi, k, outer iterations).

        If warm is False (no good k estimate yet), unshifted power iteration is
        run until k settles, so the shift lands above the fundamental mode. The
        shift k_s = k + shift is then kept for later calls unless k gets close.
        """
        c = self.core
        M, F = self.matrices()
        x = phi.ravel().copy()
        fs = F @ x
        if not warm:
            self.ks = None
        ks = self.ks
        rtol = min(1e-6, 1e-2 * tol_s)
        solve = self._solver(M - F / ks if ks else M, ks, rtol)
        lam_s = 1.0 / ks if ks else 0.0

        for it in range(1, max_outer + 1):
            lam = 1.0 / k
            x = solve((lam - lam_s) * fs, x)
            fs_new = F @ x
            r = fs_new.sum() / fs.sum()
            k_new = 1.0 / (lam_s + (lam - lam_s) / r)
            x /= r
            fs_new /= r
            err_s = np.abs(fs_new - fs).max() / np.abs(fs_new).max()
            dk = abs(k_new - k)
            k, fs = k_new, fs_new
            if (ks is None and dk < 1e-5 and it > 3) or (ks and k > ks - 0.25 * self.shift):
                # switch on (or move) the Wielandt shift
                self.ks = ks = k + self.shift
                lam_s = 1.0 / ks
                solve = self._solver(M - F / ks, ks, rtol)
            elif ks and dk < tol_k and err_s < tol_s:
                break

        # normalize: total fission production = number of nodes
        x *= c.N / fs.sum()
        return x.reshape(c.N, c.G), k, it
