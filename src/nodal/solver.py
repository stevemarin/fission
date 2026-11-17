"""Nonlinear NEM / pCMFD eigenvalue driver.

One nodal iteration:
  1. hand the CMFD solution (node fluxes, k, partial currents) to the nodal solver;
  2. compute transverse-leakage moments from the CMFD face currents;
  3. sweep the one-node NEM problems (partial-current response) over all
     groups and directions at fixed node fluxes and k;
  4. set the pCMFD corrections Dhat+- so CMFD reproduces the NEM partial currents;
  5. re-solve the CMFD eigenproblem.
Repeat until k and the node fluxes stop changing.
"""

import time
from dataclasses import dataclass, field

import numpy as np

from . import nem
from .cmfd import CMFD
from .geometry import Core


@dataclass
class Result:
    k: float
    phi: np.ndarray          # node-average flux (N, G)
    power: np.ndarray        # node power density, volume-weighted fuel mean = 1 (N,)
    nodal_iterations: int
    outer_iterations: int
    converged: bool
    seconds: float
    history: list = field(default_factory=list)  # (k, dk, dphi) per nodal update


class NodalSolver:
    """Steady-state k-eigenvalue solver.

    method='nem':  CMFD-accelerated NEM (production-code style).
    method='fdm':  plain coarse-mesh finite difference (no nodal correction),
                   useful to show how much the nodal method buys at coarse mesh.
    """

    def __init__(self, core: Core, method: str = "nem"):
        if method not in ("nem", "fdm"):
            raise ValueError("method must be 'nem' or 'fdm'")
        self.core = core
        self.method = method
        self.cmfd = CMFD(core)
        N, G = core.N, core.G
        self.a = np.zeros((3, N, G, 4))     # NEM coefficients a1..a4 per direction
        self.m1 = np.zeros((3, N, G))       # flux moments per direction
        self.m2 = np.zeros((3, N, G))
        self.jin = np.zeros((N, 3, 2, G))   # incoming partial currents [node, dir, -/+ face, g]
        self.jout = np.zeros((N, 3, 2, G))  # outgoing partial currents
        self.phi_nodal = None               # node fluxes the last nodal sweep used

    def solve(self, tol_k=1e-7, tol_phi=1e-6, max_nodal=200, sweeps=1,
              verbose=False) -> Result:
        t0 = time.perf_counter()
        c = self.core
        phi, k, n_out = self.cmfd.eigensolve(np.ones((c.N, c.G)), 1.0, warm=False)
        if verbose:
            print(f"  FDM start   k = {k:.6f}  ({n_out} outers)")
        history, it, converged = [], 0, self.method == "fdm"

        if self.method == "nem":
            # boundary faces start from the finite-difference surface flux
            psm, psp = self.cmfd.fd_surface_flux(phi)
            self.jout[:, :, 0] = 0.25 * psm
            self.jout[:, :, 1] = 0.25 * psp
            self.phi_nodal = phi
            dk = dphi = 1.0
            for it in range(1, max_nodal + 1):
                self._nodal_update(phi, k, sweeps)
                # no point converging CMFD far below the current nonlinear change
                phi_new, k_new, n = self.cmfd.eigensolve(
                    phi, k, tol_k=max(0.1 * tol_k, min(1e-5, 0.01 * dk)),
                    tol_s=max(0.1 * tol_phi, min(1e-4, 0.01 * dphi)))
                n_out += n
                dk = abs(k_new - k)
                dphi = np.abs(phi_new - phi).max() / np.abs(phi_new).max()
                phi, k = phi_new, k_new
                history.append((k, dk, dphi))
                if verbose:
                    print(f"  nodal {it:3d}   k = {k:.7f}  dk = {dk:.1e}  dphi = {dphi:.1e}"
                          f"  ({n} outers)")
                if not np.isfinite(k):
                    raise RuntimeError(f"nodal iteration diverged at iteration {it}")
                if dk < tol_k and dphi < tol_phi:
                    converged = True
                    break

        return Result(k=k, phi=phi, power=self.power(phi), nodal_iterations=it,
                      outer_iterations=n_out, converged=converged,
                      seconds=time.perf_counter() - t0, history=history)

    def power(self, phi):
        c = self.core
        p = (c.kappa_sigma_f * phi).sum(axis=1)
        fuel = c.kappa_sigma_f.sum(axis=1) > 0
        return p * c.vol[fuel].sum() / (p * c.vol).sum()

    # ------------------------------------------------------------------ NEM step
    def _transverse_leakage_moments(self, Jm, Jp):
        """L1, L2 (3, N, G): moments of the quadratic transverse leakage per direction."""
        c = self.core
        leak = (Jp - Jm) / c.h[:, :, None]  # net leakage per unit length, (N, 3, G)
        total = leak.sum(axis=1)
        L1 = np.zeros((3, c.N, c.G))
        L2 = np.zeros((3, c.N, c.G))
        for d in c.active_dirs:
            TL = total - leak[:, d]
            h0 = c.h[:, d]
            nm, np_ = c.nbr[:, d, 0], c.nbr[:, d, 1]
            hm = np.where(nm >= 0, c.h[nm, d], h0) / h0
            hp = np.where(np_ >= 0, c.h[np_, d], h0) / h0
            # default stencil: both neighbours (a missing one -> mirror image of the node)
            La = np.where((nm >= 0)[:, None], TL[nm], TL)
            Lb = np.where((np_ >= 0)[:, None], TL[np_], TL)
            a_lo, a_hi = -0.5 - hm, -0.5 + 0 * hm
            b_lo, b_hi = 0.5 + 0 * hp, 0.5 + hp
            # at a non-reflective edge the mirror is wrong (leakage falls off toward
            # a vacuum): use the next node inward instead, a one-sided quadratic
            nmm = np.where(nm >= 0, c.nbr[nm, d, 0], -1)
            npp = np.where(np_ >= 0, c.nbr[np_, d, 1], -1)
            fix_m = (nm < 0) & ~c.reflective[:, d, 0] & (npp >= 0)
            fix_p = (np_ < 0) & ~c.reflective[:, d, 1] & (nmm >= 0)
            hpp = c.h[npp[fix_m], d] / h0[fix_m]
            La[fix_m] = TL[npp[fix_m]]
            a_lo[fix_m], a_hi[fix_m] = 0.5 + hp[fix_m], 0.5 + hp[fix_m] + hpp
            hmm = c.h[nmm[fix_p], d] / h0[fix_p]
            Lb[fix_p] = TL[nmm[fix_p]]
            b_lo[fix_p], b_hi[fix_p] = -0.5 - hm[fix_p] - hmm, -0.5 - hm[fix_p]
            L1[d], L2[d] = nem.leakage_moments(TL, La, Lb, a_lo[:, None], a_hi[:, None],
                                               b_lo[:, None], b_hi[:, None])
        return L1, L2

    def _partial_currents_from_cmfd(self, phi):
        """Hand the CMFD solution to the nodal solver as incoming partial currents.

        Interfaces take the pCMFD partial currents directly. Boundary faces keep
        their last nodal value, rescaled by the change in the node flux.
        """
        c = self.core
        ratio = phi / self.phi_nodal
        for d in c.active_dirs:
            L, R = c.interfaces[d].L, c.interfaces[d].R
            Jplus, Jminus = self.cmfd.partial_currents(phi, d)
            self.jin[R, d, 0] = Jplus
            self.jin[L, d, 1] = Jminus
            b = c.boundaries[d]
            n, s = b.node, (b.side > 0).astype(int)
            self.jin[n, d, s] = b.refl[:, None] * self.jout[n, d, s] * ratio[n]

    def _nodal_update(self, phi, k, sweeps):
        """Run NEM sweeps at fixed node fluxes and k, then refresh the pCMFD corrections."""
        c, cm = self.core, self.cmfd
        self._partial_currents_from_cmfd(phi)
        self.phi_nodal = phi
        Jm, Jp = cm.node_currents(phi)
        L1, L2 = self._transverse_leakage_moments(Jm, Jp)
        for _ in range(sweeps):
            for g in range(c.G):
                # source coupling of group g to every group's flux moments (N, G)
                coef = c.sigma_s[:, :, g].copy()
                coef[:, g] = 0.0
                coef += c.chi[:, g, None] * c.nu_sigma_f / k
                D, sr, f, ph = c.D[:, g], c.sigma_r[:, g], c.adf[:, g], phi[:, g]
                for d in c.active_dirs:
                    S1 = (coef * self.m1[d]).sum(axis=1) - L1[d, :, g]
                    S2 = (coef * self.m2[d]).sum(axis=1) - L2[d, :, g]
                    a, jo_m, jo_p = nem.one_node(D, c.h[:, d], sr, f, ph, S1, S2,
                                                 self.jin[:, d, 0, g], self.jin[:, d, 1, g])
                    self.a[d, :, g] = a
                    self.m1[d, :, g], self.m2[d, :, g] = nem.flux_moments(a)
                    self.jout[:, d, 0, g] = jo_m
                    self.jout[:, d, 1, g] = jo_p
                    # exchange: what leaves one node enters its neighbour
                    L, R = c.interfaces[d].L, c.interfaces[d].R
                    self.jin[R, d, 0, g] = jo_p[L]
                    self.jin[L, d, 1, g] = jo_m[R]
                    b = c.boundaries[d]
                    s = (b.side > 0).astype(int)
                    self.jin[b.node, d, s, g] = b.refl * self.jout[b.node, d, s, g]

        # nonlinear coupling corrections from the nodal partial currents
        for d in c.active_dirs:
            L, R = c.interfaces[d].L, c.interfaces[d].R
            cm.set_partial_correction(d, phi, self.jout[L, d, 1], self.jout[R, d, 0])
            b = c.boundaries[d]
            n, s = b.node, (b.side > 0).astype(int)
            Jout = self.jout[n, d, s] - self.jin[n, d, s]
            cm.Dh_b[d][:] = Jout / phi[n] - cm.Dt_b[d]
