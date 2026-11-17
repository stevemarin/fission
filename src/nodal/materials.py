"""Homogenized multigroup cross-section sets."""

from dataclasses import dataclass, field

import numpy as np


@dataclass
class Material:
    """Few-group homogenized cross sections for one node composition.

    sigma_s[g, g'] is scattering FROM group g TO group g'. The within-group
    term (diagonal) is ignored, since it cancels out of the diffusion equation.
    """

    D: np.ndarray
    sigma_a: np.ndarray
    nu_sigma_f: np.ndarray
    sigma_s: np.ndarray
    chi: np.ndarray | None = None
    kappa_sigma_f: np.ndarray | None = None
    adf: np.ndarray | None = None  # assembly discontinuity factors (one per group)
    name: str = ""

    G: int = field(init=False)

    def __post_init__(self):
        self.D = np.atleast_1d(np.asarray(self.D, float))
        self.G = G = self.D.size
        self.sigma_a = np.atleast_1d(np.asarray(self.sigma_a, float))
        self.nu_sigma_f = np.atleast_1d(np.asarray(self.nu_sigma_f, float))
        self.sigma_s = np.asarray(self.sigma_s, float).reshape(G, G)
        if self.chi is None:
            self.chi = np.eye(G)[0]
        self.chi = np.asarray(self.chi, float)
        if self.kappa_sigma_f is None:
            self.kappa_sigma_f = self.nu_sigma_f.copy()
        self.kappa_sigma_f = np.asarray(self.kappa_sigma_f, float)
        if self.adf is None:
            self.adf = np.ones(G)
        self.adf = np.asarray(self.adf, float)
        for name in ("sigma_a", "nu_sigma_f", "chi", "kappa_sigma_f", "adf"):
            if getattr(self, name).shape != (G,):
                raise ValueError(f"{name} must have {G} entries")

    @property
    def sigma_r(self) -> np.ndarray:
        """Removal cross section: absorption plus out-scatter (self-scatter excluded)."""
        out = self.sigma_s.sum(axis=1) - np.diag(self.sigma_s)
        return self.sigma_a + out

    @classmethod
    def two_group(cls, D1, D2, sa1, sa2, s12, nsf1, nsf2, ksf1=None, ksf2=None,
                  adf=None, name=""):
        ksf = None if ksf1 is None else [ksf1, ksf2]
        return cls(D=[D1, D2], sigma_a=[sa1, sa2], nu_sigma_f=[nsf1, nsf2],
                   sigma_s=[[0.0, s12], [0.0, 0.0]], chi=[1.0, 0.0],
                   kappa_sigma_f=ksf, adf=adf, name=name)
