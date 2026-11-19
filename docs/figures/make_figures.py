"""Generate the figures used in docs/theory/*.md.

    uv run python docs/figures/make_figures.py        (or: make figures)

Every figure is computed with the project's own code from ENDF/B-VIII.1 data,
so they double as visual checks of the implementation. Output: SVG files next
to this script. Dense curves are rasterized inside the SVG to keep files small.
"""

from __future__ import annotations

import math
import os
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from fission.data import EndfLibrary
from fission.endf import Evaluation, fast
from fission.endf.cache import cached_pointwise
from fission.endf.channel import penetrability_shift_phase
from fission.endf.distributions.angular import angular_distribution
from fission.endf.distributions.energy import energy_distribution
from fission.endf.distributions.kinematics import kalbach_mann
from fission.endf.distributions.mf6 import energy_angle
from fission.endf.distributions.nu import (
    delayed_groups,
    multiplicity_distribution,
    nu_bar,
)
from fission.endf.distributions.sab import scattering_law
from fission.endf.interpolation import interpolate
from fission.endf.pointwise import PointwiseTable
from fission.endf.reconstruct import adaptive_grid, reconstruct
from fission.endf.slbw import potential_scattering, slbw
from fission.processing.cache import cached_doppler
from fission.processing.doppler import (
    broaden_table,
    broadened_at,
    doppler_width,
)
from fission.processing.groups import XMAS_172, equal_lethargy
from fission.processing.library import cached_groups
from fission.processing.thermal import FreeGas, Tabulated, thermal_kernel
from fission.processing.weights import ThermalSlowingDownFission

OUT = Path(__file__).parent
PREVIEW = (
    Path(os.environ["FIG_PREVIEW_DIR"]) if os.environ.get("FIG_PREVIEW_DIR") else None
)
LIB = EndfLibrary()

plt.rcParams.update(
    {
        "figure.figsize": (7.0, 4.2),
        "figure.dpi": 110,
        "savefig.dpi": 150,
        "font.size": 10,
        "axes.grid": True,
        "grid.alpha": 0.3,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "legend.frameon": False,
        "svg.fonttype": "none",
        "svg.hashsalt": "fission",  # stable element ids, so reruns give clean diffs
    }
)
C = ["#1f5fa8", "#d1495b", "#2e8b57", "#e09f3e", "#6a4c93", "#555555"]


def save(fig, name):
    fig.tight_layout()
    fig.savefig(OUT / f"{name}.svg", metadata={"Date": None})
    if PREVIEW:  # optional PNG copies for a quick look
        fig.savefig(PREVIEW / f"{name}.png")
    plt.close(fig)
    print("wrote", name)


def evaluation(nuclide):
    return Evaluation.from_file(LIB.neutron(nuclide))


def logspace(lo, hi, n):
    return np.exp(np.linspace(math.log(lo), math.log(hi), n))


# ---------------------------------------------------------------- chapter 1
def fig_h1_capture():
    table = cached_pointwise("H1")
    e = np.array(table.energies)
    fig, ax = plt.subplots()
    ax.loglog(e, table.xs[2], color=C[0], label="elastic (MT2)")
    ax.loglog(
        e, table.xs[102], color=C[1], label=r"capture (MT102), $^1$H(n,$\gamma$)$^2$H"
    )
    ref = table(102, 0.0253) * np.sqrt(0.0253 / e)
    ax.loglog(
        e,
        ref,
        "--",
        color=C[5],
        lw=0.8,
        label=r"$1/v$: $\sigma_\gamma(E_0)\sqrt{E_0/E}$",
    )
    ax.axvline(0.0253, color=C[5], lw=0.6, ls=":")
    ax.annotate(
        "0.0253 eV\n(2200 m/s)",
        (0.0253, 3e-6),
        xytext=(4, 0),
        textcoords="offset points",
        fontsize=8,
        color=C[5],
    )
    ax.set(
        xlabel="neutron energy E (eV)",
        ylabel="cross section (b)",
        ylim=(1e-6, 1e2),
        title="H-1: smooth cross sections, capture follows the 1/v law",
    )
    ax.legend(loc="lower left")
    save(fig, "h1_capture")


def fig_u238_overview():
    table = cached_pointwise("U238")
    ev = evaluation("U238")
    iso = ev.resonances[0]
    e = np.array(table.energies)
    fig, ax = plt.subplots(figsize=(7.5, 4.4))
    ax.loglog(e, table.xs[1], color=C[0], lw=0.4, label="total", rasterized=True)
    ax.loglog(e, table.xs[102], color=C[1], lw=0.4, label="capture", rasterized=True)
    ax.loglog(e, table.xs[18], color=C[2], lw=0.6, label="fission", rasterized=True)
    for rng, label, color in (
        (iso.resolved, "resolved", "#d9e6f5"),
        (iso.unresolved, "unresolved", "#f5e6d9"),
    ):
        ax.axvspan(max(1e-05, rng.el), rng.eh, color=color, alpha=0.6, lw=0)
        ax.text(
            math.sqrt(max(rng.el, 1e-5) * rng.eh) if label == "unresolved" else 1.0,
            3e4,
            label + " range",
            fontsize=8,
            ha="center",
        )
    ax.set(
        xlabel="neutron energy E (eV)",
        ylabel="cross section (b)",
        ylim=(1e-4, 1e5),
        xlim=(1e-5, 3e7),
        title="U-238 (0 K), ENDF/B-VIII.1, reconstructed",
    )
    ax.legend(loc="lower left", ncols=3)
    save(fig, "u238_overview")


# ---------------------------------------------------------------- chapter 2
def fig_interpolation_laws():
    x1, x2, y1, y2 = 1.0, 10.0, 8.0, 0.5
    xs = np.linspace(x1, x2, 300)
    names = {
        1: "1 histogram",
        2: "2 lin-lin",
        3: "3 lin-log (y lin in ln x)",
        4: "4 log-lin (ln y lin in x)",
        5: "5 log-log (power law)",
    }
    fig, ax = plt.subplots(figsize=(6.5, 4.0))
    for law, color in zip(names, C, strict=False):
        ax.plot(
            xs,
            [interpolate(law, x, x1, x2, y1, y2) for x in xs],
            color=color,
            label=names[law],
        )
    ax.plot([x1, x2], [y1, y2], "ko", ms=5, zorder=5)
    ax.set(
        xlabel="x",
        ylabel="y",
        title="The five ENDF interpolation laws between the same two points",
    )
    ax.legend()
    save(fig, "interpolation_laws")


# ---------------------------------------------------------------- chapter 3
def fig_hard_sphere():
    rho = logspace(1e-3, 20, 400)
    fig, (a, b) = plt.subplots(1, 2, figsize=(8.0, 3.8))
    for l in range(4):
        p = [penetrability_shift_phase(l, r, r)[0] for r in rho]
        phi = [penetrability_shift_phase(l, r, r)[2] for r in rho]
        a.loglog(rho, p, color=C[l], label=f"$l={l}$")
        b.semilogx(rho, np.sin(phi) ** 2, color=C[l], label=f"$l={l}$")
    a.set(
        xlabel=r"$\rho = ka$",
        ylabel=r"$P_l(\rho)$",
        ylim=(1e-12, 30),
        title=r"penetrability: $P_l \sim \rho^{2l+1}$ for small $\rho$",
    )
    b.set(
        xlabel=r"$\hat\rho = k\,AP$",
        ylabel=r"$\sin^2\varphi_l$",
        title="hard-sphere phase shift",
    )
    a.legend()
    rho_u = 2.1968e-3 * 236 / 237 * math.sqrt(1e3) * 0.840  # k a, U-238 channel radius
    for ax in (a, b):
        ax.axvline(rho_u, color=C[5], lw=0.7, ls=":")
    a.annotate(
        "U-238 at 1 keV",
        (rho_u, 1e-9),
        xytext=(4, 0),
        textcoords="offset points",
        fontsize=8,
    )
    save(fig, "hard_sphere")


# ---------------------------------------------------------------- chapter 4
def u238_667():
    rr = evaluation("U238").resonances[0].resolved
    group = rr.groups[0]
    res = min(group.resonances, key=lambda r: abs(r.energy - 6.674))
    return rr, group, res


def fig_slbw_667():
    rr, group, res = u238_667()
    e = np.linspace(6.3, 7.05, 1500)
    cap = [slbw(rr, x, only=[res]).capture for x in e]
    el = [slbw(rr, x, only=[res]).elastic for x in e]
    pot = potential_scattering(rr, group, 6.674)
    fig, ax = plt.subplots()
    ax.semilogy(e, cap, color=C[1], label="capture")
    ax.semilogy(e, el, color=C[0], label="elastic")
    ax.axhline(
        pot, color=C[5], ls="--", lw=0.8, label=r"potential scattering $4\pi a^2$"
    )
    i = int(np.argmin(el))
    ax.annotate(
        "interference dip",
        (e[i], el[i]),
        xytext=(-90, 10),
        textcoords="offset points",
        arrowprops={"arrowstyle": "->", "lw": 0.7},
        fontsize=8,
    )
    ax.set(
        xlabel="E (eV)",
        ylabel="cross section (b)",
        ylim=(0.05, 5e4),
        title="U-238, the 6.67 eV resonance alone (SLBW, 0 K)",
    )
    ax.legend(loc="upper left")
    save(fig, "slbw_667")


# ---------------------------------------------------------------- chapter 5
def fig_mlbw_zr90():
    rr = evaluation("Zr90").resonances[0].resolved
    model = fast.mlbw(rr)
    e = np.linspace(12_900, 13_800, 2000)
    fig, ax = plt.subplots()
    ax.plot(
        e,
        [slbw(rr, x).elastic for x in e],
        color=C[1],
        label="SLBW (levels independent)",
    )
    ax.plot(e, [model(x)[0] for x in e], color=C[0], label="MLBW (levels interfere)")
    ax.axhline(0, color="k", lw=0.6)
    ax.set(
        xlabel="E (eV)",
        ylabel="elastic cross section (b)",
        ylim=(-15, 80),
        title="Zr-90, two s-wave resonances at 13.37 and 13.44 keV",
    )
    ax.legend(loc="upper left")
    save(fig, "mlbw_zr90")


# ---------------------------------------------------------------- chapter 6
def fig_rm_vs_slbw():
    u238 = evaluation("U238").resonances[0].resolved
    u235 = evaluation("U235").resonances[0].resolved
    rm238, rm235 = fast.reich_moore(u238), fast.reich_moore(u235)
    fig, (a, b) = plt.subplots(1, 2, figsize=(8.6, 3.8))
    e = np.linspace(6.35, 6.72, 800)
    a.plot(e, [slbw(u238, x).elastic for x in e], color=C[1], label="SLBW")
    a.plot(e, [rm238(x)[0] for x in e], color=C[0], label="Reich-Moore")
    a.axhline(0, color="k", lw=0.6)
    a.set(
        xlabel="E (eV)",
        ylabel="elastic (b)",
        ylim=(-4, 30),
        title="U-238 elastic below 6.67 eV",
    )
    a.legend(loc="upper left")
    e = logspace(0.01, 5.0, 1500)
    b.loglog(e, [slbw(u235, x).fission for x in e], color=C[1], label="SLBW")
    b.loglog(e, [rm235(x)[2] for x in e], color=C[0], label="Reich-Moore")
    b.set(xlabel="E (eV)", ylabel="fission (b)", title="U-235 fission")
    b.legend(loc="lower left")
    save(fig, "rm_vs_slbw")


# ---------------------------------------------------------------- chapter 7
def fig_adaptive_grid():
    ev = evaluation("U238")
    table = reconstruct(ev, 6.45, 6.9)
    e = np.array(table.energies)
    fig, (a, b) = plt.subplots(
        2, 1, figsize=(7.0, 5.2), sharex=True, gridspec_kw={"height_ratios": [3, 1.3]}
    )
    a.semilogy(e, table.capture, "-", color=C[1], lw=0.8)
    a.semilogy(e, table.capture, ".", color=C[1], ms=2.5)
    a.set(
        ylabel="capture (b)",
        title=f"Adaptive grid on the 6.67 eV resonance: {len(e)} points, tol = 0.1%",
    )
    mid = 0.5 * (e[1:] + e[:-1])
    b.semilogy(mid, np.diff(e), ".", color=C[0], ms=2.5)
    b.set(xlabel="E (eV)", ylabel=r"spacing $\Delta E$ (eV)")
    save(fig, "adaptive_grid")


def fig_points_per_decade():
    table = cached_pointwise("U238")
    e = np.array(table.energies)
    edges = 10.0 ** np.arange(-5, 8)
    counts, _ = np.histogram(e, bins=edges)
    fig, ax = plt.subplots(figsize=(6.5, 3.4))
    ax.stairs(counts, edges, fill=True, color=C[0], alpha=0.85)
    for left, n in zip(edges[:-1], counts, strict=True):
        ax.text(left * 3.16, n * 1.15, f"{n:,}", ha="center", fontsize=7)
    ax.set(
        xscale="log",
        yscale="log",
        xlabel="energy (eV)",
        ylabel="points per decade",
        ylim=(10, 1e6),
        title=f"U-238 table: {len(e):,} points, where they go",
    )
    save(fig, "points_per_decade")


# ---------------------------------------------------------------- chapter 8
TEMPS = (293.6, 600.0, 900.0, 1200.0)


def _window(table, mt, lo, hi):
    e = np.array(table.energies)
    k = (e >= lo) & (e <= hi)
    return e[k], np.array(table.xs[mt])[k]


def fig_doppler_667():
    fig, ax = plt.subplots(figsize=(7.0, 4.0))
    curves = [("0 K", cached_pointwise("U238"))] + [
        (f"{t:g} K", cached_doppler("U238", t)) for t in TEMPS
    ]
    for (label, table), color in zip(curves, C, strict=False):
        e, s = _window(table, 102, 6.2, 7.15)
        ax.semilogy(e, s, color=color, lw=1.2, label=label)
    ax.set(
        xlabel="E (eV)",
        ylabel="capture (b)",
        ylim=(3, 4e4),
        title="U-238 6.67 eV resonance: Doppler broadening (free gas)",
    )
    ax.legend(loc="upper right")
    save(fig, "doppler_667")


def fig_thermal_rise():
    """H-1 elastic for a free-gas target: the 1/v rise at low energy."""
    fig, ax = plt.subplots(figsize=(6.5, 3.8))
    e = np.array(logspace(2e-5, 10.0, 400))
    cold = cached_pointwise("H1")
    ax.loglog(e, [cold(2, x) for x in e], color=C[0], label="0 K")
    for t, color in zip((293.6, 900.0), C[1:], strict=False):
        hot = cached_doppler("H1", t)
        ax.loglog(e, [hot(2, x) for x in e], color=color, label=f"{t:g} K")
    ax.set(
        xlabel="E (eV)",
        ylabel="elastic (b)",
        title="H-1 elastic, free-gas target: the low-energy rise",
    )
    ax.legend()
    save(fig, "thermal_rise")


def fig_corner_rounding():
    """A ramp between two plateaus: why extra seeds are needed."""
    cold = PointwiseTable([90.0, 100.0, 110.0, 120.0], {18: [50.0, 50.0, 10.0, 10.0]})
    awr, temp = (
        233.0,
        293.6,
    )  # Gamma_D(100 eV) = 0.21 eV, much less than the 10 eV interval

    def f(v):
        return tuple(broadened_at(cold, awr, temp, [v])[18])

    naive_e, naive_v = adaptive_grid(f, cold.energies, 1e-3, 1e-5)
    good = broaden_table(cold, awr, temp)
    fine = np.linspace(95.0, 115.0, 4001)
    exact = np.array(broadened_at(cold, awr, temp, fine)[18])
    naive = np.interp(fine, naive_e, [v[0] for v in naive_v])
    fig, (l, r) = plt.subplots(1, 2, figsize=(8.6, 3.6))
    l.plot(fine, [cold(18, v) for v in fine], color="0.6", lw=1, ls="--", label="0 K")
    l.plot(fine, exact, color=C[0], lw=1.6, label="broadened (exact)")
    l.plot(
        naive_e,
        [v[0] for v in naive_v],
        "o",
        color=C[1],
        ms=4,
        label="midpoint test only",
    )
    l.set(
        xlim=(95, 115),
        xlabel="E (eV)",
        ylabel=r"$\sigma$ (b)",
        title=f"ramp, $\\Gamma_D$ = {doppler_width(100, awr, temp):.2f} eV",
    )
    l.legend(fontsize=8)
    k = (fine > 99.0) & (fine < 101.5)
    r.plot(
        fine[k], 100 * (naive[k] / exact[k] - 1), color=C[1], label="midpoint test only"
    )
    good_k = np.interp(fine[k], good.energies, good.xs[18])
    r.plot(
        fine[k], 100 * (good_k / exact[k] - 1), color=C[2], label="with Doppler seeds"
    )
    r.axhspan(-0.1, 0.1, color="0.85", label="tolerance 0.1%")
    r.set(
        xlabel="E (eV)",
        ylabel="table error (%)",
        title="error near the corner at 100 eV",
    )
    r.legend(fontsize=8)
    save(fig, "corner_rounding")


def fig_shielded_ri():
    def shielded(t, s0):
        e, g = _window(t, 102, 0.5, 2e4)
        _, tot = _window(t, 1, 0.5, 2e4)
        w = 1.0 / e if s0 is None else s0 / ((tot + s0) * e)
        return np.trapezoid(g * w, e)

    tables = [(0.0, cached_pointwise("U238"))] + [
        (t, cached_doppler("U238", t)) for t in TEMPS
    ]
    fig, ax = plt.subplots(figsize=(6.5, 3.8))
    for s0, color, label in [
        (None, C[5], r"dilute ($\sigma_0 \to \infty$)"),
        (1000.0, C[0], r"$\sigma_0$ = 1000 b"),
        (100.0, C[2], r"$\sigma_0$ = 100 b"),
        (50.0, C[1], r"$\sigma_0$ = 50 b"),
    ]:
        ri = np.array([shielded(t, s0) for _, t in tables])
        ax.plot([t for t, _ in tables], ri / ri[1], "o-", color=color, label=label)
    ax.set(
        xlabel="temperature (K)",
        ylabel="capture RI / RI(293.6 K)",
        title="U-238 capture resonance integral vs. temperature",
    )
    ax.legend()
    save(fig, "shielded_ri")


# ---------------------------------------------------------------- chapter 9
def fig_nu_multiplicity():
    u235, u238 = evaluation("U235"), evaluation("U238")
    fig, (a, b) = plt.subplots(1, 2, figsize=(8.6, 3.6))
    e = np.array(logspace(2e-5, 2e7, 600))
    for ev, name, color in ((u235, "U-235", C[0]), (u238, "U-238", C[1])):
        a.semilogx(
            e, [nu_bar(ev, 452)(x) for x in e], color=color, label=f"{name} total"
        )
        a.semilogx(
            e,
            [100 * nu_bar(ev, 455)(x) for x in e],
            color=color,
            ls="--",
            label=f"{name} delayed ×100",
        )
    a.set(
        xlabel="incident energy (eV)", ylabel=r"$\bar\nu$", title="neutrons per fission"
    )
    a.legend(fontsize=8)
    width = 0.38
    for k, (ev, name, energy, color) in enumerate(
        ((u235, "U-235, thermal", 0.0253, C[0]), (u238, "U-238, 2 MeV", 2e6, C[1]))
    ):
        md = multiplicity_distribution(ev)
        p = md.probabilities(energy)[:9]
        b.bar(np.arange(len(p)) + (k - 0.5) * width, p, width, color=color, label=name)
    b.set(
        xlabel=r"prompt neutrons $\nu$",
        ylabel=r"$P(\nu)$",
        title="multiplicity distribution",
    )
    b.legend(fontsize=8)
    save(fig, "nu_multiplicity")


def fig_fission_spectrum():
    ev = evaluation("U235")
    chi = energy_distribution(ev, 18)
    groups = delayed_groups(ev)
    e_out = np.array(logspace(1e3, 2e7, 500))
    fig, ax = plt.subplots(figsize=(6.8, 4.0))
    for energy, color, label in (
        (0.0253, C[0], "prompt, thermal fission"),
        (1.4e7, C[1], "prompt, 14 MeV fission"),
    ):
        ax.loglog(
            e_out / 1e6, [1e6 * chi(energy, x) for x in e_out], color=color, label=label
        )
    delayed = [1e6 * groups.spectra(0.0253, x) for x in e_out]
    ax.loglog(e_out / 1e6, delayed, color=C[2], label="delayed (6 groups)")
    ax.set(
        xlabel="outgoing energy (MeV)",
        ylabel=r"$\chi(E')$ (1/MeV)",
        ylim=(1e-4, 3),
        title="U-235 fission neutron spectra",
    )
    ax.legend()
    save(fig, "fission_spectrum")


def fig_elastic_angular():
    fig, axes = plt.subplots(1, 2, figsize=(8.6, 3.6), sharey=False)
    mu = np.linspace(-1, 1, 401)
    for ax, nuc in zip(axes, ("O16", "U238"), strict=True):
        ad = angular_distribution(evaluation(nuc), 2)
        for energy, color in zip((1e4, 1e6, 3e6, 1e7), C, strict=False):
            ax.semilogy(
                mu,
                [ad(energy, m) for m in mu],
                color=color,
                label=f"{energy / 1e6:g} MeV",
            )
        ax.set(xlabel=r"$\mu_{\mathrm{cm}}$", ylabel=r"$f(\mu)$")
        ax.legend(fontsize=8)
    axes[0].set_title("O-16 elastic (CM)")
    axes[1].set_title("U-238 elastic (CM)")
    save(fig, "elastic_angular")


def fig_continuum():
    ev = evaluation("U238")
    c = energy_angle(ev, 91).products[0].distribution
    fig, (a, b) = plt.subplots(1, 2, figsize=(8.6, 3.6))
    for energy, color in zip((2e6, 5e6, 1.4e7), C, strict=False):
        e_out = np.linspace(0, energy, 600)
        a.plot(
            e_out / 1e6,
            [1e6 * c.energy_pdf(energy, x) for x in e_out],
            color=color,
            label=f"E = {energy / 1e6:g} MeV",
        )
    a.set(
        xlabel="outgoing energy (MeV)",
        ylabel="f(E') (1/MeV)",
        xlim=(0, 14),
        title="U-238 continuum inelastic (MT91)",
    )
    a.legend(fontsize=8)
    mu = np.linspace(-1, 1, 201)
    for (aa, r), color in zip(
        ((0.3, 0.0), (0.3, 1.0), (2.0, 0.0), (2.0, 1.0)), C, strict=False
    ):
        b.plot(
            mu,
            [kalbach_mann(m, aa, r) for m in mu],
            color=color,
            label=f"a = {aa}, r = {r}",
        )
    b.set(
        xlabel=r"$\mu_{\mathrm{cm}}$",
        ylabel=r"$f(\mu)$",
        title="Kalbach-Mann angular shapes",
    )
    b.legend(fontsize=8)
    save(fig, "continuum")


# ---------------------------------------------------------------- chapter 10
def fig_weight_function():
    w = ThermalSlowingDownFission.at(293.6)
    e = np.array(logspace(1e-5, 2e7, 2000))
    fig, ax = plt.subplots(figsize=(6.8, 3.8))
    ax.loglog(e, [x * w(x) for x in e], color=C[0])
    for x, label in ((w.eb, "$E_b = 2kT$"), (w.ec, "$E_c = 1.5\\,T_f$")):
        ax.axvline(x, color="0.6", ls="--", lw=0.8)
        ax.text(x * 1.3, 2e-4, label, fontsize=9)
    ax.set(
        xlabel="E (eV)",
        ylabel=r"$E\,w(E)$ (flux per unit lethargy)",
        ylim=(1e-5, 3),
        title="Weighting spectrum: Maxwellian - 1/E - fission",
    )
    save(fig, "weight_function")


def fig_group_capture():
    table = cached_doppler("U238", 293.6)
    gl = cached_groups("U238", 293.6)
    e, s = _window(table, 102, 1.0, 1e3)
    fig, ax = plt.subplots(figsize=(7.0, 4.0))
    ax.loglog(e, s, color="0.65", lw=0.6, label="pointwise (293.6 K)", rasterized=True)
    edges = np.array(XMAS_172.boundaries[::-1])
    values = gl.sigma[102][::-1]
    k = (edges[:-1] >= 1.0) & (edges[1:] <= 1e3)
    ax.stairs(
        values[k],
        edges[np.r_[k, False] | np.r_[False, k]],
        color=C[1],
        lw=1.4,
        label="XMAS-172 group average",
    )
    ax.set(
        xlabel="E (eV)",
        ylabel="capture (b)",
        xlim=(1, 1e3),
        ylim=(1e-1, 4e4),
        title="U-238 capture: pointwise and infinitely dilute group values",
    )
    ax.legend(loc="lower left")
    save(fig, "group_capture")


def fig_scatter_matrices():
    fig, axes = plt.subplots(1, 2, figsize=(8.6, 4.0))
    for ax, nuc in zip(axes, ("H1", "U238"), strict=True):
        gl = cached_groups(nuc, 293.6)
        m = gl.scatter_elastic[0]
        im = ax.imshow(
            np.where(m > 0, m, np.nan),
            norm=matplotlib.colors.LogNorm(vmin=1e-6, vmax=m.max()),
            cmap="viridis",
            interpolation="nearest",
        )
        ax.set(
            xlabel="to group g'",
            ylabel="from group g",
            title=f"{nuc[0]}-{nuc[1:]} elastic $P_0$ (b)",
        )
        fig.colorbar(im, ax=ax, shrink=0.8)
    save(fig, "scatter_matrices")


def fig_chi_groups():
    gl = cached_groups("U235", 293.6)
    chi = energy_distribution(evaluation("U235"), 18)
    du = np.array([XMAS_172.lethargy_width(g) for g in range(XMAS_172.size)])
    edges = np.array(XMAS_172.boundaries[::-1])
    fig, ax = plt.subplots(figsize=(6.8, 3.8))
    e = np.array(logspace(1e3, 2e7, 400))
    ax.semilogx(
        e / 1e6,
        [x * chi(0.0253, x) for x in e],
        color="0.6",
        label="pointwise $E\\,\\chi(E)$",
    )
    ax.stairs(
        (gl.chi / du)[::-1],
        edges / 1e6,
        color=C[0],
        lw=1.4,
        label=r"$\chi_g / \Delta u_g$",
    )
    ax.set(
        xlabel="E (MeV)",
        ylabel="per unit lethargy",
        xlim=(1e-3, 20),
        title="U-235 fission spectrum in XMAS-172 groups",
    )
    ax.legend()
    save(fig, "chi_groups")


# ---------------------------------------------------------------- chapter 11
def _water(temperature):
    return Tabulated(
        scattering_law(Evaluation.from_file(LIB.thermal("HinH2O")), temperature)
    )


def fig_thermal_xs():
    e = np.array(logspace(1e-5, 9.0, 70))
    fig, ax = plt.subplots(figsize=(6.8, 4.0))
    for model, color, label in (
        (_water(293.6), C[0], "H in H$_2$O, 293.6 K"),
        (_water(600.0), C[2], "H in H$_2$O, 600 K"),
    ):
        ax.loglog(
            e,
            [thermal_kernel(model, XMAS_172, 0, x)[1] for x in e],
            color=color,
            label=label,
        )
    free = _water(293.6).law.sigma_free
    gas = FreeGas(0.99917, 293.6)
    ax.loglog(
        e,
        [free * thermal_kernel(gas, XMAS_172, 0, x)[1] for x in e],
        color=C[1],
        ls="--",
        label="free H gas, 293.6 K",
    )
    ax.axhline(free, color="0.6", lw=0.8, ls=":")
    ax.text(2e-5, free * 0.82, "free atom 20.4 b", fontsize=8, color="0.4")
    ax.set(
        xlabel="E (eV)",
        ylabel="scattering cross section per H (b)",
        title="Hydrogen in water: bound vs. free",
    )
    ax.legend()
    save(fig, "thermal_xs")


def fig_thermal_spectra():
    fine = equal_lethargy(500, 1e-4, 3.0)
    widths = np.array([fine.bounds(g)[1] - fine.bounds(g)[0] for g in range(fine.size)])
    mids = np.array(
        [np.sqrt(fine.bounds(g)[0] * fine.bounds(g)[1]) for g in range(fine.size)]
    )
    fig, axes = plt.subplots(1, 2, figsize=(8.6, 3.7))
    water, free = _water(293.6), FreeGas(0.99917, 293.6)
    scale = (
        water.law.sigma_free
    )  # the FreeGas kernel is per 1 b of free-atom cross section
    for ax, energy in zip(axes, (0.0253, 0.5), strict=True):
        kw, _ = thermal_kernel(water, fine, 0, energy)
        kf, _ = thermal_kernel(free, fine, 0, energy)
        ax.loglog(mids, kw[0] / widths, color=C[0], label="H in H$_2$O")
        ax.loglog(mids, scale * kf[0] / widths, color=C[1], ls="--", label="free H gas")
        ax.axvline(energy, color="0.6", lw=0.8, ls=":")
        ax.set(
            xlabel="E' (eV)",
            ylabel=r"$\sigma(E \to E')$ (b/eV)",
            xlim=(1e-4, 3),
            ylim=(1e-2, None),
            title=f"incident E = {energy:g} eV (293.6 K)",
        )
        ax.legend(fontsize=8)
    save(fig, "thermal_spectra")


def fig_thermal_matrix():
    gl = cached_groups("H1", 293.6, thermal="HinH2O")
    first = XMAS_172.find(1.0)
    m = gl.scatter_elastic[0][first:, first:]
    fig, ax = plt.subplots(figsize=(5.2, 4.4))
    im = ax.imshow(
        np.where(m > 1e-4, m, np.nan),
        norm=matplotlib.colors.LogNorm(vmin=1e-4, vmax=np.nanmax(m)),
        cmap="viridis",
        interpolation="nearest",
        extent=(first - 0.5, XMAS_172.size - 0.5, XMAS_172.size - 0.5, first - 0.5),
    )
    ax.plot([first, XMAS_172.size], [first, XMAS_172.size], color="w", lw=0.6)
    ax.set(
        xlabel="to group g'",
        ylabel="from group g",
        title="H in H$_2$O, $P_0$ below 1 eV (b)\n(left of the diagonal: upscattering)",
    )
    fig.colorbar(im, ax=ax, shrink=0.85)
    save(fig, "thermal_matrix")


# ---------------------------------------------------------------- chapter 12
def _pincell(shielding):
    from fission.pincell.run import run

    return run("1A", shielding)


def fig_pincell_spectrum():
    edges = np.array(XMAS_172.boundaries[::-1])
    du = np.array([XMAS_172.lethargy_width(g) for g in range(XMAS_172.size)])
    fig, ax = plt.subplots(figsize=(7.0, 4.0))
    for shielding, color, label in (
        ("none", C[1], "infinitely dilute"),
        ("nr", C[0], "self-shielded (NR)"),
    ):
        sol = _pincell(shielding)
        ax.stairs((sol.flux / du)[::-1], edges, color=color, lw=1.2, label=label)
    ax.set(
        xscale="log",
        xlabel="E (eV)",
        ylabel="flux per unit lethargy (per fission neutron)",
        xlim=(1e-4, 2e7),
        title="Homogenized VERA 1A pin cell: infinite-medium spectrum",
    )
    ax.legend()
    save(fig, "pincell_spectrum")


def fig_pincell_capture():
    edges = np.array(XMAS_172.boundaries[::-1])
    du = np.array([XMAS_172.lethargy_width(g) for g in range(XMAS_172.size)])
    fig, ax = plt.subplots(figsize=(7.0, 4.0))
    for shielding, color, label in (
        ("none", C[1], "infinitely dilute"),
        ("nr", C[0], "self-shielded (NR)"),
    ):
        sol = _pincell(shielding)
        m = sol.mixture
        ((c, lib),) = [
            (c, lib)
            for c, lib in zip(m.components, m.libraries, strict=True)
            if c.nuclide == "U238"
        ]
        rate = c.density * lib.sigma[102] * sol.flux
        ax.stairs(
            (rate / du)[::-1],
            edges,
            color=color,
            lw=1.2,
            label=f"{label}: {rate.sum():.3f} per fission neutron",
        )
    ax.set(
        xscale="log",
        yscale="log",
        xlabel="E (eV)",
        ylabel="U-238 capture per unit lethargy",
        xlim=(1e-3, 2e5),
        ylim=(1e-4, None),
        title="Where U-238 captures neutrons (VERA 1A, homogenized)",
    )
    ax.legend(fontsize=8)
    save(fig, "pincell_capture")


if __name__ == "__main__":
    for make in (
        fig_h1_capture,
        fig_u238_overview,
        fig_interpolation_laws,
        fig_hard_sphere,
        fig_slbw_667,
        fig_mlbw_zr90,
        fig_rm_vs_slbw,
        fig_adaptive_grid,
        fig_points_per_decade,
        fig_doppler_667,
        fig_thermal_rise,
        fig_corner_rounding,
        fig_shielded_ri,
        fig_nu_multiplicity,
        fig_fission_spectrum,
        fig_elastic_angular,
        fig_continuum,
        fig_weight_function,
        fig_group_capture,
        fig_scatter_matrices,
        fig_chi_groups,
        fig_thermal_xs,
        fig_thermal_spectra,
        fig_thermal_matrix,
        fig_pincell_spectrum,
        fig_pincell_capture,
    ):
        make()
