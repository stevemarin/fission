"""The core at power: nodal diffusion with thermal-hydraulic feedback (M4b).

Every node's few-group constants depend on its state (fuel temperature,
water temperature and density, boron), and the state depends on the power
the constants produce. The two are iterated (Picard):

1. evaluate each node's constants from its material's table at its state;
2. solve the core with NEM;
3. power per node from kappa Sigma_f phi, scaled to the core power;
4. per channel (a column of nodes), the coolant enthalpy rises with the
   power below: h = h_in + (power up to the node centre) / mass flow; the
   water temperature and density follow from IAPWS-IF97 (``water``);
5. per fuel node, the rod temperatures from its linear heat rate
   (``fuel``); the Doppler temperature is the state's fuel temperature;
6. repeat until the eigenvalue and every temperature stop changing.

Reflector nodes take the water state of their channel: the radial
reflector the inlet water, the top reflector the water leaving the fuel.
A channel is one node column, with the flow of an assembly scaled by the
column's share of an assembly's area (half on a quarter-core symmetry
line); there is no cross flow between channels.
"""

from __future__ import annotations

import time
from collections.abc import Callable, Sequence
from dataclasses import dataclass, field

import numpy as np

from fission.model import BRANCH_AXES, Model
from fission.nodal import nem
from fission.nodal.mesh import CoreMesh, core_mesh, with_states

from . import convection, fuel, water

EV = 1.602176634e-19  # J


@dataclass
class Result:
    k: float
    flux: np.ndarray  # [G, nz, ny, nx]
    power: np.ndarray  # W per node
    states: dict[str, np.ndarray]  # BRANCH_AXES -> [nz, ny, nx]
    centre: np.ndarray  # fuel centreline temperature per node (K), NaN outside fuel
    film: np.ndarray  # clad-to-water heat transfer coefficient (W/(m^2 K)), NaN if none
    mesh: CoreMesh
    history: list[tuple[float, float, float]] = field(
        default_factory=list
    )  # k, dT_f, dT_m
    seconds: float = 0.0


def node_power(mesh: CoreMesh, flux: np.ndarray, total: float) -> np.ndarray:
    """W per node, proportional to kappa Sigma_f phi V, summing to `total`."""
    kappa = mesh.field(lambda x: x.kappa_fission if x.kappa_fission is not None
                       else x.nu_fission)  # fmt: skip
    p = (kappa * flux).sum(axis=0) * mesh.volumes()
    return p * total / p.sum()


def channels(mesh: CoreMesh, power: np.ndarray, flow_per_area: float, h_in: float,
             pressure: float) -> tuple[np.ndarray, np.ndarray]:  # fmt: skip
    """Water temperature (K) and density (g/cm^3) at every node centre."""
    dx, dy, _ = mesh.widths
    area = dy[:, None] * dx[None, :]  # cm^2 per node column
    flow = flow_per_area * area  # kg/s
    below = np.cumsum(power, axis=0) - power / 2  # W up to each node centre
    h = h_in + below / flow[None]
    t = np.zeros_like(h)
    rho = np.zeros_like(h)
    for index in np.ndindex(h.shape):
        t[index] = water.temperature(pressure, float(h[index]))
        if t[index] > 623.15:
            raise ValueError(f"water at {t[index]:.1f} K: beyond IF97 region 1 "
                             "(boiling is not modelled)")  # fmt: skip
        rho[index] = water.density(pressure, t[index]) / 1000.0
    return t, rho


def isothermal(
    model: Model,
    temperatures: Sequence[float],
    refine: int = 1,
    axial: float | None = None,
    banks: dict[str, float] | None = None,
    boron: float | None = None,
) -> list[float]:
    """k of the core at zero power with fuel and water all at one
    temperature T (water density from IF97 at the core pressure, boron
    from the [core] state unless given, ppm), for each T. The slope of the
    reactivity is the isothermal temperature coefficient that start-up
    tests measure. The reflector keeps its fixed constants. banks: control
    rods inserted (bank -> steps withdrawn)."""
    core = model.core
    if core is None or core.thermal is None or core.state is None:
        raise ValueError("the core needs [core] state = ... and [core.thermal]")
    base = core_mesh(model, refine, axial, banks)
    if boron is None:
        boron = model.states[core.state].boron_ppm
    out = []
    for t in temperatures:
        rho = water.density(core.thermal.pressure, t) / 1000.0
        states = {
            "fuel_temperature": np.full(base.shape, t),
            "moderator_temperature": np.full(base.shape, t),
            "moderator_density": np.full(base.shape, rho),
            "boron_ppm": np.full(base.shape, boron),
        }
        out.append(nem.solve(with_states(base, states)).k)
    return out


def bank_worths(model: Model, refine: int = 1, axial: float | None = None,
                boron: float | None = None) -> tuple[float, dict[str, float]]:  # fmt: skip
    """Worth (pcm) of each control rod bank of the model, fully inserted
    alone, at zero power with fuel and water at the inlet temperature:
    rho(all rods out) - rho(bank in). Returns (k all rods out, worths)."""
    core = model.core
    if core is None or core.banks is None or core.thermal is None:
        raise ValueError("needs [core.banks] and [core.thermal]")
    t = [core.thermal.inlet_temperature]
    k_out = isothermal(model, t, refine, axial, boron=boron)[0]
    worths = {}
    for bank in core.banks.positions:
        k_in = isothermal(model, t, refine, axial, {bank: 0.0}, boron)[0]
        worths[bank] = 1e5 * (1 / k_in - 1 / k_out)
    return k_out, worths


def solve(model: Model, refine: int = 1, axial: float | None = None,
          tol_k: float = 1e-6, tol_t: float = 0.05, max_iterations: int = 50,
          relax: float = 1.0, power_level: float = 1.0, verbose: bool = False,
          banks: dict[str, float] | None = None,
          boron: float | None = None) -> Result:  # fmt: skip
    """Coupled NEM + thermal-hydraulics of the model's core at its
    operating conditions ([core] state and [core.thermal]); power_level
    scales the rated power (the flow stays at its rated value); banks:
    control rods inserted (bank -> steps withdrawn); boron (ppm) replaces
    the state's."""
    t0 = time.perf_counter()
    core = model.core
    if core is None or core.thermal is None or core.state is None:
        raise ValueError("the core needs [core] state = ... and [core.thermal]")
    th = core.thermal
    base = core_mesh(model, refine, axial, banks)
    state = model.states[core.state]
    shape = base.shape
    fuel_node = base.field(lambda x: x.nu_fission).sum(axis=0) > 0
    pitch = core.assembly_pitch
    dz = base.widths[2]
    # the fuel stack height, from the column with the most fuel nodes
    j, i = np.unravel_index(np.argmax(fuel_node.sum(axis=0)), shape[1:])
    fuel_height = float(dz[fuel_node[:, j, i]].sum())  # cm
    # power density: the assembly power over the fuel volume of one assembly
    density = power_level * th.assembly_power / (pitch * pitch * fuel_height)  # W/cm^3
    total = density * float((base.volumes() * fuel_node).sum())
    flow_per_area = th.assembly_flow / (pitch * pitch)  # kg/s per cm^2
    rods_per_area = th.rods / (pitch * pitch)
    film_on = th.flow_area is not None
    if film_on:  # clad-to-water heat transfer: the same mass flux in every channel
        assert th.flow_area is not None and th.hydraulic_diameter is not None
        assert th.rod_pitch is not None
        mass_flux = th.assembly_flow / (th.flow_area * 1e-4)  # kg/(m^2 s)
        d_h = th.hydraulic_diameter / 100.0
        p_over_d = th.rod_pitch / (2 * th.clad_outer)
    film = np.full(shape, np.nan)
    h_in = water.enthalpy(th.pressure, th.inlet_temperature)
    rho_in = water.density(th.pressure, th.inlet_temperature) / 1000.0
    states = {
        "fuel_temperature": np.full(shape, state.fuel_temperature),
        "moderator_temperature": np.full(shape, th.inlet_temperature),
        "moderator_density": np.full(shape, rho_in),
        "boron_ppm": np.full(shape, state.boron_ppm if boron is None else boron),
    }
    history = []
    centre = np.full(shape, np.nan)
    k_old = 0.0
    r = None
    power = np.zeros(shape)
    for it in range(1, max_iterations + 1):
        mesh = with_states(base, states)
        r = nem.solve(mesh, start=r)  # from the previous iteration
        power = node_power(mesh, r.flux, total)
        t_mod, rho = channels(mesh, power, flow_per_area, h_in, th.pressure)
        t_fuel = states["fuel_temperature"].copy()
        dx, dy, _ = mesh.widths
        area = dy[:, None] * dx[None, :]
        for k, j, i in zip(*np.nonzero(fuel_node), strict=True):
            q = (
                power[k, j, i]
                * th.fuel_heat_fraction
                / (rods_per_area * area[j, i] * dz[k] / 100.0)
            )  # W/m, fmt: skip
            h_film = None
            if film_on:
                h_film = convection.film_coefficient(th.pressure, float(t_mod[k, j, i]),
                                                     mass_flux, d_h, p_over_d)  # fmt: skip
                film[k, j, i] = h_film
            temps = fuel.rod(q, float(t_mod[k, j, i]), th.pellet_radius / 100,
                             th.clad_inner / 100, th.clad_outer / 100, th.gap_conductance,
                             h_film)  # fmt: skip
            t_fuel[k, j, i] = temps.doppler
            centre[k, j, i] = temps.centre
        d_fuel = float(np.max(np.abs(t_fuel - states["fuel_temperature"])[fuel_node]))
        d_mod = float(np.max(np.abs(t_mod - states["moderator_temperature"])))
        states["fuel_temperature"] = np.where(
            fuel_node, (1 - relax) * states["fuel_temperature"] + relax * t_fuel,
            states["fuel_temperature"],
        )  # fmt: skip
        states["moderator_temperature"] = t_mod
        states["moderator_density"] = rho
        history.append((r.k, d_fuel, d_mod))
        if verbose:
            print(f"  feedback {it:2d}: k = {r.k:.6f}, max change fuel {d_fuel:.2f} K, "
                  f"water {d_mod:.3f} K", flush=True)  # fmt: skip
        if abs(r.k - k_old) < tol_k and d_fuel < tol_t and d_mod < tol_t:
            break
        k_old = r.k
    assert r is not None
    assert set(states) == set(BRANCH_AXES)
    return Result(r.k, r.flux, power, states, centre, film, mesh, history,
                  time.perf_counter() - t0)  # fmt: skip


def critical_boron(
    k_of: Callable[[float], float],
    start: tuple[float, float] = (1000.0, 1300.0),
    tol: float = 0.05,
    max_iterations: int = 20,
) -> tuple[float, float, int]:
    """The boron concentration (ppm) that makes k_of(boron) = 1.

    Secant iteration on the reactivity rho = 1 - 1/k, which is close to
    linear in boron: soluble boron is a thin 1/v absorber spread through
    the water, and over a few hundred ppm its worth per ppm barely changes
    (the tables interpolate it linearly between their boron nodes). Stops
    when the step is below tol ppm. Returns (boron, k there, solves)."""
    b0, b1 = start
    r0, r1 = 1 - 1 / k_of(b0), 1 - 1 / (k := k_of(b1))
    for n in range(3, max_iterations + 3):
        b2 = b1 - r1 * (b1 - b0) / (r1 - r0)
        b0, r0 = b1, r1
        k = k_of(b2)
        b1, r1 = b2, 1 - 1 / k
        if abs(b1 - b0) < tol:
            return b1, k, n
    raise RuntimeError(
        f"no critical boron after {max_iterations} steps (last {b1:.1f} ppm)"
    )
