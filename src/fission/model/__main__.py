"""Check a model file and summarize it.

uv run python -m fission.model vera-p1            # models/vera-p1.toml
uv run python -m fission.model my.toml --state 1A # number densities at a state
"""

from __future__ import annotations

import argparse

from .homogenize import material_areas
from .load import load_model
from .materials import number_densities
from .spec import Model


def summary(model: Model, state: str | None = None) -> str:
    lines = [f"{model.name}: {model.description}".rstrip(": ")]
    lines.append(f"  materials:  {', '.join(model.materials)}")
    for p in model.pins.values():
        rings = ", ".join(
            f"{m} < {r:g}" for m, r in zip(p.materials, p.radii, strict=False)
        )
        if p.radii:
            lines.append(f"  pin {p.name}: {rings}, {p.materials[-1]} outside")
        else:
            lines.append(f"  pin {p.name}: {p.materials[-1]} only (no rings)")
    for lat in model.lattices.values():
        if lat.material is not None:
            lines.append(f"  lattice {lat.name}: homogeneous {lat.material}")
            continue
        assert lat.pins is not None
        n = len(lat.pins)
        areas = material_areas(model, lat)
        total = sum(areas.values())
        frac = ", ".join(f"{m} {a / total:.1%}" for m, a in areas.items())
        lines.append(
            f"  lattice {lat.name}: {n}x{n} at pitch {lat.pitch:g} cm; area: {frac}"
        )
    for a in model.assemblies.values():
        stack = " / ".join(f"{s.lattice} {s.height:g}" for s in a.segments)
        lines.append(f"  assembly {a.name} (bottom to top, cm): {stack}")
    if model.core is not None:
        full = model.core.full_map()
        lines.append(
            f"  core: {len(full)}x{len(full[0])} positions ({model.core.symmetry} map), "
            f"{sum(a is not None for r in full for a in r)} assemblies"
        )
    for s in model.states.values():
        lines.append(
            f"  state {s.name}: fuel {s.fuel_temperature:g} K, moderator "
            f"{s.moderator_temperature:g} K at {s.moderator_density:g} g/cm3, "
            f"{s.boron_ppm:g} ppm boron"
        )
    if state is not None:
        st = model.states[state]
        lines.append(f"  number densities at {state} (atoms/b-cm):")
        for m in model.materials.values():
            nd = number_densities(m, st)
            body = ", ".join(f"{k} {v:.5e}" for k, v in nd.items())
            lines.append(f"    {m.name} ({st.temperature(m.temperature):g} K): {body}")
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(prog="python -m fission.model", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)  # fmt: skip
    ap.add_argument("model", help="a TOML file, or a name in models/")
    ap.add_argument("--state", help="also print number densities at this state")
    args = ap.parse_args(argv)
    print(summary(load_model(args.model), args.state))


if __name__ == "__main__":
    main()
