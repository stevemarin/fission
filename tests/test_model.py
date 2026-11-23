"""The model file: parsing, checks, number densities, homogenization."""

import math
import tomllib

import pytest

from fission.data import EndfLibrary
from fission.model import (
    Isotopics,
    Material,
    ModelError,
    State,
    homogenize,
    load_model,
    material_areas,
    number_densities,
    parse,
)
from fission.model.materials import AVOGADRO, BARN, atom_fractions
from fission.pincell.run import UNSUPPORTED, components

# A toy core: every number is made up; it exercises all layers.
TOY = """
[model]
format = 1
name = "toy"

[materials.fuel]
temperature = "fuel"
atom_densities = { U235 = 7e-4, U238 = 2.2e-2, O16 = 4.6e-2 }

[materials.clad]
temperature = "moderator"
atom_densities = { Zr90 = 4.3e-2 }

[materials.water]
temperature = "moderator"
density = "moderator"
formula = { H = 2, O = 1 }
atom_percent.H = { H1 = 100.0 }
atom_percent.O = { O16 = 100.0 }
soluble_boron.atom_percent.B = { B10 = 20.0, B11 = 80.0 }
thermal = { H1 = "HinH2O" }

[pins.rod]
radii = [0.41, 0.475]
materials = ["fuel", "clad", "water"]

[pins.guide]
radii = [0.56, 0.60]
materials = ["water", "clad", "water"]

[lattices.A]
pitch = 1.26
legend = { F = "rod", G = "guide" }
map = '''
F F F
F G F
F F F
'''

[lattices.reflector]
material = "water"

[assemblies.A]
axial = [
  { lattice = "reflector", height = 20.0 },
  { lattice = "A", height = 360.0 },
  { lattice = "reflector", height = 20.0 },
]

[assemblies.R]
axial = [{ lattice = "reflector", height = 400.0 }]

[core]
assembly_pitch = 4.0
symmetry = "quarter"
legend = { A = "A", R = "R" }
map = '''
A A R
A R .
R . .
'''

[states.hzp]
fuel_temperature = 565.0
moderator_temperature = 565.0
moderator_density = 0.743
boron_ppm = 1300
"""


def toy(**patch) -> dict:
    """The toy model as a dict, with dotted-path replacements (None deletes)."""
    data = tomllib.loads(TOY)
    for path, value in patch.items():
        *parents, last = path.split("__")
        node = data
        for p in parents:
            node = node[p]
        if value is None:
            del node[last]
        else:
            node[last] = value
    return data


def test_toy_model_parses_all_layers():
    m = parse(toy())
    assert list(m.materials) == ["fuel", "clad", "water"]
    assert m.lattices["A"].width == pytest.approx(3 * 1.26)
    assert m.assemblies["A"].height == 400.0
    assert m.core is not None
    full = m.core.full_map()
    assert len(full) == 5 and all(len(r) == 5 for r in full)
    assert full[2] == ("R", "A", "A", "A", "R")  # the centre row, mirrored
    assert full[0] == (None, None, "R", None, None)
    assert sum(a == "A" for r in full for a in r) == 5
    assert m.materials["water"].thermal_law("H1") == "HinH2O"
    assert m.materials["water"].thermal_law("O16") == "free-gas"


@pytest.mark.parametrize(
    ("patch", "message"),
    [
        ({"model__format": 2}, "model.format: this code reads format 1"),
        ({"pins__rod__radius": [0.4]}, "pins.rod: unknown key(s) ['radius']"),
        ({"pins__rod__radii": [0.475, 0.41]}, "must be positive and increasing"),
        ({"pins__rod__materials": ["fuel", "clad"]}, "2 radii need 3 materials"),
        ({"pins__rod__materials": ["fuel", "steel", "water"]}, "unknown material 'steel'"),
        ({"lattices__A__pitch": 0.9}, "does not fit in the pitch 0.9 cm"),
        ({"lattices__A__map": "F F\nF G F\nF F F"}, "row 2 has 3 entries, row 1 has 2"),
        ({"lattices__A__map": "F F\nF X"}, "row 2: 'X' is not in the legend"),
        ({"lattices__A__geometry": "hex"}, "'hex' is reserved but not implemented yet"),
        ({"core__symmetry": "eighth"}, "'eighth' is reserved but not implemented yet"),
        ({"core__assembly_pitch": 3.0}, "is 3.78 cm wide, more than the assembly pitch"),
        ({"assemblies__R__axial": [{"lattice": "reflector", "height": 380.0}]},
         "same total height"),
        ({"assemblies__A__axial": [{"lattice": "B", "height": 1.0}]},
         "unknown lattice 'B'"),
        ({"materials__water__atom_percent__H": {"H1": 99.0}}, "add up to 99, not 100"),
        ({"materials__water__atom_percent__O": {"H1": 100.0}}, "is not an isotope of O"),
        ({"materials__water__formula": {"H": 2, "O": 1, "N": 1}},
         "no isotopics for element 'N'"),
        ({"materials__water__thermal": {"Zr90": "ZrinZrH"}}, "nuclide not in this material"),
        ({"materials__fuel__temperature": "clad"}, "'clad' is not one of"),
        ({"materials__fuel__density": 10.0}, "not allowed together with atom_densities"),
        ({"states__hzp__boron_ppm": "1300"}, "expected a number, got a string"),
        ({"states__hzp__fuel_temperature": None}, "missing required key 'fuel_temperature'"),
        ({"materials__fuel__shielded": "yes"}, "expected true or false, got a string"),
        ({"core__banks": {"rod": "rcca", "map": "D"}}, "unknown control rod 'rcca'"),
    ],
)  # fmt: skip
def test_mistakes_are_reported_with_their_location(patch, message):
    with pytest.raises(ModelError, match=None) as err:
        parse(toy(**patch))
    assert message in str(err.value)


# ---------------------------------------------------------------- number densities
MASSES = {
    "H1": 1.0,
    "O16": 16.0,
    "B10": 10.0,
    "B11": 11.0,
    "U235": 235.0,
    "U238": 238.0,
}


def _state(density=1.0, ppm=0.0):
    return State("s", 300.0, 300.0, density, ppm)


def _water(boron=None):
    return Material(
        "w", "moderator", density="moderator", formula={"H": 2, "O": 1},
        isotopics={"H": Isotopics("atom", {"H1": 100.0}),
                   "O": Isotopics("atom", {"O16": 100.0})},
        soluble_boron=boron,
    )  # fmt: skip


def test_compound_number_densities():
    n = number_densities(_water(), _state(1.8), MASSES.__getitem__)
    molecules = 1.8 / 18.0 * AVOGADRO * BARN  # rho N_A / M
    assert n["O16"] == pytest.approx(molecules, rel=1e-14)
    assert n["H1"] == pytest.approx(2 * molecules, rel=1e-14)


def test_soluble_boron_is_parts_per_million_by_weight():
    boron = Isotopics("atom", {"B10": 20.0, "B11": 80.0})
    n = number_densities(_water(boron), _state(1.0, 2000.0), MASSES.__getitem__)
    grams = {k: v / AVOGADRO / BARN * MASSES[k] for k, v in n.items()}  # g/cm3
    assert sum(grams.values()) == pytest.approx(1.0, rel=1e-14)
    assert (grams["B10"] + grams["B11"]) == pytest.approx(2000e-6, rel=1e-14)
    assert n["B10"] / (n["B10"] + n["B11"]) == pytest.approx(0.2, rel=1e-14)


def test_weight_and_atom_percent_are_consistent():
    by_weight = Isotopics("weight", {"U235": 3.0, "U238": 97.0})
    frac, mean = atom_fractions(by_weight, MASSES.__getitem__)
    by_atom = Isotopics("atom", {k: 100 * v for k, v in frac.items()})
    frac2, mean2 = atom_fractions(by_atom, MASSES.__getitem__)
    assert mean2 == pytest.approx(mean, rel=1e-14)
    assert frac2 == pytest.approx(frac, rel=1e-14)
    assert frac["U235"] * 235.0 / mean == pytest.approx(0.03, rel=1e-14)  # back to wt


def test_homogenization_conserves_atoms():
    m = parse(toy())
    lat = m.lattices["A"]
    areas = material_areas(m, lat)
    assert sum(areas.values()) == pytest.approx((3 * 1.26) ** 2, rel=1e-14)
    assert areas["fuel"] == pytest.approx(8 * math.pi * 0.41**2, rel=1e-14)
    comps = homogenize(m, "A", "hzp", MASSES.__getitem__)
    u238 = sum(c.density for c in comps if c.nuclide == "U238") * (3 * 1.26) ** 2
    assert u238 == pytest.approx(2.2e-2 * areas["fuel"], rel=1e-14)
    o16 = {c.region: c for c in comps if c.nuclide == "O16"}
    assert set(o16) == {"fuel", "water"}  # two components, two temperatures
    assert {c.thermal for c in comps if c.nuclide == "H1"} == {"HinH2O"}


# ---------------------------------------------------------------- VERA problem 1
# Table P1-4 of CASL-U-2012-0131-004 Rev. 4: the number densities the model
# computes from mass density, isotopics and boron must reproduce it.
TABLE_P1_4 = {
    ("fuel", "1A"): {"U234": 6.11864e-06, "U235": 7.18132e-04, "U236": 3.29861e-06,
                     "U238": 2.21546e-02, "O16": 4.57642e-02},
    ("moderator", "1A"): {"O16": 2.48112e-02, "H1": 4.96224e-02,
                          "B10": 1.07070e-05, "B11": 4.30971e-05},
    ("moderator", "1B"): {"O16": 2.20729e-02, "H1": 4.41459e-02,
                          "B10": 9.52537e-06, "B11": 3.83408e-05},
}  # fmt: skip


@pytest.fixture(scope="module")
def vera():
    if not (
        EndfLibrary().root / "archives" / "neutrons-version.VIII.1.tar.gz"
    ).exists():
        pytest.skip("ENDF/B-VIII.1 neutron data not downloaded (atomic masses)")
    return load_model("vera-p1")


def test_vera_number_densities_match_table_p1_4(vera):
    for (material, state), table in TABLE_P1_4.items():
        n = number_densities(vera.materials[material], vera.states[state])
        assert set(n) == set(table)
        for nuclide, ref in table.items():
            # Table P1-4 prints 6 digits; boron is rounded at 5e-6
            assert n[nuclide] == pytest.approx(ref, rel=1e-5), (material, nuclide)


def test_vera_cell(vera):
    lat = vera.lattices["cell"]
    areas = material_areas(vera, lat)
    assert sum(areas.values()) == pytest.approx(1.26**2, rel=1e-14)
    assert areas["fuel"] / 1.26**2 == pytest.approx(0.332, abs=5e-4)  # chapter 12
    assert areas["moderator"] / 1.26**2 == pytest.approx(0.554, abs=5e-4)
    all_comps = homogenize(vera, "cell", "1A")
    comps = components(vera, "cell", "1A")
    assert len(all_comps) - len(comps) == 3  # Fe-57, Fe-58, Cr-53 left out
    assert len(comps) == 36 and len({c.nuclide for c in comps}) == 35
    clad = vera.materials["clad"].atom_densities
    assert clad is not None
    left_out = sum(clad.get(n, 0.0) for n in UNSUPPORTED)
    assert left_out < 2e-5 and left_out / sum(clad.values()) < 5e-4
    temps = {c.region: c.temperature for c in homogenize(vera, "cell", "1C")}
    assert temps == {"fuel": 900.0, "gap": 600.0, "clad": 600.0, "moderator": 600.0}
    assert vera.states["1A"].reference is not None
    assert vera.states["1A"].reference.keff == 1.187038


def test_a_model_can_extend_another(tmp_path):
    from fission.model import load_model

    (tmp_path / "base.toml").write_text(TOML_BASE := TOY)
    (tmp_path / "child.toml").write_text(
        '[model]\nextends = "base"\nname = "child"\n'
        "[materials.fuel]\natom_densities = { U235 = 1e-3, U238 = 2.2e-2, O16 = 4.6e-2 }\n"
        "[states.hfp]\nfuel_temperature = 900.0\nmoderator_temperature = 580.0\n"
        "moderator_density = 0.7\nboron_ppm = 1000\n"
    )
    m = load_model(tmp_path / "child.toml")
    assert m.name == "child"
    assert m.materials["fuel"].atom_densities["U235"] == 1e-3  # overridden
    assert m.materials["clad"].atom_densities == {"Zr90": 4.3e-2}  # inherited
    assert set(m.states) == {"hzp", "hfp"}
    assert TOML_BASE


def test_wbn1_enrichment_zones():
    """The WBN1 loading pattern: 193 assemblies of three fuel types, each fuel
    with the uranium isotopics of VERA Table 17 (U-234 = 0.007731 w^1.0837,
    U-236 = 0.0046 w, U-238 the balance)."""
    zoned = load_model("wbn1-zoned")
    cells = [a for row in zoned.core.full_map() for a in row if a is not None]
    fuel = [a for a in cells if a in "abc"]
    assert len(fuel) == 193
    assert set(cells) == {"a", "b", "c", "side"}
    for name, enrichment in (("211", 2.11), ("262", 2.619), ("310", 3.1)):
        u = load_model(f"wbn1-fuel-{name}").materials["fuel"].isotopics["U"]
        assert u.basis == "weight"
        assert u.percent["U235"] == enrichment
        assert u.percent["U234"] == pytest.approx(
            0.007731 * enrichment**1.0837, abs=6e-5
        )
        assert u.percent["U236"] == pytest.approx(0.0046 * enrichment, abs=6e-5)
        assert sum(u.percent.values()) == pytest.approx(100.0, abs=1e-4)
