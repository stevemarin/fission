"""Human-readable names for ENDF reaction numbers (MT), incident neutrons."""

_NAMES = {
    1: "total",
    2: "elastic",
    3: "nonelastic",
    4: "(n,n') total",
    5: "(n,anything)",
    11: "(n,2nd)",
    16: "(n,2n)",
    17: "(n,3n)",
    18: "fission (total)",
    19: "(n,f) first-chance",
    20: "(n,nf) second-chance",
    21: "(n,2nf) third-chance",
    22: "(n,n alpha)",
    23: "(n,n 3alpha)",
    24: "(n,2n alpha)",
    28: "(n,np)",
    32: "(n,nd)",
    33: "(n,nt)",
    37: "(n,4n)",
    38: "(n,3nf) fourth-chance",
    41: "(n,2np)",
    91: "(n,n') continuum",
    101: "disappearance",
    102: "(n,gamma) capture",
    103: "(n,p) total",
    104: "(n,d) total",
    105: "(n,t) total",
    106: "(n,3He) total",
    107: "(n,alpha) total",
    108: "(n,2 alpha)",
    111: "(n,2p)",
    112: "(n,p alpha)",
    452: "nu-bar total",
    455: "nu-bar delayed",
    456: "nu-bar prompt",
    458: "fission energy release",
    849: "(n,alpha) continuum",
}
_LEVEL_SERIES = [  # (first MT, last MT, particle): MT = first + k is level k, last = continuum
    (600, 649, "p"),
    (650, 699, "d"),
    (700, 749, "t"),
    (750, 799, "3He"),
    (800, 849, "alpha"),
]


def mt_name(mt: int) -> str:
    if mt in _NAMES:
        return _NAMES[mt]
    if 51 <= mt <= 90:
        return f"(n,n') level {mt - 50}"
    if 875 <= mt <= 890:
        return f"(n,2n) level {mt - 875}"
    for first, last, particle in _LEVEL_SERIES:
        if first <= mt < last:
            return f"(n,{particle}) level {mt - first}"
        if mt == last:
            return f"(n,{particle}) continuum"
    return f"MT{mt}"
