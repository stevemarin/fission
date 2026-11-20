# fission

A lattice-physics code written from scratch as a **learning project**: raw
ENDF/B-VIII.1 → pointwise cross sections → 2D pin-cell transport (MOC, later
CP) → few-group constants → 3D nodal diffusion. No NJOY/OpenMC in the
pipeline (they are installed under tools/ for independent checks: `make tools`). Python package:
`fission`. An earlier NEM/CMFD nodal diffusion code lived on branch
`claude-nodal` (deleted by the user; the only copy left is the gitignored
extract in `endf-data/nodal-src`, used by `make nodal-check`).

## Working rules

- **Python first, then Rust.** Write every kernel as clean, readable Python
  (the reference), then port the hot path to Rust and add parity tests
  (agreement to ~1e-12, identical grids for adaptive algorithms). The Python
  version stays as the documentation of the method.
- **Ask before adding any Python or Rust dependency.** Current deps: numpy;
  dev: matplotlib, mypy, pytest, ruff; Rust: pyo3, numpy.
- **Go step by step and explain.** The user is learning the theory; explain
  the physics of each step, not just the code. One roadmap step at a time.
- **Each step extends the theory docs** (`docs/theory/`, see below).
- Commit only when asked.

## Commands (use `make`, it is allow-listed)

```sh
make tools        # NJOY 2016.78 + OpenMC 0.16.0 (conda-forge) in tools/env, for independent checks
make check        # lint + format-check + mypy + pytest + cargo test (run before reporting done)
make test ARGS="-k reich -x"   # pytest; slow tests are deselected by default
make test-slow    # full-resolution physics checks (minutes)
make fix          # ruff --fix + format
make data         # download pin-cell evaluations into endf-data/
make tables       # build/load cached pointwise tables for pin-cell nuclides
make doppler-tables TEMPS="293.6 900"  # broadened tables (endf-data/doppler/)
make pincell CASES=1A SHIELD=nr   # k-inf of the homogenized VERA 1 cell (none|nr; MODEL=vera-p1)
make pin STATES="1A 1B" PIN_SHIELD=sd PIN_ARGS="--rings --sph --clad --cp"  # heterogeneous VERA pin, MOC (+ CP cross-check; ~3.5 min/state at the converged defaults, PIN_ARGS+="--azim 64 --fuel-rings 3" for a quick look)
make assembly STATES=2A PIN_ARGS="--rings --sph --clad"  # VERA 17x17 lattice: k + pin powers vs KENO (~11 min/state)
make fewgroup STATES=1A PIN_ARGS="--rings --sph"  # 2-group constants (homogenize, B1, condense) -> endf-data/fewgroup/
make branches BRANCHES=vera  # state matrix (M2): constants at every branch point -> endf-data/branches/ (~76 min infinite; BRANCH_ARGS="--solver moc --rings --sph --clad")
make xs-tables BRANCHES=vera  # M3: table for the nodal code (endf-data/tables/) + tests: leave-one-out, PRESS-selected fit, check points (make branches BRANCH_ARGS="--check 10")
make nodal-check STATE=1A  # old nodal code (endf-data/nodal-src copy, tools env) on those constants
make openmc-library  # NJOY -> OpenMC HDF5 of our ENDF files (endf-data/openmc/, ~10 min; check only)
make openmc-pin STATES=1A PIN_ARGS="--rings --sph"  # OpenMC CE reference of the pin vs this code (~15 min/state)
make model MODEL=vera-p1 STATE=1A  # check a model file, print number densities
make tracks MODEL=vera-p1 AZIM=64 SPACING=0.01  # cyclic tracks + ray tracing, area check
make moc CORE=c5g7 MOC_ARGS="--azim 128"  # MOC+CMFD eigenvalue + pin powers vs C5G7 (~1 min)
make groups       # XMAS-172 libraries, free-gas thermal + H in H2O (endf-data/groups/); make library NUC=U238 T=293.6
make docs         # docs/build/theory.html (pandoc, MathML, offline)
make figures      # regenerate docs/figures/*.svg (deterministic output)
make show NUC=U235 | make scan NUC=U238 E=6.674 | make rm | make doppler | make secondary   # quick CLIs
```

- NJOY and OpenMC live in an isolated micromamba env (`tools/`, gitignored),
  never in the uv env or the pipeline. Run them with
  `tools/bin/micromamba run -r tools/mamba -p tools/env <cmd>` (Makefile
  variable `IN_TOOLS`); exchange data with checking scripts through files.
- Everything runs through `uv run`. Editing `rust/**/*.rs` or `Cargo.toml`
  makes the next `uv run` rebuild the extension (`[tool.uv] cache-keys`); no
  manual `maturin develop` needed.
- Lint/type scope is `src tests docs/figures`.

## Layout

- `src/fission/processing/`: stages that take 0 K tables as input:
  `doppler.py` (+ `cache.py`), group averaging `groupr.py` with `groups.py`
  (XMAS-172 from NJOY groupr.f90), `weights.py` (NJOY iwt=4, slope-matched
  breaks), `library.py` (GroupLibrary cache + CLI, `thermal=` option),
  `thermal.py` (free-gas and S(alpha,beta) kernels, `with_thermal`). Kept outside `endf/` so
  edits do not invalidate the 0 K table cache.
- `models/*.toml`: reactor models (M1). `vera-p1.toml` = VERA problem 1
  (CASL-U-2012-0131-004 Rev. 4: geometry, compositions, states, KENO
  references). Format: materials -> pins -> lattices (2D xy) -> assemblies
  (z stacks) -> core map, plus states; see docs ch. 13.
- `src/fission/model/`: `spec.py` (dataclasses), `load.py` (`load_model`,
  strict TOML checks with located errors; unknown keys are errors; reserved
  but unimplemented: hex geometry, eighth symmetry), `materials.py` (number
  densities from rho + formula + isotopics + ppm boron; masses = ENDF AWR x
  m_n), `homogenize.py` (area weighting -> `Component`).
- `src/fission/transport/` (B1): `geometry.py` (`Geometry.from_lattice`:
  FSRs = equal-area rings x sectors per cell, exact areas; cells numbered
  from bottom-left, model map rows are top-down), `tracks.py`
  (`cyclic_tracks`: OpenMOC-style cyclic angles/spacing, weights, sweep
  links fwd/bwd, reflective or vacuum per side), `trace.py` (segments via
  sorted crossings + midpoint classification; `track_areas`). Rust
  `raytrace.rs` is bitwise identical to Python. `moc.py` (B2): flat-source
  isotropic MOC, Gauss-Legendre polar, area-corrected segments, power
  iteration stopping on the extrapolated error d rho/(1-rho); Rust
  `moc.rs` sweep (threaded over tracks). `run.py`: core driver + C5G7
  reference (Tables 3-5). `Geometry.from_core` flattens the core map.
  `cmfd.py` (B3): pin-cell coarse mesh; sweeps tally net face currents
  (`trace.crossings`, signed sum W psi); D^ makes the coarse problem
  reproduce MOC exactly; D^ under-relaxed (RELAX = 0.7, plain CMFD diverges
  on C5G7's optically thick water cells); coarse power iteration with
  red-black SOR (Rust `cmfd.rs`); prolongation scales FSR and boundary
  angular fluxes. `solve(cmfd=True)` is the default; `--no-cmfd` in run.py.
- `models/c5g7.toml` + `c5g7-xs.toml`: macroscopic materials (`macroscopic =`,
  `[model] multigroup =` library file), pins without rings (`radii = []`),
  `fill`/`size` lattices, `[core] boundary` per side. Layout read from the
  spec's Figure 3 by pixel classification; XS parsed from its text.
- `src/fission/pincell/shielding.py` (C1): Dancoff C by black-pin MOC
  fixed source (ratio to the same-tracking unshadowed current),
  equivalence sigma_0 (Bell & Glasstone 2.122, 8.87-8.88; sigma_p = 4 pi a^2
  from MF2), exact escape probability P(Sigma) table by MOC, region flux
  P + (1-P) Sigma_p/Sigma_t (NR) or the ultra-fine slowing-down `march`
  (1 eV-24.8 keV, exact elastic kernels, moderator level depleted).
  `heterogeneous.py`: transport-corrected (outflow) macroscopic XS from
  nuclide libraries, MOC+CMFD VERA pin; shielding none|equivalence|cp|sd.
  C2: `ring_escape` (E_k by reciprocity), `collision_matrix` (P_j->i
  between pellet rings), multi-ring `march` (each ring its own source; Rust
  `slowing.rs`, parity 1e-12), `sph_factors` (one-group fixed-source MOC per
  group, mu <- <phi>/phi_MG); `run(ring_shielding=, sph=)`, ring regions
  are MOC materials `fuel/0`(centre)..`fuel/2` (`Geometry.split_rings`).
  `reference.py`: JSON spec for OpenMC + comparison by nuclide/ring/band;
  tools side `checks/openmc_library.py`, `checks/openmc_pin.py` (outside
  lint scope, run in tools env; results in endf-data/openmc/runs/).
- C4: `models/vera-p2.toml` (17x17, guide/instrument tubes, lattice key
  `gap` = inter-assembly half gap; gap points belong to the edge cells'
  outer regions, exact off-centre wedge areas; tracer Python+Rust and CMFD
  widths handle it; `from_core` rejects gaps) + `vera-p2-pin-powers.toml`
  (KENO octants). `pincell/assembly.py` (pin powers vs KENO). `transport/cp.py`:
  1D cylindrical CP (Bickley Ki3, Carlvik chords, white boundary), `pin_cell`
  cross-check (+71 pcm vs MOC on 1A). `shielding.coupled_collisions`
  (P_j->i on a Sigma_fuel x Sigma_clad grid) + `coupled_flux` (rings + clad
  marched together; march takes per-region totals/scatter, Rust too);
  `run(coupled_clad=)` / `--clad`. `moc.surface_flux` (boundary psi ->
  surface scalar flux) -> ADFs in `fewgroup` (`Solution.psi`). CMFD waits
  for positive coarse fluxes (negative transport-corrected H self-scatter
  makes early guide-tube fluxes negative). Assembly needs 64 azim / 0.03 cm
  (0.05 cm tracks: -127 pcm, thinner than the gap).
- `src/fission/fewgroup/` (C3): `homogenize.py` (flux-volume weighting of a
  heterogeneous Result, plain totals + P0/P1, chi by region production),
  `b1.py` (B1 equations, alpha(B/Sigma_t), secant search for k=1, D = J/(B phi)),
  `condense.py` (`FewGroup`, current-weighted D, `material()` = kwargs of
  nodal.Material, JSON writer), `kappa.py` (MF1/MT458 recoverable energy ER;
  kept outside fission/endf so the table cache stays valid). CLI
  `python -m fission.fewgroup`. `reference.few_group_reference` condenses
  OpenMC tallies.
- `src/fission/branches/` (M2): `[branches.<name>]` in the model file
  (`model.spec.Branches`: base state + values per axis of `BRANCH_AXES`;
  points = tensor product, C order, names `<name>-i-j-k-l`; loader checks
  increasing values containing the base value). `run_matrix` (solver
  infinite | moc, B1 + condense per point, table rewritten after each point
  and resumed if settings match, shaped libraries of a point deleted
  after it), `Table.array(key)` (grid-shaped), `coefficients` (d rho/dx at
  the base). `fewgroup.homogenize.from_mixture` = constants of the infinite
  medium (chi weighted by nuclide production; keeps k to 1e-15).
- `src/fission/tables/` (M3): `multilinear` (axis by axis, end segment
  extrapolates, single-node axis constant), `XSTable` (from a branch Table,
  fuel axis in sqrt(T) via `TRANSFORMS`, `without(axis, i)` for leave one
  out, JSON write/read), `select` (greedy terms by PRESS e/(1-h), fewest
  within 1% of best, one term set for all constants), `fit`, `k_infinity`.
  `branches.sample_states` / `run_states` = off-grid check points.
- `src/fission/pincell/`: `infinite.py` (homogeneous mixture, eigen solve,
  four factors), `run.py` (driver on a model lattice + state; `UNSUPPORTED`
  nuclides; NR self-shielding via a flux shape 1/Sigma_t(E)).
- `src/fission/data.py`: ENDF/B-VIII.1 downloader (NNDC tarballs, SHA-1 checked)
  → `endf-data/` (gitignored).
- `src/fission/endf/`: ENDF-6 reader (`records`, `interpolation`,
  `evaluation`, `resonances` for MF2), resonance physics (`channel`, `slbw`,
  `mlbw`, `reich_moore`), tables (`reconstruct` adaptive bisection,
  `pointwise` full-range tables, `cache` disk cache), `fast.py` = thin
  wrappers over the Rust kernels.
- `src/fission/endf/distributions/`: secondary neutrons: `nu` (MF1
  452/455/456, delayed groups, MF6/MT18 P(nu)), `angular` (MF4), `energy`
  (MF5 LF=1,5), `mf6` (neutron products, LAW 0/<0/1/2), `kinematics`
  (two-body CM->LAB, Kalbach-Mann), `incident` (incl. unit-base
  interpolation), `sab` (MF7 S(alpha,beta) reader). A subpackage, so it is
  outside the 0 K table-cache fingerprint (but inside the processing one).
- `rust/src/` Rust ports of A3–A5: `groupr.rs` (hat integrals), `gauss.rs`
  (Gauss–Legendre of any order), `transfer.rs` (two-body kernels + tables,
  threaded), `thermal.rs` (free gas + S(alpha,beta)), `infinite.rs` (LU +
  power iteration). Still Python: continuum (MF6 LAW=1) and fission kernels.
- `rust/src/`: PyO3 module `fission._rust` (stubs in `src/fission/_rust.pyi`).
  Rust structs read Python objects by attribute via `FromPyObject`, so the
  Python dataclasses are the interface.
- `tests/`: one file per module plus `test_rust_parity.py`;
  `helpers.endf_line` builds synthetic 80-column records.
- `docs/theory/NN-*.md`: theory chapters (Markdown + LaTeX, ENDF-102
  notation); `22-references.md` holds anchors `[]{#ref-xxx}`; cite as
  `[[Name](#ref-xxx) §x.y]`.

## Conventions and gotchas

- Units follow ENDF: eV, barns, radii in 1e-12 cm (so π/k² is in barns).
  Physical constants are CODATA 2018 (`channel.K_CONSTANT`).
- Reconstruction defaults: tol = 1e-3 relative, atol = 1e-5 b. Total MT1 is
  rebuilt as the fsum of non-redundant partials (MT5 is a partial, not
  redundant); jumps are shaded by 1e-7 and must not be bisected.
- Reich–Moore capture is computed directly (Σ Γγ|A|²), never by subtraction
  (cancellation at 1e-7).
- The table cache key includes a fingerprint of every `fission/endf/*.py` and
  the `_rust` .so: **any edit there invalidates all cached tables**, and the
  next `make tables` / `make figures` rebuilds them (~15 s; U-238 ~10 s).
  Any Rust change does the same (the .so is fingerprinted). Broadened tables
  are keyed on the 0 K key + T + params + `fission/processing` code.
- Doppler: free-gas exact kernel, weights shared by all MTs; Gauss rule
  chosen by panel width; seeds at ±1–3 Γ_D around sparse 0 K points (the
  midpoint test alone misses rounded corners whose shifts cancel); refined
  then thinned; observed error ≤ ~1.5 tol.
- Group averaging: group 0 = highest energy; S[l, g, g'] = transfer g -> g'
  (LAB P_l). All constants are Int_g sigma w K with K tabulated on a coarse
  kernel grid; only the hat integrals A[g, j] see the fine grid (Rust
  `groupr.rs`). Libraries are infinitely dilute (self-shielding is C1).
- Flux shapes: `group_library(flux_shape=(E, F))` weights every integral by
  w F (table pre-multiplied with `shape_table`); `cached_groups(...,
  flux_shape, shape_name)` caches them. Used for homogeneous NR shielding.
- Unsupported formats so far: LRF=7 (Fe-54, Fe-57), URR LRU=2 LRF=1 (Fe-58),
  Reich-Moore negative-J channel spins (Cr-53); left out of VERA clad.
- Thermal (A4): H-1 (awr < 1) is a free gas from the end of its law up to
  1 keV (OpenMC's rule). Otherwise below 4 eV the elastic rows are replaced by a free-gas kernel
  (default, normalized to the broadened elastic table, as NJOY) or by MF7
  S(alpha,beta) (`thermal="HinH2O"`, absolute). Kernel tabulated every 0.02
  in lethargy (detailed balance ~0.1%); cosine integral in s = sqrt(alpha).
  LAT=1 rescales only the table lookup; NJOY's SCT fallback rules copied.
  H in H2O exists only at the file's temperatures (no 1200 K).
- Rust refinement is threaded with std only (`adaptive_grid_parallel`,
  identical to serial). `pointwise(use_rust=True)` runs the whole full-range
  refinement in Rust (`rust/src/pointwise.rs`, exact `fsum.rs` port of
  `math.fsum`); `use_rust=False` is the all-Python reference.
- ENDF/B-VIII.1 data quirks found so far: the nu_p implied by MF6/MT18 P(nu)
  differs from MF1/MT456 by up to 1.2% (U-235) / 1.9% (U-238), so use only
  normalized P(nu) with nu-bar from MF1; U-238 chi at 20 MeV integrates to
  1.0003 (renormalize spectra in group averaging).
- Physics numbers in docs or tests must be computed by this code or taken
  from a cited source, not from memory. Known checks: U-238 thermal capture
  2.682 b, U-235 fission 586.2 b, U-238 capture RI 275.28 b.
- Docs: verify external links and section numbers before citing; prefer free
  sources and mark (free)/(not free). Python 3.14+ syntax is fine (`type`
  aliases etc.).

## Roadmap

Done: setup, ENDF reader, resonance reconstruction (0 K pointwise tables),
A1 Doppler broadening, A2 secondary-neutron data, A3 multigroup libraries,
A4 thermal scattering, A5 infinite-medium k-inf, M1 model file, B1 geometry + ray tracing, B2 MOC (C5G7), B3 CMFD, C1 self-shielding, C2 OpenMC reference + SPH + rings, C3 few-group constants (B1), C4 CP + clad + 17x17 assembly; M2 state matrix, M3 tables; theory docs ch. 0–22. Data and transport are independent tracks that meet at
self-shielding (C1).

- **A: nuclear data → multigroup library**
  - A1 Doppler broadening (done)
  - A2 rest of the reader: ν̄, χ, MF4/MF5/MF6 (done)
  - A3 group averaging: σ_g, νσ_f, χ, P0–P3 matrices (done)
  - A4 thermal scattering: free gas + S(α,β) for H in H2O (done)
  - A5 infinite-medium spectrum + k∞ of the homogenized VERA 1 cell (done)
- **B: transport**
  - B1 pin-cell geometry + cyclic ray tracing (done)
  - B2 MOC, verified on C5G7 7-group data (done; Rust sweep included)
  - B3 CMFD acceleration (done: C5G7 1242 -> ~33 iterations)  
- **C: coupling**
  - C1 self-shielding (done: equivalence, exact escape probabilities,
    ultra-fine slowing down; VERA 1A -713 pcm, Doppler trends within ~10%)
  - C2 real pin cell vs OpenMC CE on ENDF/B-VIII.1 (done: data vs KENO
    -24..+69 pcm; SPH + coupled pellet rings take 1A from -689 to -162 pcm vs
    OpenMC; 1D -348, Doppler steps 8-12% too strong; rest = U-238 capture in
    large-Gamma_n resonances, clad NR treatment a candidate)
  - C3 homogenize/condense to 2 groups + B1 → nodal (done: flux-volume
    homogenization keeps the MOC k to 4e-10; B1 critical B^2 = 3.26e-3/cm^2
    for 1A, M^2 57 cm^2; two-group constants run unchanged in branch
    `claude-nodal` via `make nodal-check`; vs OpenMC 2-group XS within 0.5%,
    thermal fuel/water flux ratio +0.4%)
  - C4 CP cross-check, clad in the slowing down, 17×17 assembly, ADFs
    (done: CP +71 pcm vs MOC; coupled clad +85..115 pcm, pin 1A-1D vs
    OpenMC -77/-62/-139/-234; VERA 2A-2D vs OpenMC -11/-46/-80/-117, pin
    powers vs KENO 0.05-0.07% RMS; ADF 2A 0.996/0.9995)
  - Doppler trend (done, ch. 19 "Moving targets"): data = NJOY to 0.01%;
    the march is target-at-rest. OpenMC 1B->1D: default (constant-sigma
    free gas < 400 kT) -1833, at rest > 4 eV -1970, exact RVS -2101; ours
    -2005. Vs RVS: 1B +17, 1D +113 pcm (step 4.6% too weak). Reference
    `spec --free-gas-cutoff 4 | --resonance-scattering --tag=-rvs`;
    `compare --result <file>`.
  - Moving targets (done, ch. 19): `pincell/moving.py` free-gas kernel with
    0 K sigma (Ouisloumen-Sanchez; split at z* into prefix/suffix log-sums,
    hat weights W_j of z sigma_0; mesh uniform in x = sqrt(AE/kT), dx 1/8),
    Rust `moving.rs` (3 s, parity 2e-15); `march(moving=[kernel|None])`
    Gauss-Seidel passes (first at rest), flux continued flat below the
    range; `run(moving_targets=)` / `--moving` (U-238, 1 eV-1 keV). Vs OpenMC
    RVS-U238 + others at rest (tag -rvsrest): step 1B->1D -2208 vs -2172+-28;
    level 1B -124, 1D -160 at default settings.
    Kinematics seam at the kernel top is physical (OpenMC has it too).
  - The 66 pcm (done, ch. 19 "Chasing the level"): OpenMC keff noise (use
    tally k = nuSigf/abs for differences) + discretization. Default MOC
    tracking 64/0.02/3 GL polar is ~190 pcm low on the 172-group pin (GL
    overestimates the thick-region inflow integral, int sin dmu, 0.46% at 3
    points; unshielded +560 pcm); 3 pellet rings over-capture in the
    large-Gamma_n resonances (20.9, 36.7 eV), which moving targets amplify.
    At 10 rings + final MOC 256/0.005/12: moving effect 1B -86 / 1D -277 vs
    OpenMC tally -70 / -267; capture added +63/+209 vs +50+-8/+201+-8.
    Converged level at rest is +170 pcm vs OpenMC, all from U-238 capture at
    4-10 eV (6.67 eV) being 2.6% low (hidden before by the polar error).
    Defaults unchanged so far.
  - The 6.67 eV deficit (done, ch. 19 "hydrogen moves too", ch. 11): the
    flux was low by ~2.7 kT/E in every region from H-1 at rest above 4 eV.
    OpenMC: targets at rest only if awr > 1, so H is S(a,b) to 10 eV and
    free gas above. Library rule now: law up to `law.emax` (water 10 eV ->
    XMAS rows to 9.906 eV), free gas for awr < 1 up to `LIGHT_FREE_GAS_TOP`
    = 1 keV, at rest above. Converged (10 rings, final MOC 256/0.005/12):
    1B -44, 1D -42 pcm vs OpenMC rest; every U-238 capture band within 1%.
    Default settings (3 rings, 64/0.02/3) vs OpenMC default runs now
    1A-1D -271/-274/-351/-447 (were -77/-62/-139/-234 by cancellation).
    Numbers in ch. 12-19 and the assembly predate the fix (not rerun).
  - Converged defaults (done, ch. 15 + 19): Tabuchi-Yamamoto polar set
    (`Polar.tabuchi_yamamoto`, 1-3 angles, OpenMOC values; default of
    `tracking(quadrature=)`; Ki3 error 3.4e-5 at 3 angles < GL 12). Pin
    defaults `heterogeneous.AZIM/SPACING/POLAR/FUEL_RINGS` = 256 / 0.02 / 3
    / 10 (~3.5 min warm, ~20 pcm from converged); CLIs take --fuel-rings.
    Assembly keeps 64 / 0.03 / 3 rings. `fixed_source` stops on the region
    tallies, drops face currents, takes `boundary=` (SPH warm start).
    1A-1D vs OpenMC default -62/-92/-85/-164; 1B/1D vs at-rest -89/-23, all
    rates within 1e-4 per fission neutron. C5G7 +12 pcm.
  - Next: M4 (needs a new nodal code), the MOC state matrix (~3 h), the
    march with moving H, or the assembly at converged tracking
- **M: model building** (describe the reactor once, as data; generate every calculation from it)
  - M1 declarative model file: TOML layers + states, number densities, homogenization (done)
  - M2 state-matrix driver (done, ch. 20): `make branches`; VERA matrix
    `[branches.vera]` around 1B, 4 x 2 x 2 x 3 = 48 points, infinite solver
    ~90 s/point (shaped libraries). Coefficients at 1B: fuel T -3.96 pcm/K,
    moderator T -1.73 pcm/K, density +22.5 pcm per 0.01 g/cm3, boron
    -5.94 pcm/ppm; density effect +1092 / +185 / -293 pcm at 0 / 1300 /
    2000 ppm. Homogeneous Doppler 1B->1D -2201 vs KENO -1422 (no spatial
    shielding). MOC run of the matrix (~3 h) not done yet.
  - M3 cross-section tables (done, ch. 20): multilinear in sqrt(Tf) (LOO
    at 900 K: 13 pcm vs 94 in T; upscatter wants T), PRESS-selected
    polynomial (13 terms: rho, Tf, B, Tm, rho B, ...). 10 check points:
    multilinear sqrt(T) max +29 pcm (boron curvature 0-1300 ppm), T +57,
    fit +4 pcm; constants within 0.35%. ADFs interpolated like XS (1 for
    the infinite medium).
  - M4 core model -> nodal input (core map + tables); the nodal code must be
    rewritten here (branch `claude-nodal` deleted)
  - Later: hex geometry, eighth symmetry, burnup histories (with depletion)
- Deferred: unresolved-range probability tables, depletion.
