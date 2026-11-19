# fission

A lattice-physics code written from scratch as a **learning project**: raw
ENDF/B-VIII.1 → pointwise cross sections → 2D pin-cell transport (MOC, later
CP) → few-group constants → 3D nodal diffusion. No NJOY/OpenMC in the
pipeline (they are installed under tools/ for independent checks: `make tools`). Python package:
`fission`. An earlier NEM/CMFD nodal diffusion code lives on branch
`claude-nodal`.

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
make pin STATES="1A 1B" PIN_SHIELD=sd PIN_ARGS="--rings --sph"  # heterogeneous VERA pin, MOC + shielded ENDF libraries
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
  notation); `19-references.md` holds anchors `[]{#ref-xxx}`; cite as
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
- Thermal (A4): below 4 eV the elastic rows are replaced by a free-gas kernel
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
A4 thermal scattering, A5 infinite-medium k-inf, M1 model file, B1 geometry + ray tracing, B2 MOC (C5G7), B3 CMFD, C1 self-shielding, C2 OpenMC reference + SPH + rings; theory docs ch. 0–19. Data and transport are independent tracks that meet at
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
  - **C3 homogenize/condense to 2 groups + B1 → nodal (next; branch `claude-nodal`)**
  - C4 CP cross-check (also pointwise fuel+gap+clad slowing down), 17×17 assembly
- **M: model building** (describe the reactor once, as data; generate every calculation from it)
  - M1 declarative model file: TOML layers + states, number densities, homogenization (done)
  - M2 state-matrix driver: lattice calculation at every branch point
    (fuel T, moderator density/T, boron), first with the infinite-medium
    solver, later MOC
  - M3 cross-section tables for the nodal code: multilinear interpolation,
    ADFs from C3; then a fitted/sparse parametrization compared with it
  - M4 core model -> nodal input (core map + tables), code on `claude-nodal`
  - Later: hex geometry, eighth symmetry, burnup histories (with depletion)
- Deferred: unresolved-range probability tables, depletion.
