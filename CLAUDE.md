# fission

A lattice-physics code written from scratch as a **learning project**: raw
ENDF/B-VIII.1 → pointwise cross sections → 2D pin-cell transport (MOC, later
CP) → few-group constants → 3D nodal diffusion. No NJOY/OpenMC in the
pipeline (they may be used as independent checks). Python package:
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
make check        # lint + format-check + mypy + pytest + cargo test (run before reporting done)
make test ARGS="-k reich -x"   # pytest; slow tests are deselected by default
make test-slow    # full-resolution physics checks (minutes)
make fix          # ruff --fix + format
make data         # download pin-cell evaluations into endf-data/
make tables       # build/load cached pointwise tables for pin-cell nuclides
make doppler-tables TEMPS="293.6 900"  # broadened tables (endf-data/doppler/)
make docs         # docs/build/theory.html (pandoc, MathML, offline)
make figures      # regenerate docs/figures/*.svg (deterministic output)
make show NUC=U235 | make scan NUC=U238 E=6.674 | make rm | make doppler | make secondary   # quick CLIs
```

- Everything runs through `uv run`. Editing `rust/**/*.rs` or `Cargo.toml`
  makes the next `uv run` rebuild the extension (`[tool.uv] cache-keys`); no
  manual `maturin develop` needed.
- Lint/type scope is `src tests docs/figures`.

## Layout

- `src/fission/processing/`: stages that take 0 K tables as input
  (`doppler.py`; later group averaging). Kept outside `endf/` so edits do
  not invalidate the table cache.
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
  interpolation). A subpackage, so it is outside the table-cache fingerprint.
- `rust/src/`: PyO3 module `fission._rust` (stubs in `src/fission/_rust.pyi`).
  Rust structs read Python objects by attribute via `FromPyObject`, so the
  Python dataclasses are the interface.
- `tests/`: one file per module plus `test_rust_parity.py`;
  `helpers.endf_line` builds synthetic 80-column records.
- `docs/theory/NN-*.md`: theory chapters (Markdown + LaTeX, ENDF-102
  notation); `11-references.md` holds anchors `[]{#ref-xxx}`; cite as
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
A1 Doppler broadening, A2 secondary-neutron data; theory docs ch. 0–11. Data and transport are independent tracks that meet at
self-shielding (C1).

- **A: nuclear data → multigroup library**
  - A1 Doppler broadening (done)
  - A2 rest of the reader: ν̄, χ, MF4/MF5/MF6 (done)
  - **A3 group averaging (next)** (σ_g, νσ_f, χ, P0/P1 matrices)
  - A4 thermal scattering (free gas, then S(α,β))
  - A5 infinite-medium spectrum + k∞ of a homogenized pin cell
- **B: transport**
  - B1 pin-cell geometry + cyclic ray tracing
  - B2 MOC, verified on C5G7 7-group data
  - B3 Rust sweep
- **C: coupling**
  - C1 self-shielding (equivalence/Bondarenko + Dancoff, or subgroup)
  - C2 real pin cell vs VERA problem 1
  - C3 homogenize/condense to 2 groups + B1 → nodal (branch `claude-nodal`)
  - C4 CP cross-check, 17×17 assembly
- Deferred: unresolved-range probability tables, depletion.
