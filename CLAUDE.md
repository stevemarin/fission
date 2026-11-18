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
make docs         # docs/build/theory.html (pandoc, MathML, offline)
make figures      # regenerate docs/figures/*.svg (deterministic output)
make show NUC=U235 | make scan NUC=U238 E=6.674 | make rm ...   # quick CLIs
```

- Everything runs through `uv run`. Editing `rust/**/*.rs` or `Cargo.toml`
  makes the next `uv run` rebuild the extension (`[tool.uv] cache-keys`); no
  manual `maturin develop` needed.
- Lint/type scope is `src tests docs/figures`.

## Layout

- `src/fission/data.py`: ENDF/B-VIII.1 downloader (NNDC tarballs, SHA-1 checked)
  → `endf-data/` (gitignored).
- `src/fission/endf/`: ENDF-6 reader (`records`, `interpolation`,
  `evaluation`, `resonances` for MF2), resonance physics (`channel`, `slbw`,
  `mlbw`, `reich_moore`), tables (`reconstruct` adaptive bisection,
  `pointwise` full-range tables, `cache` disk cache), `fast.py` = thin
  wrappers over the Rust kernels.
- `rust/src/`: PyO3 module `fission._rust` (stubs in `src/fission/_rust.pyi`).
  Rust structs read Python objects by attribute via `FromPyObject`, so the
  Python dataclasses are the interface.
- `tests/`: one file per module plus `test_rust_parity.py`;
  `helpers.endf_line` builds synthetic 80-column records.
- `docs/theory/NN-*.md`: theory chapters (Markdown + LaTeX, ENDF-102
  notation); `09-references.md` holds anchors `[]{#ref-xxx}`; cite as
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
  next `make tables` / `make figures` rebuilds them (~2–3 min, U-238 ~90 s).
- Physics numbers in docs or tests must be computed by this code or taken
  from a cited source, not from memory. Known checks: U-238 thermal capture
  2.682 b, U-235 fission 586.2 b, U-238 capture RI 275.28 b.
- Docs: verify external links and section numbers before citing; prefer free
  sources and mark (free)/(not free). Python 3.14+ syntax is fine (`type`
  aliases etc.).

## Roadmap

0 setup → 1 ENDF reader → 2 resonance reconstruction (0 K) → **3 pin-cell
geometry + cyclic ray tracing (next)** → 4 MOC → 5 Doppler broadening →
6 group constants (+ S(α,β)) → 7 first pin cell (dilute) → 8 self-shielding,
compare with VERA problem 1 → 9 homogenize/condense, B1, feed nodal →
10 CP cross-check, 17×17 assembly.
