//! Compute kernels for `fission`, exposed to Python as `fission._rust`.
//!
//! Every kernel here is a port of a Python reference implementation in
//! `src/fission`; tests/test_rust_parity.py checks that the two agree.

pub mod channel;
pub mod cmfd;
pub mod complex;
pub mod diffusion;
pub mod doppler;
pub mod fsum;
pub mod fuel;
pub mod gauss;
pub mod groupr;
pub mod infinite;
pub mod mlbw;
pub mod moving;
pub mod moc;
pub mod pointwise;
pub mod raytrace;
pub mod reconstruct;
pub mod reich_moore;
pub mod slowing;
pub mod spin_groups;
pub mod tab1;
pub mod tables;
pub mod thermal;
pub mod transfer;
pub mod water;

use numpy::{IntoPyArray, PyArray1, PyReadonlyArray1};
use pyo3::prelude::*;
use spin_groups::{PotentialTerm, SpinGroupIn};

/// A resolved-range formalism: (elastic, capture, fission) in barns at an energy.
pub trait Resonances: Sync {
    fn eval(&self, energy: f64) -> [f64; 3];
}

impl Resonances for reich_moore::ReichMoore {
    fn eval(&self, energy: f64) -> [f64; 3] {
        reich_moore::ReichMoore::eval(self, energy)
    }
}

impl Resonances for mlbw::Mlbw {
    fn eval(&self, energy: f64) -> [f64; 3] {
        mlbw::Mlbw::eval(self, energy)
    }
}

type Array<'py> = Bound<'py, PyArray1<f64>>;
type IntArray<'py> = Bound<'py, PyArray1<i64>>;
/// Face crossings from Python: (start, position, face, sign, number of faces).
type FacesIn<'py> = (
    PyReadonlyArray1<'py, i64>,
    PyReadonlyArray1<'py, i64>,
    PyReadonlyArray1<'py, i64>,
    PyReadonlyArray1<'py, f64>,
    usize,
);

fn evaluate_many<'py, R: Resonances>(
    model: &R,
    py: Python<'py>,
    energies: PyReadonlyArray1<'py, f64>,
) -> PyResult<(Array<'py>, Array<'py>, Array<'py>)> {
    let energies = energies.as_slice()?;
    let mut out = [
        Vec::with_capacity(energies.len()),
        Vec::with_capacity(energies.len()),
        Vec::with_capacity(energies.len()),
    ];
    for &e in energies {
        let v = model.eval(e);
        for j in 0..3 {
            out[j].push(v[j]);
        }
    }
    let [a, b, c] = out;
    Ok((a.into_pyarray(py), b.into_pyarray(py), c.into_pyarray(py)))
}

/// Reich-Moore resolved-range cross sections (port of `fission.endf.reich_moore`).
#[pyclass(module = "fission._rust", name = "ReichMoore", frozen)]
struct PyReichMoore {
    inner: reich_moore::ReichMoore,
}

#[pymethods]
impl PyReichMoore {
    /// Build from the Python `spin_groups(rr)`, `potential_terms(rr)` and target spin.
    #[new]
    fn new(spin_groups: Vec<SpinGroupIn>, potentials: Vec<PotentialTerm>, spi: f64) -> Self {
        PyReichMoore {
            inner: reich_moore::ReichMoore::new(spin_groups, potentials, spi),
        }
    }

    /// (elastic, capture, fission) in barns at one energy.
    fn __call__(&self, energy: f64) -> (f64, f64, f64) {
        let [e, c, f] = self.inner.eval(energy);
        (e, c, f)
    }

    /// Vectorized evaluation: three arrays (elastic, capture, fission).
    fn evaluate<'py>(
        &self,
        py: Python<'py>,
        energies: PyReadonlyArray1<'py, f64>,
    ) -> PyResult<(Array<'py>, Array<'py>, Array<'py>)> {
        evaluate_many(&self.inner, py, energies)
    }
}

/// Multi-Level Breit-Wigner cross sections (port of `fission.endf.mlbw`).
#[pyclass(module = "fission._rust", name = "Mlbw", frozen)]
struct PyMlbw {
    inner: mlbw::Mlbw,
}

#[pymethods]
impl PyMlbw {
    /// Build from the Python `spin_groups(rr)`, `potential_terms(rr)` and target spin.
    #[new]
    fn new(spin_groups: Vec<SpinGroupIn>, potentials: Vec<PotentialTerm>, spi: f64) -> Self {
        PyMlbw {
            inner: mlbw::Mlbw::new(spin_groups, potentials, spi),
        }
    }

    /// (elastic, capture, fission) in barns at one energy.
    fn __call__(&self, energy: f64) -> (f64, f64, f64) {
        let [e, c, f] = self.inner.eval(energy);
        (e, c, f)
    }

    /// Vectorized evaluation: three arrays (elastic, capture, fission).
    fn evaluate<'py>(
        &self,
        py: Python<'py>,
        energies: PyReadonlyArray1<'py, f64>,
    ) -> PyResult<(Array<'py>, Array<'py>, Array<'py>)> {
        evaluate_many(&self.inner, py, energies)
    }
}

/// Either formalism, as accepted by `reconstruct`.
#[derive(FromPyObject)]
enum Model<'py> {
    ReichMoore(PyRef<'py, PyReichMoore>),
    Mlbw(PyRef<'py, PyMlbw>),
}

impl Model<'_> {
    fn resonances(&self) -> &dyn Resonances {
        match self {
            Model::ReichMoore(m) => &m.inner,
            Model::Mlbw(m) => &m.inner,
        }
    }
}

/// Adaptive reconstruction of resonance + background cross sections
/// (port of `fission.endf.reconstruct.adaptive_grid` with `resolved_range_xs`).
///
/// `model` is a ReichMoore or Mlbw object; `backgrounds` holds three optional
/// Python `Tabulated1D` objects (elastic, capture, fission MF3 sections).
#[pyfunction(name = "reconstruct")]
#[pyo3(signature = (model, backgrounds, seeds, tol=1e-3, atol=1e-5, min_rel_width=1e-9))]
#[allow(clippy::type_complexity)]
fn py_reconstruct<'py>(
    py: Python<'py>,
    model: Model<'py>,
    backgrounds: Vec<Option<tab1::Tab1>>,
    seeds: Vec<f64>,
    tol: f64,
    atol: f64,
    min_rel_width: f64,
) -> PyResult<(Array<'py>, Array<'py>, Array<'py>, Array<'py>)> {
    if backgrounds.len() != 3 || seeds.len() < 2 {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "need three background slots and at least two seeds",
        ));
    }
    let resonances = model.resonances();
    let f = |e: f64| {
        let mut v = resonances.eval(e);
        for (j, bg) in backgrounds.iter().enumerate() {
            if let Some(t) = bg {
                v[j] += t.eval(e);
            }
        }
        v
    };
    // the physics does not touch Python objects: let other Python threads run
    let (energies, values) = py.detach(|| {
        let none = reconstruct::Unsplittable::new();
        reconstruct::adaptive_grid_parallel(f, &seeds, tol, atol, min_rel_width, &none)
    });
    let pick = |j: usize| values.iter().map(|v| v[j]).collect::<Vec<f64>>();
    Ok((
        energies.into_pyarray(py),
        pick(0).into_pyarray(py),
        pick(1).into_pyarray(py),
        pick(2).into_pyarray(py),
    ))
}

// ---------------------------------------------------- full-range tables
/// Full-range pointwise table (port of the refinement in `pointwise()`).
///
/// `tables[i]` is the MF3 `Tabulated1D` of value slot i (None for the total
/// or a missing reaction); `resonance_slots` are the slots of (elastic,
/// capture, fission); slot `total` is the exact sum of `partials`.
/// Returns the grid and one array per slot.
#[pyfunction]
#[pyo3(signature = (model, el, eh, resonance_slots, tables, total, partials, seeds,
                    unsplittable, tol=1e-3, atol=1e-5, min_rel_width=1e-9))]
#[allow(clippy::too_many_arguments)]
fn pointwise_table<'py>(
    py: Python<'py>,
    model: Option<Model<'py>>,
    el: f64,
    eh: f64,
    resonance_slots: [Option<usize>; 3],
    tables: Vec<Option<tab1::Tab1>>,
    total: usize,
    partials: Vec<usize>,
    seeds: Vec<f64>,
    unsplittable: Vec<(f64, f64)>,
    tol: f64,
    atol: f64,
    min_rel_width: f64,
) -> PyResult<(Array<'py>, Vec<Array<'py>>)> {
    let n = tables.len();
    let in_range = |s: &usize| *s < n;
    if seeds.len() < 2
        || !in_range(&total)
        || !partials.iter().all(in_range)
        || !resonance_slots.iter().flatten().all(in_range)
    {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "bad slots or fewer than two seeds",
        ));
    }
    let reactions = pointwise::Reactions {
        tables,
        resonance: model.as_ref().map(|m| pointwise::ResonancePart {
            model: m.resonances(),
            el,
            eh,
            slots: resonance_slots,
        }),
        total,
        partials,
    };
    let keep: reconstruct::Unsplittable = unsplittable
        .iter()
        .map(|(a, b)| (a.to_bits(), b.to_bits()))
        .collect();
    let (grid, rows) = py.detach(|| {
        reconstruct::adaptive_grid_parallel(
            |e| reactions.eval(e),
            &seeds,
            tol,
            atol,
            min_rel_width,
            &keep,
        )
    });
    Ok((grid.into_pyarray(py), columns_to_arrays(py, n, &rows)))
}

/// Group-averaging hat integrals (port of `processing.groupr.hat_integrals`).
/// Returns (groups, hat indices, values per column) of the nonzero entries.
#[pyfunction]
#[allow(clippy::type_complexity)]
fn group_hat_integrals(
    py: Python<'_>,
    x: Vec<f64>,
    columns: Vec<Vec<f64>>,
    weight: Vec<f64>,
    boundaries: Vec<f64>,
    z: Vec<f64>,
    max_ratio: f64,
) -> PyResult<(Vec<usize>, Vec<usize>, Vec<Vec<f64>>)> {
    let w = groupr::Weight::from_params(&weight)
        .ok_or_else(|| pyo3::exceptions::PyValueError::new_err("bad weight parameters"))?;
    if x.len() < 2
        || z.len() < 2
        || boundaries.len() < 2
        || columns.iter().any(|c| c.len() != x.len())
    {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "bad grids or columns",
        ));
    }
    let acc = py.detach(|| groupr::hat_integrals(&x, &columns, w, &boundaries, &z, max_ratio));
    let mut gs = Vec::with_capacity(acc.len());
    let mut js = Vec::with_capacity(acc.len());
    let mut vs = Vec::with_capacity(acc.len());
    for ((g, j), v) in acc {
        gs.push(g);
        js.push(j);
        vs.push(v);
    }
    Ok((gs, js, vs))
}

// ---------------------------------------------------- two-body transfer kernels
type PanelIn = (i64, Vec<f64>, Vec<f64>);

fn angular_from(
    energies: Vec<f64>,
    breaks: Vec<usize>,
    laws: Vec<i64>,
    panels: Vec<PanelIn>,
) -> PyResult<transfer::Angular> {
    if energies.is_empty() {
        return Ok(transfer::Angular::isotropic());
    }
    let panels = panels
        .into_iter()
        .map(|(kind, a, b)| {
            if kind == 0 {
                transfer::Panel::Legendre(a)
            } else {
                transfer::Panel::Tabulated(a, b)
            }
        })
        .collect();
    transfer::Angular::new(energies, breaks, laws, panels)
        .map_err(pyo3::exceptions::PyValueError::new_err)
}

/// K_l(E -> g') of a two-body reaction, flattened [l][g'] (port of `two_body_kernel`).
#[pyfunction]
#[allow(clippy::too_many_arguments)]
fn two_body_kernel<'py>(
    py: Python<'py>,
    energies: Vec<f64>,
    breaks: Vec<usize>,
    laws: Vec<i64>,
    panels: Vec<PanelIn>,
    awr: f64,
    q: f64,
    boundaries: Vec<f64>,
    lmax: usize,
    energy: f64,
) -> PyResult<Array<'py>> {
    let ang = angular_from(energies, breaks, laws, panels)?;
    Ok(transfer::two_body_kernel(&ang, awr, q, &boundaries, lmax, energy).into_pyarray(py))
}

/// Tabulated two-body kernel: (energies, rows) (port of `build_kernel` + `two_body_kernel`).
#[pyfunction]
#[pyo3(signature = (energies, breaks, laws, panels, awr, q, boundaries, lmax, seeds, sigma=None, adaptive=true, tol=1e-3, atol=1e-5))]
#[allow(clippy::too_many_arguments)]
fn two_body_table<'py>(
    py: Python<'py>,
    energies: Vec<f64>,
    breaks: Vec<usize>,
    laws: Vec<i64>,
    panels: Vec<PanelIn>,
    awr: f64,
    q: f64,
    boundaries: Vec<f64>,
    lmax: usize,
    seeds: Vec<f64>,
    sigma: Option<tab1::Tab1>,
    adaptive: bool,
    tol: f64,
    atol: f64,
) -> PyResult<(Array<'py>, Vec<Array<'py>>)> {
    let ang = angular_from(energies, breaks, laws, panels)?;
    if seeds.len() < 2 {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "need at least two seeds",
        ));
    }
    let (z, rows) = py.detach(|| {
        transfer::two_body_table(
            &ang,
            awr,
            q,
            &boundaries,
            lmax,
            &seeds,
            sigma.as_ref(),
            adaptive,
            tol,
            atol,
        )
    });
    Ok((
        z.into_pyarray(py),
        rows.into_iter().map(|r| r.into_pyarray(py)).collect(),
    ))
}

// ---------------------------------------------------- thermal kernels
/// A thermal model from Python: kind 0 = free gas [awr, kT, sigma_b];
/// kind 1 = tabulated [kT, A0, sigma_b, scale, T_eff/T] + alpha, beta, ln S.
fn thermal_model(
    kind: i64,
    params: Vec<f64>,
    alpha: Vec<f64>,
    beta: Vec<f64>,
    ln_s: Vec<f64>,
) -> PyResult<thermal::Model> {
    let bad = || pyo3::exceptions::PyValueError::new_err("bad thermal model parameters");
    match kind {
        0 if params.len() == 3 => Ok(thermal::Model::FreeGas {
            awr: params[0],
            kt: params[1],
            sigma_bound: params[2],
        }),
        1 if params.len() == 5
            && alpha.len() >= 2
            && beta.len() >= 2
            && ln_s.len() == alpha.len() * beta.len() =>
        {
            Ok(thermal::Model::Tabulated {
                kt: params[0],
                a0: params[1],
                sigma_bound: params[2],
                scale: params[3],
                t_ratio: params[4],
                alpha,
                beta,
                ln_s,
            })
        }
        _ => Err(bad()),
    }
}

/// Thermal kernels at many incident energies: per energy the flattened
/// K_l(E -> g') followed by the total (port of `thermal_kernel`).
#[pyfunction]
#[allow(clippy::too_many_arguments)]
fn thermal_table<'py>(
    py: Python<'py>,
    kind: i64,
    params: Vec<f64>,
    alpha: Vec<f64>,
    beta: Vec<f64>,
    ln_s: Vec<f64>,
    boundaries: Vec<f64>,
    lmax: usize,
    energies: Vec<f64>,
) -> PyResult<Vec<Array<'py>>> {
    let model = thermal_model(kind, params, alpha, beta, ln_s)?;
    let rows = py.detach(|| {
        transfer::parallel_map(&energies, |e| {
            let (mut k, total) = thermal::thermal_kernel(&model, &boundaries, lmax, e);
            k.push(total);
            k
        })
    });
    Ok(rows.into_iter().map(|r| r.into_pyarray(py)).collect())
}

/// k-infinity by power iteration (port of `pincell.infinite.solve`):
/// returns (k, flux, iterations).
#[pyfunction]
#[pyo3(signature = (total, scatter, fission, tol=1e-13))]
fn infinite_medium<'py>(
    py: Python<'py>,
    total: Vec<f64>,
    scatter: Vec<f64>,
    fission: Vec<f64>,
    tol: f64,
) -> PyResult<(f64, Array<'py>, usize)> {
    let n = total.len();
    if scatter.len() != n * n || fission.len() != n * n {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "matrices must be n x n",
        ));
    }
    let (k, phi, it) = py
        .detach(|| infinite::k_infinity(&total, &scatter, &fission, tol))
        .map_err(pyo3::exceptions::PyValueError::new_err)?;
    Ok((k, phi.into_pyarray(py), it))
}

/// Segments of many tracks through a lattice (port of
/// `transport.trace.trace_track`): (regions, lengths, start per track).
#[pyfunction]
#[allow(clippy::too_many_arguments)]
#[pyo3(signature = (pitch, nx, ny, radii, sectors, offsets, x0, y0, x1, y1, gap=0.0))]
fn trace_tracks<'py>(
    py: Python<'py>,
    pitch: f64,
    nx: usize,
    ny: usize,
    radii: Vec<Vec<f64>>,
    sectors: Vec<usize>,
    offsets: Vec<usize>,
    x0: PyReadonlyArray1<'py, f64>,
    y0: PyReadonlyArray1<'py, f64>,
    x1: PyReadonlyArray1<'py, f64>,
    y1: PyReadonlyArray1<'py, f64>,
    gap: f64,
) -> PyResult<(IntArray<'py>, Array<'py>, IntArray<'py>)> {
    let cells = nx * ny;
    if radii.len() != cells || sectors.len() != cells || offsets.len() != cells {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "radii, sectors and offsets need one entry per cell",
        ));
    }
    if sectors.contains(&0) {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "sectors must be >= 1",
        ));
    }
    let (x0, y0, x1, y1) = (
        x0.as_slice()?,
        y0.as_slice()?,
        x1.as_slice()?,
        y1.as_slice()?,
    );
    if y0.len() != x0.len() || x1.len() != x0.len() || y1.len() != x0.len() {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "track arrays differ in length",
        ));
    }
    let lattice = raytrace::Lattice {
        pitch,
        gap,
        nx,
        ny,
        radii,
        sectors,
        offsets,
    };
    let (regions, lengths, start) = py.detach(|| raytrace::trace_all(&lattice, x0, y0, x1, y1));
    Ok((
        regions.into_pyarray(py),
        lengths.into_pyarray(py),
        start.into_pyarray(py),
    ))
}

/// One MOC sweep over all tracks (port of `transport.moc.sweep`):
/// returns (tally per region x group, next incoming fluxes, leakage, net
/// currents per face x group; empty without `faces`).
#[pyfunction]
#[pyo3(signature = (region, length, start, next_fwd, fwd_forward, next_bwd, bwd_forward, weights, sin, sigma_t, source, psi_in, groups, faces=None))]
#[allow(clippy::too_many_arguments)]
fn moc_sweep<'py>(
    py: Python<'py>,
    region: PyReadonlyArray1<'py, i64>,
    length: PyReadonlyArray1<'py, f64>,
    start: PyReadonlyArray1<'py, i64>,
    next_fwd: PyReadonlyArray1<'py, i64>,
    fwd_forward: PyReadonlyArray1<'py, i64>,
    next_bwd: PyReadonlyArray1<'py, i64>,
    bwd_forward: PyReadonlyArray1<'py, i64>,
    weights: PyReadonlyArray1<'py, f64>,
    sin: PyReadonlyArray1<'py, f64>,
    sigma_t: PyReadonlyArray1<'py, f64>,
    source: PyReadonlyArray1<'py, f64>,
    psi_in: PyReadonlyArray1<'py, f64>,
    groups: usize,
    faces: Option<FacesIn<'py>>,
) -> PyResult<(Array<'py>, Array<'py>, f64, Array<'py>)> {
    let tracks = moc::Tracks {
        region: region.as_slice()?,
        length: length.as_slice()?,
        start: start.as_slice()?,
        next_fwd: next_fwd.as_slice()?,
        fwd_forward: fwd_forward.as_slice()?,
        next_bwd: next_bwd.as_slice()?,
        bwd_forward: bwd_forward.as_slice()?,
        weights: weights.as_slice()?,
        sin: sin.as_slice()?,
        groups,
    };
    let n = tracks.start.len().saturating_sub(1);
    let (sigma_t, source, psi_in) = (sigma_t.as_slice()?, source.as_slice()?, psi_in.as_slice()?);
    let bad = |m: &str| Err(pyo3::exceptions::PyValueError::new_err(m.to_string()));
    if groups == 0 || sigma_t.len() % groups != 0 || source.len() != sigma_t.len() {
        return bad("sigma_t and source must be regions x groups");
    }
    if tracks.region.len() != tracks.length.len()
        || tracks.next_fwd.len() != n
        || tracks.fwd_forward.len() != n
        || tracks.next_bwd.len() != n
        || tracks.bwd_forward.len() != n
        || tracks.weights.len() != n * tracks.sin.len()
        || psi_in.len() != n * 2 * tracks.sin.len() * groups
    {
        return bad("track arrays have inconsistent lengths");
    }
    let regions = (sigma_t.len() / groups) as i64;
    if tracks.region.iter().any(|&r| r < 0 || r >= regions) {
        return bad("segment region out of range");
    }
    let faces = match &faces {
        None => None,
        Some((start, position, face, sign, n_faces)) => {
            let f = moc::Faces {
                start: start.as_slice()?,
                position: position.as_slice()?,
                face: face.as_slice()?,
                sign: sign.as_slice()?,
                n_faces: *n_faces,
            };
            if f.start.len() != n + 1
                || f.position.len() != f.face.len()
                || f.sign.len() != f.face.len()
                || f.face.iter().any(|&x| x < 0 || x as usize >= f.n_faces)
            {
                return bad("inconsistent face crossings");
            }
            Some(f)
        }
    };
    let (tally, next, leakage, currents) =
        py.detach(|| moc::sweep(&tracks, faces.as_ref(), sigma_t, source, psi_in));
    Ok((
        tally.into_pyarray(py),
        next.into_pyarray(py),
        leakage,
        currents.into_pyarray(py),
    ))
}

/// Coarse-mesh eigenvalue solve (port of `transport.cmfd.solve_coarse`):
/// arrays flattened (groups, ny, nx); inscatter (g', g, ny, nx).
/// Returns (k, phi, outer iterations).
#[pyfunction]
#[allow(clippy::too_many_arguments)]
#[pyo3(signature = (groups, ny, nx, diag, east, west, north, south, inscatter, fission, chi, phi, k, tol=1e-10, max_outer=5000, inner=4, omega=1.5))]
fn cmfd_solve<'py>(
    py: Python<'py>,
    groups: usize,
    ny: usize,
    nx: usize,
    diag: PyReadonlyArray1<'py, f64>,
    east: PyReadonlyArray1<'py, f64>,
    west: PyReadonlyArray1<'py, f64>,
    north: PyReadonlyArray1<'py, f64>,
    south: PyReadonlyArray1<'py, f64>,
    inscatter: PyReadonlyArray1<'py, f64>,
    fission: PyReadonlyArray1<'py, f64>,
    chi: PyReadonlyArray1<'py, f64>,
    phi: PyReadonlyArray1<'py, f64>,
    k: f64,
    tol: f64,
    max_outer: usize,
    inner: usize,
    omega: f64,
) -> PyResult<(f64, Array<'py>, usize)> {
    let n = groups * ny * nx;
    let op = cmfd::Operator {
        groups,
        ny,
        nx,
        diag: diag.as_slice()?,
        east: east.as_slice()?,
        west: west.as_slice()?,
        north: north.as_slice()?,
        south: south.as_slice()?,
        inscatter: inscatter.as_slice()?,
        fission: fission.as_slice()?,
        chi: chi.as_slice()?,
    };
    let phi = phi.as_slice()?;
    let sizes = [
        op.diag.len(),
        op.east.len(),
        op.west.len(),
        op.north.len(),
        op.south.len(),
        op.fission.len(),
        op.chi.len(),
        phi.len(),
    ];
    if n == 0 || sizes.iter().any(|&l| l != n) || op.inscatter.len() != n * groups {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "arrays must be groups x ny x nx (inscatter groups x groups x ny x nx)",
        ));
    }
    let settings = cmfd::Settings {
        tol,
        max_outer,
        inner,
        omega,
    };
    let (k, phi, outer) = py.detach(|| cmfd::solve(&op, phi, k, &settings));
    Ok((k, phi.into_pyarray(py), outer))
}

/// Red-black SOR sweeps of the nodal diffusion operator for one group (port
/// of `fission.nodal.diffusion._sweep`, `count` times both colours). Arrays
/// are flattened [nz, ny, nx]; coef = x lower, x upper, y lower, y upper,
/// z lower, z upper. Returns the new flux.
#[pyfunction]
#[allow(clippy::too_many_arguments)]
fn diffusion_sweeps<'py>(
    py: Python<'py>,
    nz: usize,
    ny: usize,
    nx: usize,
    diag: PyReadonlyArray1<'py, f64>,
    coef: [PyReadonlyArray1<'py, f64>; 6],
    inside: PyReadonlyArray1<'py, u8>,
    b: PyReadonlyArray1<'py, f64>,
    phi: PyReadonlyArray1<'py, f64>,
    count: usize,
    omega: f64,
) -> PyResult<Array<'py>> {
    let n = nz * ny * nx;
    let c: Vec<&[f64]> = coef.iter().map(|a| a.as_slice()).collect::<Result<_, _>>()?;
    let (diag, inside, b) = (diag.as_slice()?, inside.as_slice()?, b.as_slice()?);
    let mut phi = phi.as_slice()?.to_vec();
    let lengths = [diag.len(), inside.len(), b.len(), phi.len()];
    if lengths.iter().copied().chain(c.iter().map(|a| a.len())).any(|l| l != n) {
        return Err(pyo3::exceptions::PyValueError::new_err("arrays must have nz * ny * nx entries"));
    }
    let grid = diffusion::Grid { nz, ny, nx };
    let op = diffusion::Operator {
        diag,
        coef: [(c[0], c[1]), (c[2], c[3]), (c[4], c[5])],
        inside,
    };
    py.detach(|| diffusion::sweeps(&grid, &op, b, &mut phi, count, omega));
    Ok(phi.into_pyarray(py))
}

/// Ultra-fine slowing-down march of R regions (`shielding.march`).
/// Returns (psi flattened (n, R), psi_M (n)).

type KernelOut<'py> = (f64, f64, f64, usize, usize, Bound<'py, PyArray1<i64>>,
                       Bound<'py, PyArray1<i64>>, Array<'py>, Array<'py>);

/// Free-gas elastic kernel with the 0 K cross section (port of
/// `fission.pincell.moving.free_gas_kernel`): returns the MovingKernel fields
/// after (nuclide, awr, temperature).
#[pyfunction]
#[allow(clippy::too_many_arguments)]
fn free_gas_kernel<'py>(
    py: Python<'py>,
    energies: PyReadonlyArray1<'py, f64>,
    sigma0: PyReadonlyArray1<'py, f64>,
    awr: f64,
    temperature: f64,
    e_low: f64,
    e_top: f64,
    dx: f64,
    dz: f64,
) -> PyResult<KernelOut<'py>> {
    let (e, s) = (energies.as_slice()?, sigma0.as_slice()?);
    let k = py
        .detach(|| moving::free_gas_kernel(e, s, awr, temperature, e_low, e_top, dx, dz))
        .map_err(pyo3::exceptions::PyValueError::new_err)?;
    Ok((k.scale, k.base, k.dx, k.top, k.nodes, k.starts.into_pyarray(py),
        k.rows.into_pyarray(py), k.values.into_pyarray(py), k.norm.into_pyarray(py)))
}

/// A `MovingKernel` read by attribute.
#[derive(FromPyObject)]
struct PyKernel<'py> {
    scale: f64,
    base: f64,
    dx: f64,
    top: usize,
    nodes: usize,
    starts: PyReadonlyArray1<'py, i64>,
    rows: PyReadonlyArray1<'py, i64>,
    values: PyReadonlyArray1<'py, f64>,
}

#[pyfunction]
#[allow(clippy::too_many_arguments)]
#[pyo3(signature = (energies, total, scatter, alphas, escape, collide, volumes, history,
                    start, share, outside, kernels=Vec::new(), tol=1e-10, max_passes=200))]
fn slowing_down_march<'py>(
    py: Python<'py>,
    energies: PyReadonlyArray1<'py, f64>,
    total: PyReadonlyArray1<'py, f64>,
    scatter: PyReadonlyArray1<'py, f64>,
    alphas: PyReadonlyArray1<'py, f64>,
    escape: PyReadonlyArray1<'py, f64>,
    collide: PyReadonlyArray1<'py, f64>,
    volumes: PyReadonlyArray1<'py, f64>,
    history: PyReadonlyArray1<'py, f64>,
    start: usize,
    share: f64,
    outside: PyReadonlyArray1<'py, f64>,
    kernels: Vec<Option<PyKernel<'py>>>,
    tol: f64,
    max_passes: usize,
) -> PyResult<(Array<'py>, Array<'py>)> {
    let mut refs = Vec::with_capacity(kernels.len());
    for k in &kernels {
        refs.push(match k {
            None => None,
            Some(k) => Some(slowing::Kernel {
                scale: k.scale,
                base: k.base,
                dx: k.dx,
                top: k.top,
                nodes: k.nodes,
                starts: k.starts.as_slice()?,
                rows: k.rows.as_slice()?,
                values: k.values.as_slice()?,
            }),
        });
    }
    let inp = slowing::Input {
        energies: energies.as_slice()?,
        total: total.as_slice()?,
        scatter: scatter.as_slice()?,
        alphas: alphas.as_slice()?,
        escape: escape.as_slice()?,
        collide: collide.as_slice()?,
        volumes: volumes.as_slice()?,
        history: history.as_slice()?,
        start,
        share,
        outside: outside.as_slice()?,
        kernels: refs,
        tol,
        max_passes,
    };
    let (psi, levels) = py
        .detach(|| slowing::march(&inp))
        .map_err(pyo3::exceptions::PyValueError::new_err)?;
    Ok((psi.into_pyarray(py), levels.into_pyarray(py)))
}

/// The n-point Gauss-Legendre rule (for parity tests with numpy's leggauss).
#[pyfunction]
fn gauss_legendre<'py>(py: Python<'py>, n: usize) -> (Array<'py>, Array<'py>) {
    let (x, w) = gauss::gauss_legendre(n);
    (x.into_pyarray(py), w.into_pyarray(py))
}

/// Correctly rounded sum (port of `math.fsum`, for parity tests).
#[pyfunction(name = "fsum")]
fn py_fsum(values: Vec<f64>) -> f64 {
    fsum::fsum(values)
}

// ------------------------------------------------------------ Doppler
fn columns_to_arrays<'py>(py: Python<'py>, n: usize, rows: &[Vec<f64>]) -> Vec<Array<'py>> {
    (0..n)
        .map(|m| {
            rows.iter()
                .map(|r| r[m])
                .collect::<Vec<f64>>()
                .into_pyarray(py)
        })
        .collect()
}

fn check_columns(energies: &[f64], columns: &[Vec<f64>]) -> PyResult<()> {
    if energies.len() < 2 || columns.iter().any(|c| c.len() != energies.len()) {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "need at least two energies and one value per energy in every column",
        ));
    }
    Ok(())
}

/// Broadening weights at reduced speed y: (index of the first weight, weights).
#[pyfunction]
fn doppler_weights<'py>(py: Python<'py>, x: Vec<f64>, y: f64) -> (usize, Array<'py>) {
    let (start, w) = doppler::broadening_weights(&x, y);
    (start, w.into_pyarray(py))
}

/// Broadened columns at the energies `at` (port of `broadened_at`).
#[pyfunction]
fn doppler_evaluate<'py>(
    py: Python<'py>,
    energies: Vec<f64>,
    columns: Vec<Vec<f64>>,
    awr: f64,
    temperature: f64,
    at: Vec<f64>,
) -> PyResult<Vec<Array<'py>>> {
    check_columns(&energies, &columns)?;
    let rows = py.detach(|| {
        let cold = doppler::ColdTable::new(&energies, &columns, awr, temperature);
        at.iter().map(|&e| cold.eval(e)).collect::<Vec<_>>()
    });
    Ok(columns_to_arrays(py, columns.len(), &rows))
}

/// Indices kept by thinning; `columns` hold one array per reaction (port of `thin`).
#[pyfunction]
fn doppler_thin(
    energies: Vec<f64>,
    columns: Vec<Vec<f64>>,
    tol: f64,
    atol: f64,
) -> PyResult<Vec<usize>> {
    check_columns(&energies, &columns)?;
    let rows: Vec<Vec<f64>> = (0..energies.len())
        .map(|k| columns.iter().map(|c| c[k]).collect())
        .collect();
    Ok(doppler::thin(&energies, &rows, tol, atol))
}

/// Broaden a table: evaluate at `seeds`, refine, thin.
#[pyfunction]
#[pyo3(signature = (energies, columns, awr, temperature, seeds, tol=1e-3, atol=1e-5, min_rel_width=1e-9))]
#[allow(clippy::too_many_arguments)]
fn doppler_broaden_table<'py>(
    py: Python<'py>,
    energies: Vec<f64>,
    columns: Vec<Vec<f64>>,
    awr: f64,
    temperature: f64,
    seeds: Vec<f64>,
    tol: f64,
    atol: f64,
    min_rel_width: f64,
) -> PyResult<(Array<'py>, Vec<Array<'py>>)> {
    check_columns(&energies, &columns)?;
    if seeds.len() < 2 {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "need at least two seeds",
        ));
    }
    let (grid, rows) = py.detach(|| {
        let cold = doppler::ColdTable::new(&energies, &columns, awr, temperature);
        doppler::broaden_table(&cold, &seeds, tol, atol, min_rel_width)
    });
    Ok((
        grid.into_pyarray(py),
        columns_to_arrays(py, columns.len(), &rows),
    ))
}

/// Evaluate a Python `Tabulated1D` at many points (for parity tests).
#[pyfunction]
fn tab1_eval<'py>(
    py: Python<'py>,
    table: tab1::Tab1,
    xs: PyReadonlyArray1<'py, f64>,
) -> PyResult<Array<'py>> {
    let out: Vec<f64> = xs.as_slice()?.iter().map(|&x| table.eval(x)).collect();
    Ok(out.into_pyarray(py))
}

/// The wave-number constant (for parity tests).
#[pyfunction]
fn k_constant() -> f64 {
    channel::k_constant()
}

/// Multilinear interpolation of a table at many points (port of
/// `fission.tables.multilinear`): nodes per axis, values flattened
/// [n_0, ..., n_d-1, m], points flattened [p, d]. Returns [p, m] flattened.
#[pyfunction]
fn table_multilinear<'py>(
    py: Python<'py>,
    nodes: Vec<PyReadonlyArray1<'py, f64>>,
    values: PyReadonlyArray1<'py, f64>,
    m: usize,
    points: PyReadonlyArray1<'py, f64>,
) -> PyResult<Array<'py>> {
    let axes: Vec<&[f64]> = nodes.iter().map(|a| a.as_slice()).collect::<Result<_, _>>()?;
    let (values, points) = (values.as_slice()?, points.as_slice()?);
    let size: usize = axes.iter().map(|a| a.len()).product::<usize>() * m;
    if axes.is_empty() || axes.iter().any(|a| a.is_empty()) || values.len() != size
        || points.len() % axes.len() != 0
    {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "values must be the grid times m, points a multiple of the axis count",
        ));
    }
    let out = py.detach(|| tables::multilinear(&axes, values, m, points));
    Ok(out.into_pyarray(py))
}

/// Water temperatures (K) of enthalpies h (J/kg) at pressure p (Pa), IF97
/// region 1 by Newton from `guess` (`feedback.water.temperature`); NaN where
/// Newton does not converge.
#[pyfunction]
fn water_temperatures<'py>(
    py: Python<'py>,
    p: f64,
    h: PyReadonlyArray1<'py, f64>,
    guess: f64,
) -> PyResult<Array<'py>> {
    let h = h.as_slice()?;
    let out: Vec<f64> = py.detach(|| {
        h.iter().map(|&x| water::temperature(p, x, guess).unwrap_or(f64::NAN)).collect()
    });
    Ok(out.into_pyarray(py))
}

/// Water densities (kg/m^3) at p (Pa) and temperatures t (K).
#[pyfunction]
fn water_densities<'py>(
    py: Python<'py>,
    p: f64,
    t: PyReadonlyArray1<'py, f64>,
) -> PyResult<Array<'py>> {
    let out: Vec<f64> = t.as_slice()?.iter().map(|&x| water::density(p, x)).collect();
    Ok(out.into_pyarray(py))
}

/// Dittus-Boelter clad-to-water coefficients (W/(m^2 K)) at p (Pa) and
/// water temperatures t (K) (`feedback.convection.film_coefficient`).
#[pyfunction]
#[pyo3(signature = (p, t, mass_flux, diameter, pitch_over_diameter=None))]
fn water_film_coefficients<'py>(
    py: Python<'py>,
    p: f64,
    t: PyReadonlyArray1<'py, f64>,
    mass_flux: f64,
    diameter: f64,
    pitch_over_diameter: Option<f64>,
) -> PyResult<Array<'py>> {
    let t = t.as_slice()?;
    let out: Vec<f64> = py.detach(|| {
        t.iter()
            .map(|&x| water::film_coefficient(p, x, mass_flux, diameter, pitch_over_diameter))
            .collect()
    });
    Ok(out.into_pyarray(py))
}

/// Fuel rod temperatures for linear heat rates q (W/m), coolant
/// temperatures (K) and optional film coefficients (`feedback.fuel.rod`);
/// radii in metres. Returns (centre, pellet surface, clad inner, clad outer).
#[pyfunction]
#[allow(clippy::too_many_arguments)]
#[pyo3(signature = (q, coolant, pellet_radius, clad_inner, clad_outer, gap_conductance, film=None))]
fn fuel_rods<'py>(
    py: Python<'py>,
    q: PyReadonlyArray1<'py, f64>,
    coolant: PyReadonlyArray1<'py, f64>,
    pellet_radius: f64,
    clad_inner: f64,
    clad_outer: f64,
    gap_conductance: f64,
    film: Option<PyReadonlyArray1<'py, f64>>,
) -> PyResult<(Array<'py>, Array<'py>, Array<'py>, Array<'py>)> {
    let (q, coolant) = (q.as_slice()?, coolant.as_slice()?);
    let film = film.as_ref().map(|f| f.as_slice()).transpose()?;
    if coolant.len() != q.len() || film.is_some_and(|f| f.len() != q.len()) {
        return Err(pyo3::exceptions::PyValueError::new_err("q, coolant and film differ in length"));
    }
    let rod = fuel::Rod { pellet_radius, clad_inner, clad_outer, gap_conductance };
    let mut out: [Vec<f64>; 4] = Default::default();
    for (n, (&qn, &c)) in q.iter().zip(coolant).enumerate() {
        let temps = fuel::rod(&rod, qn, c, film.map(|f| f[n]));
        for (o, t) in out.iter_mut().zip(temps) {
            o.push(t);
        }
    }
    let [a, b, c, d] = out;
    Ok((a.into_pyarray(py), b.into_pyarray(py), c.into_pyarray(py), d.into_pyarray(py)))
}

#[pymodule]
fn _rust(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    m.add_class::<PyReichMoore>()?;
    m.add_class::<PyMlbw>()?;
    m.add_function(wrap_pyfunction!(py_reconstruct, m)?)?;
    m.add_function(wrap_pyfunction!(tab1_eval, m)?)?;
    m.add_function(wrap_pyfunction!(k_constant, m)?)?;
    m.add_function(wrap_pyfunction!(pointwise_table, m)?)?;
    m.add_function(wrap_pyfunction!(py_fsum, m)?)?;
    m.add_function(wrap_pyfunction!(two_body_kernel, m)?)?;
    m.add_function(wrap_pyfunction!(two_body_table, m)?)?;
    m.add_function(wrap_pyfunction!(thermal_table, m)?)?;
    m.add_function(wrap_pyfunction!(infinite_medium, m)?)?;
    m.add_function(wrap_pyfunction!(gauss_legendre, m)?)?;
    m.add_function(wrap_pyfunction!(trace_tracks, m)?)?;
    m.add_function(wrap_pyfunction!(moc_sweep, m)?)?;
    m.add_function(wrap_pyfunction!(cmfd_solve, m)?)?;
    m.add_function(wrap_pyfunction!(diffusion_sweeps, m)?)?;
    m.add_function(wrap_pyfunction!(slowing_down_march, m)?)?;
    m.add_function(wrap_pyfunction!(free_gas_kernel, m)?)?;
    m.add_function(wrap_pyfunction!(group_hat_integrals, m)?)?;
    m.add_function(wrap_pyfunction!(doppler_weights, m)?)?;
    m.add_function(wrap_pyfunction!(doppler_evaluate, m)?)?;
    m.add_function(wrap_pyfunction!(doppler_thin, m)?)?;
    m.add_function(wrap_pyfunction!(doppler_broaden_table, m)?)?;
    m.add_function(wrap_pyfunction!(table_multilinear, m)?)?;
    m.add_function(wrap_pyfunction!(water_temperatures, m)?)?;
    m.add_function(wrap_pyfunction!(water_densities, m)?)?;
    m.add_function(wrap_pyfunction!(water_film_coefficients, m)?)?;
    m.add_function(wrap_pyfunction!(fuel_rods, m)?)?;
    Ok(())
}
