//! Compute kernels for `fission`, exposed to Python as `fission._rust`.
//!
//! Every kernel here is a port of a Python reference implementation in
//! `src/fission`; tests/test_rust_parity.py checks that the two agree.

pub mod channel;
pub mod complex;
pub mod doppler;
pub mod fsum;
pub mod mlbw;
pub mod pointwise;
pub mod reconstruct;
pub mod reich_moore;
pub mod spin_groups;
pub mod tab1;

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
    m.add_function(wrap_pyfunction!(doppler_weights, m)?)?;
    m.add_function(wrap_pyfunction!(doppler_evaluate, m)?)?;
    m.add_function(wrap_pyfunction!(doppler_thin, m)?)?;
    m.add_function(wrap_pyfunction!(doppler_broaden_table, m)?)?;
    Ok(())
}
