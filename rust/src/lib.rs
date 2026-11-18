//! Compute kernels for `fission`, exposed to Python as `fission._rust`.
//!
//! Every kernel here is a port of a Python reference implementation in
//! `src/fission`; tests/test_rust_parity.py checks that the two agree.

pub mod channel;
pub mod complex;
pub mod mlbw;
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
    let (energies, values) =
        py.detach(|| reconstruct::adaptive_grid(f, &seeds, tol, atol, min_rel_width));
    let pick = |j: usize| values.iter().map(|v| v[j]).collect::<Vec<f64>>();
    Ok((
        energies.into_pyarray(py),
        pick(0).into_pyarray(py),
        pick(1).into_pyarray(py),
        pick(2).into_pyarray(py),
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
    Ok(())
}
