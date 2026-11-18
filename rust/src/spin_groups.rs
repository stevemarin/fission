//! Inputs shared by the resonance formalisms, read directly from the Python
//! objects built by `nodal.endf.reich_moore.spin_groups` / `potential_terms`.

use crate::channel::{penetrability_shift_phase, wave_number};
use std::f64::consts::PI;

/// Read from a Python `Resonance` by attribute name.
#[derive(Clone, Debug, pyo3::FromPyObject)]
pub struct ResonanceIn {
    pub energy: f64,
    pub gn: f64,
    pub gg: f64,
    pub gfa: f64,
    pub gfb: f64,
}

/// Read from a Python `SpinGroup`.
#[derive(Clone, Debug, pyo3::FromPyObject)]
pub struct SpinGroupIn {
    pub l: usize,
    pub j: f64,
    pub awri: f64,
    pub a: f64,
    pub ap: f64,
    pub resonances: Vec<ResonanceIn>,
}

/// Read from a Python `PotentialTerm`.
#[derive(Clone, Debug, pyo3::FromPyObject)]
pub struct PotentialTerm {
    pub l: usize,
    pub weight: f64,
    pub awri: f64,
    pub a: f64,
    pub ap: f64,
}

impl PotentialTerm {
    pub fn eval(&self, energy: f64) -> f64 {
        let k = wave_number(energy, self.awri);
        let (_, _, phi) = penetrability_shift_phase(self.l, k * self.a, k * self.ap);
        4.0 * PI / (k * k) * self.weight * phi.sin().powi(2)
    }
}
