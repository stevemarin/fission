//! Full-range pointwise tables: the function refined by `pointwise()` in
//! `fission/endf/pointwise.py` (MF3 backgrounds + resonance part, total as
//! the exact sum of the partials).

use crate::Resonances;
use crate::fsum::fsum;
use crate::tab1::Tab1;

/// The resonance part: a formalism, its energy range and the value slots
/// of (elastic, capture, fission).
pub struct ResonancePart<'a> {
    pub model: &'a dyn Resonances,
    pub el: f64,
    pub eh: f64,
    pub slots: [Option<usize>; 3],
}

/// Everything needed to evaluate all reactions at one energy.
pub struct Reactions<'a> {
    pub tables: Vec<Option<Tab1>>, // MF3 table per slot (None: total or absent)
    pub resonance: Option<ResonancePart<'a>>,
    pub total: usize,
    pub partials: Vec<usize>,
}

impl Reactions<'_> {
    /// Mirrors `f` inside `pointwise()`.
    pub fn eval(&self, e: f64) -> Vec<f64> {
        let mut v: Vec<f64> = self
            .tables
            .iter()
            .map(|t| t.as_ref().map_or(0.0, |t| t.eval(e)))
            .collect();
        if let Some(r) = &self.resonance
            && r.el <= e
            && e <= r.eh
        {
            let x = r.model.eval(e);
            for (slot, value) in r.slots.iter().zip(x) {
                if let Some(s) = slot {
                    v[*s] += value;
                }
            }
        }
        v[self.total] = fsum(self.partials.iter().map(|&s| v[s]));
        v
    }
}
