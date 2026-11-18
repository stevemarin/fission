//! Doppler broadening with the free-gas kernel, plus thinning.
//! Port of `fission/processing/doppler.py`; see the theory there.

use crate::reconstruct::{Unsplittable, adaptive_grid_parallel};

pub const BOLTZMANN_EV_PER_K: f64 = 8.617333262e-5;
pub const Z_CUTOFF: f64 = 6.0;
pub const PANEL_WIDTH: f64 = 0.5;
// Gauss-Legendre rules on [-1, 1], identical to GAUSS_RULES in doppler.py
const NODES_2: [f64; 2] = [-0.5773502691896257, 0.5773502691896257];
const WEIGHTS_2: [f64; 2] = [1.0, 1.0];
const NODES_4: [f64; 4] = [
    -0.8611363115940526,
    -0.33998104358485626,
    0.33998104358485626,
    0.8611363115940526,
];
const WEIGHTS_4: [f64; 4] = [
    0.34785484513745357,
    0.6521451548625464,
    0.6521451548625464,
    0.34785484513745357,
];
pub const NODES_8: [f64; 8] = [
    -0.9602898564975362,
    -0.7966664774136267,
    -0.525532409916329,
    -0.18343464249564978,
    0.18343464249564978,
    0.525532409916329,
    0.7966664774136267,
    0.9602898564975362,
];
pub const WEIGHTS_8: [f64; 8] = [
    0.10122853629037706,
    0.22238103445337443,
    0.3137066458778869,
    0.36268378337836166,
    0.36268378337836166,
    0.3137066458778869,
    0.22238103445337443,
    0.10122853629037706,
];

/// Fewest Gauss points for a panel of width h (mirrors `gauss_rule`).
fn gauss_rule(h: f64) -> (&'static [f64], &'static [f64]) {
    if h <= 0.01 {
        (&NODES_2, &WEIGHTS_2)
    } else if h <= 0.1 {
        (&NODES_4, &WEIGHTS_4)
    } else {
        (&NODES_8, &WEIGHTS_8)
    }
}

/// alpha = AWR / (k_B T) in 1/eV.
pub fn alpha(awr: f64, temperature: f64) -> f64 {
    awr / (BOLTZMANN_EV_PER_K * temperature)
}

/// K(x, y) = exp(-(x-y)^2) - exp(-(x+y)^2).
#[inline]
pub fn kernel(x: f64, y: f64) -> f64 {
    let d = x - y;
    let z = 4.0 * x * y;
    // for z > 40, exp(-z) < 5e-18 and 1 - exp(-z) rounds to exactly 1.0
    let tail = if z > 40.0 { 1.0 } else { -(-z).exp_m1() };
    (-(d * d)).exp() * tail
}

/// Panels of at most PANEL_WIDTH covering [a, b]: (count, width).
fn panels(a: f64, b: f64) -> (usize, f64) {
    let n = (((b - a) / PANEL_WIDTH).ceil() as usize).max(1);
    (n, (b - a) / n as f64)
}

/// Int_a^b f(x) K(x, y) dx (mirrors `_integrate`).
fn integrate(f: impl Fn(f64) -> f64, a: f64, b: f64, y: f64) -> f64 {
    let (n, h) = panels(a, b);
    let mut total = 0.0;
    for p in 0..n {
        let mid = a + (p as f64 + 0.5) * h;
        let (nodes, weights) = gauss_rule(h);
        for (t, w) in nodes.iter().zip(weights) {
            let x = mid + 0.5 * h * t;
            total += w * f(x) * kernel(x, y);
        }
    }
    0.5 * h * total
}

/// Contributions of the end values of segment [xl, xr], over [a, b] (mirrors `_segment_weights`).
fn segment_weights(xl: f64, xr: f64, a: f64, b: f64, y: f64) -> (f64, f64) {
    let d = (xr - xl) * (xr + xl);
    let (n, h) = panels(a, b);
    let (mut left, mut right) = (0.0, 0.0);
    for p in 0..n {
        let mid = a + (p as f64 + 0.5) * h;
        let (nodes, weights) = gauss_rule(h);
        for (t, w) in nodes.iter().zip(weights) {
            let x = mid + 0.5 * h * t;
            let f = w * x * x * kernel(x, y) / d;
            left += f * (xr - x) * (xr + x);
            right += f * (x - xl) * (x + xl);
        }
    }
    (0.5 * h * left, 0.5 * h * right)
}

/// Weights for sigma(y) = sum_k w[k] sigma_0[start + k] (mirrors `broadening_weights`).
/// The nodes that contribute are always contiguous.
pub fn broadening_weights(x: &[f64], y: f64) -> (usize, Vec<f64>) {
    let (lo, hi) = ((y - Z_CUTOFF).max(0.0), y + Z_CUTOFF);
    let n = x.len();
    let start = x.partition_point(|&v| v <= lo).saturating_sub(1);
    let mut w: Vec<f64> = Vec::new();
    let at = |i: usize, v: f64, w: &mut Vec<f64>| {
        if w.len() <= i - start {
            w.resize(i - start + 1, 0.0);
        }
        w[i - start] += v;
    };

    if lo < x[0] {
        let x0 = x[0];
        at(0, integrate(|t| x0 * t, lo, x0.min(hi), y), &mut w);
    }
    let mut i = start;
    while i + 1 < n && x[i] < hi {
        let (a, b) = (x[i].max(lo), x[i + 1].min(hi));
        if b > a {
            let (left, right) = segment_weights(x[i], x[i + 1], a, b, y);
            at(i, left, &mut w);
            at(i + 1, right, &mut w);
        }
        i += 1;
    }
    if hi > x[n - 1] {
        at(n - 1, integrate(|t| t * t, x[n - 1].max(lo), hi, y), &mut w);
    }

    let norm = 1.0 / (y * y * std::f64::consts::PI.sqrt());
    for v in w.iter_mut() {
        *v *= norm;
    }
    (start, w)
}

/// A 0 K table (reduced speeds and value columns) ready for broadening.
pub struct ColdTable<'a> {
    alpha: f64,
    x: Vec<f64>,
    columns: &'a [Vec<f64>],
}

impl<'a> ColdTable<'a> {
    pub fn new(energies: &[f64], columns: &'a [Vec<f64>], awr: f64, temperature: f64) -> Self {
        let alpha = alpha(awr, temperature);
        let x = energies.iter().map(|&e| (alpha * e).sqrt()).collect();
        ColdTable { alpha, x, columns }
    }

    /// Broadened value of every column at one energy.
    pub fn eval(&self, energy: f64) -> Vec<f64> {
        let (start, w) = broadening_weights(&self.x, (self.alpha * energy).sqrt());
        self.columns
            .iter()
            .map(|col| {
                w.iter()
                    .zip(&col[start..start + w.len()])
                    .map(|(wi, s)| wi * s)
                    .sum()
            })
            .collect()
    }
}

/// Do all points strictly between i and j lie within tolerance of the chord? (`_chord_holds`)
fn chord_holds(e: &[f64], v: &[Vec<f64>], i: usize, j: usize, tol: f64, atol: f64) -> bool {
    (i + 1..j).all(|k| {
        let t = (e[k] - e[i]) / (e[j] - e[i]);
        v[k].iter()
            .zip(&v[i])
            .zip(&v[j])
            .all(|((&vk, &vi), &vj)| (vk - (vi + t * (vj - vi))).abs() <= tol * vk.abs() + atol)
    })
}

/// Indices of the points kept by greedy thinning (mirrors `thin`).
pub fn thin(energies: &[f64], values: &[Vec<f64>], tol: f64, atol: f64) -> Vec<usize> {
    let n = energies.len();
    let mut keep = vec![0];
    let mut i = 0;
    while i + 1 < n {
        let mut j = i + 1;
        while j + 1 < n && chord_holds(energies, values, i, j + 1, tol, atol) {
            j += 1;
        }
        keep.push(j);
        i = j;
    }
    keep
}

/// Broaden, refine and thin (mirrors the Python branch of `broaden_table`).
/// `seeds` come from `doppler_seeds`. Returns (energies, rows).
pub fn broaden_table(
    cold: &ColdTable,
    seeds: &[f64],
    tol: f64,
    atol: f64,
    min_rel_width: f64,
) -> (Vec<f64>, Vec<Vec<f64>>) {
    let (grid, values) = adaptive_grid_parallel(
        |e| cold.eval(e),
        seeds,
        tol,
        atol,
        min_rel_width,
        &Unsplittable::new(),
    );
    let keep = thin(&grid, &values, tol, atol);
    let energies = keep.iter().map(|&k| grid[k]).collect();
    let rows = keep.into_iter().map(|k| values[k].clone()).collect();
    (energies, rows)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn one_over_v_is_unchanged() {
        let e: Vec<f64> = (0..=7000)
            .map(|i| 1e-5 * 10f64.powf(i as f64 / 1000.0))
            .collect();
        let col = vec![e.iter().map(|v| 1.0 / v.sqrt()).collect::<Vec<f64>>()];
        let cold = ColdTable::new(&e, &col, 1.0, 600.0);
        for probe in [1e-4, 0.0253, 3.0] {
            let got = cold.eval(probe)[0];
            assert!((got * probe.sqrt() - 1.0).abs() < 1e-5, "{probe}: {got}");
        }
    }

    #[test]
    fn thinning_a_straight_line_keeps_the_ends() {
        let e = [1.0, 2.0, 3.0, 4.0];
        let v: Vec<Vec<f64>> = e.iter().map(|x| vec![2.0 * x, 5.0]).collect();
        assert_eq!(thin(&e, &v, 1e-3, 0.0), vec![0, 3]);
    }
}
