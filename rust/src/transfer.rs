//! Two-body transfer kernels (elastic and discrete-level scattering).
//! Port of `two_body_kernel` and helpers in `fission/processing/groupr.py`.

use std::sync::OnceLock;

use crate::gauss::{gauss_legendre, legendre, legendre_series};
use crate::groupr::find_group;
use crate::reconstruct::{Unsplittable, adaptive_grid_parallel};
use crate::tab1::Tab1;

/// One incident energy of an angular distribution.
#[derive(Clone, Debug)]
pub enum Panel {
    /// a_1..a_NL (a_0 = 1)
    Legendre(Vec<f64>),
    /// f(mu) linear between the points
    Tabulated(Vec<f64>, Vec<f64>),
}

/// A CM angular distribution (MF4 or MF6 LAW=2); no energies = isotropic.
#[derive(Clone, Debug)]
pub struct Angular {
    pub energies: Vec<f64>,
    pub breaks: Vec<usize>,
    pub laws: Vec<i64>,
    pub panels: Vec<Panel>,
}

impl Angular {
    pub fn isotropic() -> Self {
        Angular {
            energies: vec![],
            breaks: vec![],
            laws: vec![],
            panels: vec![],
        }
    }

    /// Accepts the interpolation codes `IncidentGrid.bracket` supports.
    pub fn new(
        energies: Vec<f64>,
        breaks: Vec<usize>,
        laws: Vec<i64>,
        panels: Vec<Panel>,
    ) -> Result<Self, String> {
        if laws.iter().any(|&l| ![1, 2, 3, 21, 22, 23].contains(&l)) {
            return Err(format!(
                "unsupported incident-energy interpolation {laws:?}"
            ));
        }
        if energies.len() != panels.len() {
            return Err("one panel per incident energy".into());
        }
        Ok(Angular {
            energies,
            breaks,
            laws,
            panels,
        })
    }

    fn is_isotropic(&self) -> bool {
        self.energies.is_empty()
    }

    /// (i, t) as `IncidentGrid.bracket` (the unit-base flag is not needed here).
    fn bracket(&self, e: f64) -> (usize, f64) {
        let es = &self.energies;
        if es.len() == 1 || e <= es[0] {
            return (0, 0.0);
        }
        if e >= es[es.len() - 1] {
            return (es.len() - 2, 1.0);
        }
        let i = es.partition_point(|&v| v <= e) - 1;
        let law = self.laws[self.breaks.partition_point(|&b| b < i + 2)];
        let base = if law > 20 { law - 20 } else { law };
        let t = match base {
            1 => 0.0,
            2 => (e - es[i]) / (es[i + 1] - es[i]),
            _ => (e / es[i]).ln() / (es[i + 1] / es[i]).ln(),
        };
        (i, t)
    }

    fn panel_values(panel: &Panel, mu: &[f64]) -> Vec<f64> {
        match panel {
            Panel::Legendre(c) => {
                let mut a = Vec::with_capacity(c.len() + 1);
                a.push(1.0);
                a.extend_from_slice(c);
                mu.iter().map(|&m| legendre_series(&a, m)).collect()
            }
            Panel::Tabulated(x, y) => mu.iter().map(|&m| interp(m, x, y)).collect(),
        }
    }

    /// f(mu, E) at many cosines (mirrors `_angular_values`).
    pub fn values(&self, e: f64, mu: &[f64]) -> Vec<f64> {
        if self.is_isotropic() {
            return vec![0.5; mu.len()];
        }
        let (i, t) = self.bracket(e);
        let f1 = Self::panel_values(&self.panels[i], mu);
        if t == 0.0 {
            return f1;
        }
        let f2 = Self::panel_values(&self.panels[(i + 1).min(self.panels.len() - 1)], mu);
        f1.iter()
            .zip(f2)
            .map(|(a, b)| (1.0 - t) * a + t * b)
            .collect()
    }

    /// Breaks of [m1, m2] at tabulated cosines and a Gauss order (mirrors `_cosine_pieces`).
    fn cosine_pieces(&self, e: f64, m1: f64, m2: f64) -> (Vec<f64>, usize) {
        let mut breaks = vec![m1, m2];
        let mut order = 16;
        if !self.is_isotropic() {
            let (i, _) = self.bracket(e);
            let j = (i + 1).min(self.panels.len() - 1);
            let distinct: &[usize] = if i == j { &[i][..] } else { &[i, j][..] };
            for &k in distinct {
                match &self.panels[k] {
                    Panel::Legendre(c) => order = order.max(c.len().div_ceil(2) + 4),
                    Panel::Tabulated(x, _) => {
                        breaks.extend(x.iter().copied().filter(|&m| m1 < m && m < m2))
                    }
                }
            }
        }
        breaks.sort_by(f64::total_cmp);
        breaks.dedup();
        (breaks, order)
    }
}

/// numpy.interp: linear, clamped to the end values outside [x_0, x_last].
fn interp(m: f64, x: &[f64], y: &[f64]) -> f64 {
    let n = x.len();
    if m <= x[0] {
        return y[0];
    }
    if m >= x[n - 1] {
        return y[n - 1];
    }
    let i = x.partition_point(|&v| v <= m) - 1;
    y[i] + (y[i + 1] - y[i]) * (m - x[i]) / (x[i + 1] - x[i])
}

/// Gauss-Legendre rules, computed once.
fn rule(n: usize) -> &'static (Vec<f64>, Vec<f64>) {
    static RULES: OnceLock<Vec<(Vec<f64>, Vec<f64>)>> = OnceLock::new();
    let rules = RULES.get_or_init(|| (0..=256).map(|k| gauss_legendre(k.max(1))).collect());
    &rules[n.min(256)]
}

/// K_l(E -> g') row-major as [l][g'] (mirrors `two_body_kernel`).
pub fn two_body_kernel(
    ang: &Angular,
    awr: f64,
    q: f64,
    boundaries: &[f64],
    lmax: usize,
    energy: f64,
) -> Vec<f64> {
    let g_count = boundaries.len() - 1;
    let mut k = vec![0.0; (lmax + 1) * g_count];
    let e_cm = awr / (awr + 1.0) * (awr / (awr + 1.0) * energy + q); // residual mass = AWR
    if e_cm <= 0.0 {
        return k;
    }
    let a1 = awr + 1.0;
    let c0 = e_cm + energy / (a1 * a1);
    let c1 = 2.0 * (energy * e_cm).sqrt() / a1;
    let (Some(top), Some(bottom)) = (
        find_group(boundaries, (c0 + c1).min(boundaries[0])),
        find_group(boundaries, (c0 - c1).max(boundaries[g_count])),
    ) else {
        return k;
    };
    for gp in top..=bottom {
        let (lo, hi) = (boundaries[gp + 1], boundaries[gp]);
        let m1 = (-1f64).max((lo - c0) / c1);
        let m2 = 1f64.min((hi - c0) / c1);
        if m2 <= m1 {
            continue;
        }
        let (breaks, order) = ang.cosine_pieces(energy, m1, m2);
        let (x, w) = rule(order);
        for pair in breaks.windows(2) {
            let (a, b) = (pair[0], pair[1]);
            let mu: Vec<f64> = x
                .iter()
                .map(|t| 0.5 * (a + b) + 0.5 * (b - a) * t)
                .collect();
            let f = ang.values(energy, &mu);
            for ((m, fv), wv) in mu.iter().zip(&f).zip(w) {
                let e_lab = c0 + c1 * m;
                let mu_lab = (e_cm / e_lab).sqrt() * m + (energy / e_lab).sqrt() / a1;
                let p = legendre(lmax, mu_lab.clamp(-1.0, 1.0));
                let weight = fv * 0.5 * (b - a) * wv;
                for (l, pl) in p.iter().enumerate() {
                    k[l * g_count + gp] += pl * weight;
                }
            }
        }
    }
    k
}

/// Evaluate f at every point on all cores, keeping the order.
pub fn parallel_map<T: Send, F: Fn(f64) -> T + Sync>(points: &[f64], f: F) -> Vec<T> {
    let threads = std::thread::available_parallelism().map_or(1, |n| n.get());
    let chunk = points.len().div_ceil(threads).max(1);
    std::thread::scope(|s| {
        let handles: Vec<_> = points
            .chunks(chunk)
            .map(|c| s.spawn(|| c.iter().map(|&e| f(e)).collect::<Vec<T>>()))
            .collect();
        handles
            .into_iter()
            .flat_map(|h| h.join().expect("worker panicked"))
            .collect()
    })
}

/// A tabulated two-body kernel over incident energy (mirrors `build_kernel`
/// with `two_body_kernel`, and `_fold` with the level cross section):
/// adaptive refinement from the seeds, or just the seeds.
#[allow(clippy::too_many_arguments)]
pub fn two_body_table(
    ang: &Angular,
    awr: f64,
    q: f64,
    boundaries: &[f64],
    lmax: usize,
    seeds: &[f64],
    sigma: Option<&Tab1>,
    adaptive: bool,
    tol: f64,
    atol: f64,
) -> (Vec<f64>, Vec<Vec<f64>>) {
    let f = |e: f64| {
        let mut k = two_body_kernel(ang, awr, q, boundaries, lmax, e);
        if let Some(t) = sigma {
            let s = t.eval(e);
            for v in k.iter_mut() {
                *v *= s;
            }
        }
        k
    };
    if adaptive {
        adaptive_grid_parallel(f, seeds, tol, atol, 1e-6, &Unsplittable::new())
    } else {
        (seeds.to_vec(), parallel_map(seeds, f))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn isotropic_scattering_conserves_and_has_two_thirds_a() {
        let b: Vec<f64> = (0..=40)
            .map(|i| 1e7 * (1e-7f64).powf(i as f64 / 40.0))
            .collect();
        let awr = 11.9;
        let k = two_body_kernel(&Angular::isotropic(), awr, 0.0, &b, 1, 2.5e6);
        let g = b.len() - 1;
        let p0: f64 = k[..g].iter().sum();
        let p1: f64 = k[g..].iter().sum();
        assert!((p0 - 1.0).abs() < 1e-12, "{p0}");
        assert!((p1 - 2.0 / (3.0 * awr)).abs() < 1e-8, "{p1}");
    }
}
