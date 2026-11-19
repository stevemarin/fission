//! Thermal scattering kernels (free gas and tabulated S(alpha, beta)).
//! Port of `thermal_kernel` and the models in `fission/processing/thermal.py`.

use crate::doppler::{NODES_4, NODES_8, WEIGHTS_4, WEIGHTS_8};
use crate::gauss::legendre;
use crate::groupr::find_group;

const SQRT_4PI: f64 = 3.5449077018110318; // sqrt(4 pi)

/// A thermal scattering model (mirrors `FreeGas` and `Tabulated`).
pub enum Model {
    FreeGas {
        awr: f64,
        kt: f64,
        sigma_bound: f64,
    },
    Tabulated {
        kt: f64,
        a0: f64,
        sigma_bound: f64,
        scale: f64,   // kT / kT0 for LAT = 1, else 1
        t_ratio: f64, // T_eff / T
        alpha: Vec<f64>,
        beta: Vec<f64>,
        ln_s: Vec<f64>, // [beta][alpha]
    },
}

impl Model {
    fn kt(&self) -> f64 {
        match self {
            Model::FreeGas { kt, .. } | Model::Tabulated { kt, .. } => *kt,
        }
    }

    fn mass(&self) -> f64 {
        match self {
            Model::FreeGas { awr, .. } => *awr,
            Model::Tabulated { a0, .. } => *a0,
        }
    }

    fn sigma_bound(&self) -> f64 {
        match self {
            Model::FreeGas { sigma_bound, .. } | Model::Tabulated { sigma_bound, .. } => {
                *sigma_bound
            }
        }
    }

    /// exp(-beta/2) S(alpha, beta) (mirrors `s_asym`).
    pub fn s_asym(&self, alpha: f64, beta: f64) -> f64 {
        match self {
            Model::FreeGas { .. } => {
                let a = alpha.max(1e-300);
                (-(a + beta).powi(2) / (4.0 * a)).exp() / (SQRT_4PI * a.sqrt())
            }
            Model::Tabulated {
                a0,
                scale,
                t_ratio,
                alpha: ga,
                beta: gb,
                ln_s,
                ..
            } => {
                let (na, nb) = (ga.len(), gb.len());
                let s = |ib: usize, ia: usize| ln_s[ib * na + ia];
                let a = alpha * scale;
                let b = beta.abs() * scale;
                let mut inside = a <= ga[na - 1] && b <= gb[nb - 1];
                let ac = a.clamp(ga[0], ga[na - 1]);
                let ia = ga
                    .partition_point(|&v| v <= ac)
                    .saturating_sub(1)
                    .min(na - 2);
                let ib = gb
                    .partition_point(|&v| v <= b)
                    .saturating_sub(1)
                    .min(nb - 2);
                let ta = (ac - ga[ia]) / (ga[ia + 1] - ga[ia]);
                let tb = ((b - gb[ib]) / (gb[ib + 1] - gb[ib])).clamp(0.0, 1.0);
                let mut ln = (1.0 - tb) * ((1.0 - ta) * s(ib, ia) + ta * s(ib, ia + 1))
                    + tb * ((1.0 - ta) * s(ib + 1, ia) + ta * s(ib + 1, ia + 1));
                let c_liq = if s(0, 0) > s(0, 1) {
                    (s(0, 0) - s(1, 0)) * ga[0] / (gb[1] * gb[1])
                } else {
                    0.0
                };
                let small = a < ga[0] && b <= 0.2 && c_liq > 0.0;
                if small {
                    let am = a.max(1e-300);
                    ln = s(0, 0) + 0.5 * (ga[0] / am).ln() - c_liq * b * b / am;
                }
                let tab = if ln <= -225.0 {
                    0.0
                } else {
                    (ln.max(-700.0) - beta / 2.0).exp()
                };
                let corners = s(ib, ia)
                    .min(s(ib, ia + 1))
                    .min(s(ib + 1, ia).min(s(ib + 1, ia + 1)));
                let exhausted = (a * a0 >= 30.0 || b >= 30.0) && corners <= -225.0 && !small;
                inside = inside && !exhausted;
                if inside {
                    return tab;
                }
                let am = alpha.max(1e-300);
                let arg =
                    (am - beta.abs()).powi(2) / (4.0 * am * t_ratio) + (beta.abs() + beta) / 2.0;
                (-arg).exp() / (SQRT_4PI * (am * t_ratio).sqrt())
            }
        }
    }
}

/// Gauss nodes and weights in E' (mirrors `_outgoing_nodes`).
fn outgoing_nodes(model: &Model, energy: f64, lo: f64, boundaries: &[f64]) -> (Vec<f64>, Vec<f64>) {
    let kt = model.kt();
    let hi = energy + 40.0 * kt;
    let width = 2.0 * (energy * kt / model.mass()).sqrt();
    let mut points = vec![lo, hi, energy];
    points.extend(boundaries.iter().copied().filter(|&b| lo < b && b < hi));
    let mut scales = vec![width];
    scales.extend((-4..2).map(|k| kt * 10f64.powi(k)));
    for scale in scales {
        for c in [0.1, 0.3, 1.0, 3.0, 10.0] {
            for e in [energy - c * scale, energy + c * scale] {
                if lo < e && e < hi {
                    points.push(e);
                }
            }
        }
    }
    let decades = (hi / lo).log10();
    let n = (60.0 * decades) as usize;
    points.extend((0..=n).map(|k| lo * (hi / lo).powf(k as f64 / (60.0 * decades))));
    points.sort_by(f64::total_cmp);
    points.dedup();
    let mut nodes = Vec::with_capacity(4 * points.len());
    let mut weights = Vec::with_capacity(4 * points.len());
    for w in points.windows(2) {
        let (a, b) = (w[0], w[1]);
        for (t, wt) in NODES_4.iter().zip(WEIGHTS_4.iter()) {
            nodes.push(0.5 * (a + b) + 0.5 * (b - a) * t);
            weights.push(0.5 * (b - a) * wt);
        }
    }
    (nodes, weights)
}

/// K_l(E -> g') flattened [l][g'], and the total scattering cross section
/// (mirrors `thermal_kernel`).
pub fn thermal_kernel(
    model: &Model,
    boundaries: &[f64],
    lmax: usize,
    energy: f64,
) -> (Vec<f64>, f64) {
    let g_count = boundaries.len() - 1;
    let (kt, a) = (model.kt(), model.mass());
    let lo = boundaries[g_count].min(energy) * 1e-3;
    let (e_out, w_out) = outgoing_nodes(model, energy, lo, boundaries);
    let root = (a * kt).sqrt();
    let s_min: Vec<f64> = e_out
        .iter()
        .map(|&e| (energy.sqrt() - e.sqrt()).abs() / root)
        .collect();
    let s_max: Vec<f64> = e_out
        .iter()
        .map(|&e| (energy.sqrt() + e.sqrt()) / root)
        .collect();
    let widest = s_min
        .iter()
        .zip(&s_max)
        .map(|(lo, hi)| hi - lo)
        .fold(0.0, f64::max);
    let n_pieces = (widest.ceil() as usize).max(1);
    let pref = model.sigma_bound() * a / (2.0 * energy);
    let mut k = vec![0.0; (lmax + 1) * g_count];
    let mut total = 0.0;
    for (n, (&ep, &wo)) in e_out.iter().zip(&w_out).enumerate() {
        let beta = (ep - energy) / kt;
        let piece = (s_max[n] - s_min[n]) / n_pieces as f64;
        let mut acc = vec![0.0; lmax + 1];
        for p in 0..n_pieces {
            let start = s_min[n] + piece * p as f64;
            for (t, w) in NODES_8.iter().zip(WEIGHTS_8.iter()) {
                let s = start + 0.5 * piece * (t + 1.0);
                let alpha = s * s;
                let mu = (energy + ep - a * kt * alpha) / (2.0 * (energy * ep).sqrt());
                let f = s * model.s_asym(alpha, beta) * (0.5 * piece * w);
                for (l, pl) in legendre(lmax, mu.clamp(-1.0, 1.0)).iter().enumerate() {
                    acc[l] += pl * f;
                }
            }
        }
        total += pref * acc[0] * wo;
        if let Some(g) = find_group(boundaries, ep) {
            for (l, v) in acc.iter().enumerate() {
                k[l * g_count + g] += pref * v * wo;
            }
        }
    }
    (k, total)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// (1 + 1/2y^2) erf(y) + exp(-y^2)/(y sqrt(pi)) for y = 1 (chapter 8)
    #[test]
    fn free_gas_total_matches_closed_form() {
        let (awr, kt): (f64, f64) = (15.858, 0.0253);
        let energy = awr.recip() * kt; // y = 1
        let model = Model::FreeGas {
            awr,
            kt,
            sigma_bound: ((awr + 1.0) / awr).powi(2),
        };
        let b: Vec<f64> = (0..=60)
            .map(|i| 10.0 * (1e-6f64).powf(i as f64 / 60.0))
            .collect();
        let (_, total) = thermal_kernel(&model, &b, 0, energy);
        let exact = 1.5 * 0.8427007929497149 + (-1f64).exp() / std::f64::consts::PI.sqrt();
        let expected = exact * model.sigma_bound() * (awr / (awr + 1.0)).powi(2);
        assert!(
            (total / expected - 1.0).abs() < 1e-4,
            "{total} vs {expected}"
        );
    }
}
