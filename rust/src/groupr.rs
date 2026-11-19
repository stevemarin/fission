//! Group averaging: the hat integrals A[g, j] = Int_g sigma_c w h_j dE.
//! Port of `hat_integrals` in `fission/processing/groupr.py`.

use std::collections::BTreeMap;

use crate::doppler::{NODES_8, WEIGHTS_8};

/// Weighting spectrum (mirrors `processing/weights.py`).
#[derive(Clone, Copy, Debug)]
pub enum Weight {
    Constant,
    InverseE,
    /// Maxwellian (E <= eb) -- 1/E -- fission spectrum (E > ec): NJOY iwt = 4.
    ThermalSlowingDownFission {
        eb: f64,
        kt: f64,
        ab: f64,
        ec: f64,
        tf: f64,
        ac: f64,
    },
}

impl Weight {
    /// From [kind, parameters...] as built by `_weight_params`.
    pub fn from_params(p: &[f64]) -> Option<Self> {
        match (p.first().copied()? as i64, p.len()) {
            (0, _) => Some(Weight::Constant),
            (1, _) => Some(Weight::InverseE),
            (2, 7) => Some(Weight::ThermalSlowingDownFission {
                eb: p[1],
                kt: p[2],
                ab: p[3],
                ec: p[4],
                tf: p[5],
                ac: p[6],
            }),
            _ => None,
        }
    }

    pub fn breakpoints(&self) -> Vec<f64> {
        match *self {
            Weight::ThermalSlowingDownFission { eb, ec, .. } => vec![eb, ec],
            _ => vec![],
        }
    }

    pub fn eval(&self, e: f64) -> f64 {
        match *self {
            Weight::Constant => 1.0,
            Weight::InverseE => 1.0 / e,
            Weight::ThermalSlowingDownFission {
                eb,
                kt,
                ab,
                ec,
                tf,
                ac,
            } => {
                if e <= eb {
                    ab * e * (-e / kt).exp()
                } else if e <= ec {
                    1.0 / e
                } else {
                    ac * e.sqrt() * (-e / tf).exp()
                }
            }
        }
    }
}

/// Index i with v[i] <= x < v[i+1] (Python `bisect_right(v, x) - 1`).
pub fn bisect_right_minus_one(v: &[f64], x: f64) -> usize {
    v.partition_point(|&a| a <= x).saturating_sub(1)
}

/// Group containing e for descending boundaries (mirrors `GroupStructure.find`).
pub fn find_group(boundaries: &[f64], e: f64) -> Option<usize> {
    let n = boundaries.len() - 1;
    if !(boundaries[n] <= e && e <= boundaries[0]) {
        return None;
    }
    // ascending position, as bisect_right on the reversed list
    let below = boundaries.iter().filter(|&&b| b <= e).count(); // boundaries <= e
    let i = (below - 1).min(n - 1);
    Some(n - 1 - i)
}

/// A[(g, j)][c] for linear sigma columns on x and hat functions on z.
pub fn hat_integrals(
    x: &[f64],
    columns: &[Vec<f64>],
    weight: Weight,
    boundaries: &[f64],
    z: &[f64],
    max_ratio: f64,
) -> BTreeMap<(usize, usize), Vec<f64>> {
    let lo = x[0].max(z[0]).max(boundaries[boundaries.len() - 1]);
    let hi = x[x.len() - 1].min(z[z.len() - 1]).min(boundaries[0]);
    let mut points: Vec<f64> = x
        .iter()
        .chain(z)
        .chain(boundaries)
        .copied()
        .chain(weight.breakpoints())
        .filter(|&e| lo <= e && e <= hi)
        .collect();
    points.sort_by(f64::total_cmp);
    points.dedup();

    let mut acc: BTreeMap<(usize, usize), Vec<f64>> = BTreeMap::new();
    let nc = columns.len();
    for w in points.windows(2) {
        let (a, b) = (w[0], w[1]);
        if b <= a {
            continue;
        }
        let n = if a > 0.0 {
            (((b / a).ln() / max_ratio.ln()).ceil() as usize).max(1)
        } else {
            1
        };
        for k in 0..n {
            let (pa, pb) = if a > 0.0 {
                let r = b / a;
                let pa = a * r.powf(k as f64 / n as f64);
                let pb = if k == n - 1 {
                    b
                } else {
                    a * r.powf((k + 1) as f64 / n as f64)
                };
                (pa, pb)
            } else {
                let pa = a + (b - a) * k as f64 / n as f64;
                let pb = if k == n - 1 {
                    b
                } else {
                    a + (b - a) * (k + 1) as f64 / n as f64
                };
                (pa, pb)
            };
            piece(x, columns, weight, boundaries, z, pa, pb, nc, &mut acc);
        }
    }
    acc
}

#[allow(clippy::too_many_arguments)]
fn piece(
    x: &[f64],
    columns: &[Vec<f64>],
    weight: Weight,
    boundaries: &[f64],
    z: &[f64],
    a: f64,
    b: f64,
    nc: usize,
    acc: &mut BTreeMap<(usize, usize), Vec<f64>>,
) {
    let mid = 0.5 * (a + b);
    let i = bisect_right_minus_one(x, mid).min(x.len() - 2);
    let j = bisect_right_minus_one(z, mid).min(z.len() - 2);
    let Some(g) = find_group(boundaries, mid) else {
        return;
    };
    let (x0, x1, z0, z1) = (x[i], x[i + 1], z[j], z[j + 1]);
    let (mut m00, mut m01, mut m10, mut m11) = (0.0, 0.0, 0.0, 0.0);
    for (t, wt) in NODES_8.iter().zip(WEIGHTS_8.iter()) {
        let e = 0.5 * (a + b) + 0.5 * (b - a) * t;
        let f = wt * weight.eval(e);
        let u = (e - x0) / (x1 - x0);
        let v = (e - z0) / (z1 - z0);
        m00 += f * (1.0 - u) * (1.0 - v);
        m01 += f * (1.0 - u) * v;
        m10 += f * u * (1.0 - v);
        m11 += f * u * v;
    }
    let h = 0.5 * (b - a);
    {
        let left = acc.entry((g, j)).or_insert_with(|| vec![0.0; nc]);
        for (c, col) in columns.iter().enumerate() {
            left[c] += h * (col[i] * m00 + col[i + 1] * m10);
        }
    }
    let right = acc.entry((g, j + 1)).or_insert_with(|| vec![0.0; nc]);
    for (c, col) in columns.iter().enumerate() {
        right[c] += h * (col[i] * m01 + col[i + 1] * m11);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn flux_of_one_over_e_is_the_lethargy_width() {
        let b = [1e3, 1e2, 1.0];
        let a = hat_integrals(
            &[1.0, 1e3],
            &[vec![1.0, 1.0]],
            Weight::InverseE,
            &b,
            &[1.0, 1e3],
            1.2,
        );
        let g0: f64 = a
            .iter()
            .filter(|((g, _), _)| *g == 0)
            .map(|(_, v)| v[0])
            .sum();
        assert!((g0 - 10f64.ln()).abs() < 1e-12, "{g0}");
    }
}
