//! Ultra-fine slowing-down march (port of `fission.pincell.shielding.march`).
//!
//! Flux per unit lethargy psi_r(u) in R regions of one mixture, marched down
//! in energy; the elastic source of every nuclide is its window integral
//! C_i(u) - C_i(u - Delta_i) of f_i = Sigma_s,i psi e^u (trapezoid rule, the
//! same at both ends of the window). One region uses the closed form
//! psi = ((1 - P) known + P Sigma_t psi_M) / (Sigma_t - (1 - P) coef); several
//! regions solve Sigma_t V_i psi_i - coef sum_j V_j P_ji psi_j =
//! sum_j V_j P_ji known_j + Sigma_t V_i E_i psi_M at every point.

pub struct Input<'a> {
    pub energies: &'a [f64],   // (n), descending
    pub total: &'a [f64],      // (n)
    pub scatter: &'a [f64],    // (nuclides, n)
    pub alphas: &'a [f64],     // (nuclides)
    pub escape: &'a [f64],     // (n, R)
    pub collide: &'a [f64],    // (n, R, R): P_j->i at [k, j, i]; empty for one region
    pub volumes: &'a [f64],    // (R)
    pub history: &'a [f64],    // (start, R)
    pub start: usize,
    pub share: f64,
    pub outside: &'a [f64],    // (n) or empty
}

/// Solve a x = b in place (Gaussian elimination, partial pivoting); r is small.
fn solve_small(a: &mut [f64], b: &mut [f64], r: usize) {
    for col in 0..r {
        let mut pivot = col;
        for row in col + 1..r {
            if a[row * r + col].abs() > a[pivot * r + col].abs() {
                pivot = row;
            }
        }
        if pivot != col {
            for c in 0..r {
                a.swap(col * r + c, pivot * r + c);
            }
            b.swap(col, pivot);
        }
        for row in col + 1..r {
            let factor = a[row * r + col] / a[col * r + col];
            for c in col..r {
                a[row * r + c] -= factor * a[col * r + c];
            }
            b[row] -= factor * b[col];
        }
    }
    for row in (0..r).rev() {
        let mut sum = b[row];
        for c in row + 1..r {
            sum -= a[row * r + c] * b[c];
        }
        b[row] = sum / a[row * r + row];
    }
}

/// Returns (psi (n, R), psi_M (n)).
pub fn march(inp: &Input) -> Result<(Vec<f64>, Vec<f64>), String> {
    let n = inp.energies.len();
    let regions = inp.volumes.len();
    let m = inp.alphas.len();
    let single = inp.collide.is_empty();
    let start = inp.start;
    if inp.escape.len() != n * regions || inp.history.len() != start * regions {
        return Err("escape / history shapes do not match".into());
    }
    if !single && inp.collide.len() != n * regions * regions {
        return Err("collide must be (n, R, R)".into());
    }
    let u: Vec<f64> = inp.energies.iter().map(|e| (inp.energies[0] / e).ln()).collect();
    let deltas: Vec<f64> = inp.alphas.iter().map(|a| (1.0 / a).ln()).collect();
    let widest = deltas.iter().cloned().fold(f64::MIN, f64::max);
    if u[start - 1] < widest {
        return Err("the history must span the widest scattering window (ln 1/alpha)".into());
    }
    let factors: Vec<f64> = inp.alphas.iter().map(|a| 1.0 / (1.0 - a)).collect();
    let vsum: f64 = inp.volumes.iter().sum();
    let weights: Vec<f64> = inp.volumes.iter().map(|v| v / vsum).collect();
    let sc = |i: usize, k: usize| inp.scatter[i * n + k];

    let mut psi = vec![0.0; n * regions];
    psi[..start * regions].copy_from_slice(inp.history);
    // f[i][k, r] and cumulative C[i][k, r]
    let mut f = vec![0.0; m * n * regions];
    let mut cum = vec![0.0; m * n * regions];
    for i in 0..m {
        for k in 0..start {
            let eu = u[k].exp();
            for r in 0..regions {
                f[(i * n + k) * regions + r] = sc(i, k) * psi[k * regions + r] * eu;
            }
        }
        for k in 1..start {
            let h = u[k] - u[k - 1];
            for r in 0..regions {
                let idx = (i * n + k) * regions + r;
                let prev = (i * n + k - 1) * regions + r;
                cum[idx] = cum[prev] + 0.5 * h * (f[idx] + f[prev]);
            }
        }
    }
    let mut levels = vec![1.0; n];
    let mut pointer = vec![0usize; m];
    let mut level = 1.0;
    let mut known = vec![0.0; regions];
    let mut a = vec![0.0; regions * regions];
    let mut b = vec![0.0; regions];
    for k in start..n {
        let h = u[k] - u[k - 1];
        let eu = u[k].exp();
        known.iter_mut().for_each(|x| *x = 0.0);
        let mut coefficient = 0.0;
        for i in 0..m {
            let lower = u[k] - deltas[i];
            let mut j = pointer[i];
            while u[j + 1] < lower {
                j += 1;
            }
            pointer[i] = j;
            let fi = |kk: usize, r: usize| f[(i * n + kk) * regions + r];
            let ci = |kk: usize, r: usize| cum[(i * n + kk) * regions + r];
            let weight_k;
            if j + 1 < k {
                let t = (lower - u[j]) / (u[j + 1] - u[j]);
                for r in 0..regions {
                    let g_lower = fi(j, r) + t * (fi(j + 1, r) - fi(j, r));
                    let c_lower = ci(j, r) + 0.5 * (lower - u[j]) * (fi(j, r) + g_lower);
                    let window = ci(k - 1, r) + 0.5 * h * fi(k - 1, r) - c_lower;
                    known[r] += factors[i] * window / eu;
                }
                weight_k = 0.5 * h;
            } else {
                let t = (lower - u[k - 1]) / h;
                let width = u[k] - lower;
                for r in 0..regions {
                    let window = 0.5 * width * (1.0 - t) * fi(k - 1, r);
                    known[r] += factors[i] * window / eu;
                }
                weight_k = 0.5 * width * (1.0 + t);
            }
            coefficient += factors[i] * weight_k * sc(i, k);
        }
        if inp.share > 0.0 || !inp.outside.is_empty() {
            let mut scatter_prev = 0.0;
            for i in 0..m {
                scatter_prev += sc(i, k - 1);
            }
            let mut mean = 0.0;
            for r in 0..regions {
                mean += weights[r] * psi[(k - 1) * regions + r];
            }
            let mut absorbed = inp.share * (inp.total[k - 1] - scatter_prev) * mean;
            if !inp.outside.is_empty() {
                absorbed += inp.outside[k - 1] * level;
            }
            level -= h * absorbed;
        }
        levels[k] = level;
        let t = inp.total[k];
        if single {
            let p = inp.escape[k];
            psi[k] = ((1.0 - p) * known[0] + p * t * level) / (t - (1.0 - p) * coefficient);
        } else {
            let v = inp.volumes;
            for i in 0..regions {
                let mut rhs = t * v[i] * inp.escape[k * regions + i] * level;
                for j in 0..regions {
                    // V_j P_j->i
                    let transfer = inp.collide[(k * regions + j) * regions + i] * v[j];
                    a[i * regions + j] = -coefficient * transfer;
                    rhs += transfer * known[j];
                }
                a[i * regions + i] += t * v[i];
                b[i] = rhs;
            }
            solve_small(&mut a, &mut b, regions);
            psi[k * regions..(k + 1) * regions].copy_from_slice(&b);
        }
        for i in 0..m {
            for r in 0..regions {
                let idx = (i * n + k) * regions + r;
                let prev = (i * n + k - 1) * regions + r;
                f[idx] = sc(i, k) * psi[k * regions + r] * eu;
                cum[idx] = cum[prev] + 0.5 * h * (f[prev] + f[idx]);
            }
        }
    }
    for l in levels.iter_mut().take(start) {
        *l = 1.0;
    }
    Ok((psi, levels))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn small_solve_matches_known_solution() {
        let mut a = vec![2.0, 1.0, 0.0, 1.0, 3.0, 1.0, 0.0, 1.0, 4.0];
        let x = [1.0, -2.0, 0.5];
        let mut b: Vec<f64> = (0..3).map(|i| (0..3).map(|j| a[i * 3 + j] * x[j]).sum()).collect();
        solve_small(&mut a, &mut b, 3);
        for i in 0..3 {
            assert!((b[i] - x[i]).abs() < 1e-14);
        }
    }

    #[test]
    fn pure_scattering_keeps_the_asymptotic_flux() {
        let n = 60000; // as the Python test: the trapezoid error goes as h^2
        let energies: Vec<f64> = (0..n).map(|k| 2e3 * (1e-3f64).powf(k as f64 / (n - 1) as f64)).collect();
        let scatter = vec![0.3; n];
        let total = scatter.clone();
        let escape = vec![0.0; n];
        let u_hist = (energies[0] / energies[1]).ln();
        let start = (0.3 / u_hist).ceil() as usize + 2;
        let history = vec![1.0; start];
        let (psi, _) = march(&Input {
            energies: &energies, total: &total, scatter: &scatter, alphas: &[0.983],
            escape: &escape, collide: &[], volumes: &[1.0], history: &history, start,
            share: 0.0, outside: &[],
        })
        .unwrap();
        let worst = psi.iter().map(|p| (p - 1.0).abs()).fold(0.0, f64::max);
        assert!(worst < 1e-5, "worst {worst}");
    }
}
