//! Ultra-fine slowing-down march (port of `fission.pincell.shielding.march`).
//!
//! Flux per unit lethargy psi_r(u) in R regions (each with its own mixture), marched down
//! in energy; the elastic source of every nuclide is its window integral
//! C_i(u) - C_i(u - Delta_i) of f_i = Sigma_s,i psi e^u (trapezoid rule, the
//! same at both ends of the window). One region uses the closed form
//! psi = ((1 - P) known + P Sigma_t psi_M) / (Sigma_t - (1 - P) coef); several
//! regions solve Sigma_t,i V_i psi_i - sum_j V_j P_ji coef_j psi_j =
//! sum_j V_j P_ji known_j + Sigma_t,i V_i E_i psi_M at every point.
//!
//! Nuclides with a moving-target kernel (`moving.rs`) scatter through it below
//! its top energy: S += E sum_m P(E_m -> E) Phi_m, with the hat moments Phi_m
//! of Sigma_s psi. Upscattering makes the source depend on psi below u, so the
//! march is repeated (Gauss-Seidel, first pass at rest) until psi converges.

pub struct Kernel<'a> {
    pub scale: f64,
    pub base: f64,
    pub dx: f64,
    pub top: usize,
    pub nodes: usize,
    pub starts: &'a [i64],
    pub rows: &'a [i64],
    pub values: &'a [f64],
}

pub struct Input<'a> {
    pub energies: &'a [f64],   // (n), descending
    pub total: &'a [f64],      // (n, R)
    pub scatter: &'a [f64],    // (nuclides, n, R)
    pub alphas: &'a [f64],     // (nuclides)
    pub escape: &'a [f64],     // (n, R)
    pub collide: &'a [f64],    // (n, R, R): P_j->i at [k, j, i]; empty for one region
    pub volumes: &'a [f64],    // (R)
    pub history: &'a [f64],    // (start, R)
    pub start: usize,
    pub share: f64,
    pub outside: &'a [f64],    // (n) or empty
    pub kernels: Vec<Option<Kernel<'a>>>, // per nuclide, or empty
    pub tol: f64,
    pub max_passes: usize,
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

/// Where every march point sits on a kernel's mesh (Python `_Targets`).
struct Targets<'a> {
    kernel: &'a Kernel<'a>,
    p: Vec<f64>,
    on: Vec<f64>,
    below: Vec<bool>,
    cell: Vec<usize>,
}

impl<'a> Targets<'a> {
    fn new(kernel: &'a Kernel<'a>, energies: &[f64], start: usize) -> Result<Self, String> {
        let n = energies.len();
        let x_top = kernel.base + kernel.dx * kernel.top as f64;
        let e_top = x_top * x_top / kernel.scale;
        let p: Vec<f64> = energies.iter().map(|&e| ((kernel.scale * e).sqrt() - kernel.base) / kernel.dx).collect();
        let top = e_top * (1.0 + 1e-12);
        let on: Vec<f64> = energies.iter().map(|&e| if e >= e_top * (1.0 - 1e-12) { 1.0 } else { 0.0 }).collect();
        let mut below = vec![false; n];
        for k in 1..n {
            below[k] = k >= start && energies[k - 1] <= top;
        }
        let cell: Vec<usize> = p
            .iter()
            .map(|&q| ((q + 1e-9).floor() as i64).clamp(0, kernel.nodes as i64 - 2) as usize)
            .collect();
        for k in 0..n {
            if below[k] && p[k - 1] - cell[k] as f64 > 1.0 + 1e-6 {
                return Err("the kernel's mesh energies must be grid points".into());
            }
        }
        Ok(Targets { kernel, p, on, below, cell })
    }

    /// s(E_k) per region (Python `_moving_source`), with the per-pass cache.
    fn source(&self, k: usize, phi_new: &[f64], phi_old: &[f64], cache: &mut [f64], done_flags: &mut [bool],
              regions: usize, out: &mut [f64]) {  // fmt: skip
        out.iter_mut().for_each(|v| *v = 0.0);
        let kernel = self.kernel;
        let p = self.p[k];
        if p < 0.0 || p > (kernel.nodes - 1) as f64 {
            return;
        }
        let j = self.cell[k];
        let theta = (p - j as f64).max(0.0);
        let above = self.p[k - 1];
        let q = (above + 1e-9).floor() as i64;
        let first = if (above - above.round()).abs() < 1e-9 { q + 1 } else { q + 2 };
        for (node, weight) in [(j, 1.0 - theta), (j + 1, theta)] {
            if !done_flags[node] {
                let (a, b) = (kernel.starts[node] as usize, kernel.starts[node + 1] as usize);
                for r in 0..regions {
                    let mut sum = 0.0;
                    for e in a..b {
                        let m = kernel.rows[e] as usize;
                        let phi = if kernel.rows[e] >= first { phi_new[m * regions + r] } else { phi_old[m * regions + r] };
                        sum += kernel.values[e] * phi;
                    }
                    cache[node * regions + r] = sum;
                }
                done_flags[node] = true;
            }
            for r in 0..regions {
                out[r] += weight * cache[node * regions + r];
            }
        }
    }
}

/// Hat moments below the last grid point with Sigma_s psi continued flat
/// (Python `_below_the_range`).
fn below_the_range(tg: &Targets, e0: f64, g: &[f64], phi: &mut [f64], regions: usize) {
    let kernel = tg.kernel;
    let mut top = tg.p[tg.p.len() - 1];
    let mut node = (top + 1e-9).floor() as i64;
    if (top - top.round()).abs() < 1e-9 {
        node = top.round() as i64 - 1;
    }
    while node >= 0 {
        let e_hi = (kernel.base + top * kernel.dx).powi(2) / kernel.scale;
        let e_lo = (kernel.base + node as f64 * kernel.dx).powi(2) / kernel.scale;
        let h = (e0 / e_lo).ln() - (e0 / e_hi).ln();
        let t0 = top - node as f64;
        let nd = node as usize;
        for r in 0..regions {
            phi[nd * regions + r] += 0.5 * h * ((1.0 - t0) * g[r] + g[r]);
            phi[(nd + 1) * regions + r] += 0.5 * h * (t0 * g[r]);
        }
        top = node as f64;
        node -= 1;
    }
}

struct PassOut {
    psi: Vec<f64>,
    levels: Vec<f64>,
    phi: Vec<Vec<f64>>,
}

fn march_pass(inp: &Input, targets: &[Option<Targets>], phi_old: &[Vec<f64>], rest: bool) -> PassOut {
    let n = inp.energies.len();
    let regions = inp.volumes.len();
    let m = inp.alphas.len();
    let single = inp.collide.is_empty();
    let start = inp.start;
    let u: Vec<f64> = inp.energies.iter().map(|e| (inp.energies[0] / e).ln()).collect();
    let deltas: Vec<f64> = inp.alphas.iter().map(|a| (1.0 / a).ln()).collect();
    let factors: Vec<f64> = inp.alphas.iter().map(|a| 1.0 / (1.0 - a)).collect();
    let vsum: f64 = inp.volumes.iter().sum();
    let weights: Vec<f64> = inp.volumes.iter().map(|v| v / vsum).collect();
    let sc = |i: usize, k: usize, r: usize| inp.scatter[(i * n + k) * regions + r];
    let tot = |k: usize, r: usize| inp.total[k * regions + r];
    let ones = vec![1.0; n];
    let on: Vec<&[f64]> = targets
        .iter()
        .map(|t| match t {
            Some(t) if !rest => t.on.as_slice(),
            _ => ones.as_slice(),
        })
        .collect();

    let mut psi = vec![0.0; n * regions];
    psi[..start * regions].copy_from_slice(inp.history);
    // f[i][k, r] and cumulative C[i][k, r]
    let mut f = vec![0.0; m * n * regions];
    let mut cum = vec![0.0; m * n * regions];
    for i in 0..m {
        for k in 0..start {
            let eu = u[k].exp();
            for r in 0..regions {
                f[(i * n + k) * regions + r] = sc(i, k, r) * psi[k * regions + r] * eu;
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
    let mut phi_new: Vec<Vec<f64>> = phi_old.iter().map(|p| vec![0.0; p.len()]).collect();
    let mut caches: Vec<Vec<f64>> = phi_old.iter().map(|p| vec![0.0; p.len()]).collect();
    let mut cached: Vec<Vec<bool>> = targets
        .iter()
        .map(|t| t.as_ref().map_or(Vec::new(), |t| vec![false; t.kernel.nodes]))
        .collect();
    let mut moving = vec![0.0; regions];
    let mut levels = vec![1.0; n];
    let mut pointer = vec![0usize; m];
    let mut level = 1.0;
    let mut known = vec![0.0; regions];
    let mut coefficient = vec![0.0; regions];
    let mut a = vec![0.0; regions * regions];
    let mut b = vec![0.0; regions];
    for k in start..n {
        let h = u[k] - u[k - 1];
        let eu = u[k].exp();
        known.iter_mut().for_each(|x| *x = 0.0);
        coefficient.iter_mut().for_each(|x| *x = 0.0);
        for i in 0..m {
            let lower = u[k] - deltas[i];
            let mut j = pointer[i];
            while u[j + 1] < lower {
                j += 1;
            }
            pointer[i] = j;
            let o = on[i];
            let fi = |kk: usize, r: usize| f[(i * n + kk) * regions + r];
            let ci = |kk: usize, r: usize| cum[(i * n + kk) * regions + r];
            let weight_k;
            if j + 1 < k {
                let t = (lower - u[j]) / (u[j + 1] - u[j]);
                for r in 0..regions {
                    let g_lower = fi(j, r) + t * (fi(j + 1, r) - fi(j, r));
                    let c_lower = ci(j, r) + o[j + 1] * 0.5 * (lower - u[j]) * (fi(j, r) + g_lower);
                    let window = ci(k - 1, r) + o[k] * 0.5 * h * fi(k - 1, r) - c_lower;
                    known[r] += factors[i] * window / eu;
                }
                weight_k = o[k] * 0.5 * h;
            } else {
                let t = (lower - u[k - 1]) / h;
                let width = u[k] - lower;
                for r in 0..regions {
                    let window = o[k] * 0.5 * width * (1.0 - t) * fi(k - 1, r);
                    known[r] += factors[i] * window / eu;
                }
                weight_k = o[k] * 0.5 * width * (1.0 + t);
            }
            for r in 0..regions {
                coefficient[r] += factors[i] * weight_k * sc(i, k, r);
            }
            if let Some(tg) = &targets[i] {
                if !rest {
                    tg.source(k, &phi_new[i], &phi_old[i], &mut caches[i], &mut cached[i], regions, &mut moving);
                    for r in 0..regions {
                        known[r] += inp.energies[k] * moving[r];
                    }
                }
            }
        }
        if inp.share > 0.0 || !inp.outside.is_empty() {
            let mut mean = 0.0;
            for r in 0..regions {
                let mut scatter_prev = 0.0;
                for i in 0..m {
                    scatter_prev += sc(i, k - 1, r);
                }
                mean += weights[r] * (tot(k - 1, r) - scatter_prev) * psi[(k - 1) * regions + r];
            }
            let mut absorbed = inp.share * mean;
            if !inp.outside.is_empty() {
                absorbed += inp.outside[k - 1] * level;
            }
            level -= h * absorbed;
        }
        levels[k] = level;
        if single {
            let (p, t) = (inp.escape[k], tot(k, 0));
            psi[k] = ((1.0 - p) * known[0] + p * t * level) / (t - (1.0 - p) * coefficient[0]);
        } else {
            let v = inp.volumes;
            for i in 0..regions {
                let t = tot(k, i);
                let mut rhs = t * v[i] * inp.escape[k * regions + i] * level;
                for j in 0..regions {
                    // V_j P_j->i
                    let transfer = inp.collide[(k * regions + j) * regions + i] * v[j];
                    a[i * regions + j] = -coefficient[j] * transfer;
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
                f[idx] = sc(i, k, r) * psi[k * regions + r] * eu;
                cum[idx] = cum[prev] + on[i][k] * 0.5 * h * (f[prev] + f[idx]);
            }
            if let Some(tg) = &targets[i] {
                if tg.below[k] {
                    let j = tg.cell[k];
                    let (t0, t1) = (tg.p[k - 1] - j as f64, tg.p[k] - j as f64);
                    for r in 0..regions {
                        let g0 = sc(i, k - 1, r) * psi[(k - 1) * regions + r];
                        let g1 = sc(i, k, r) * psi[k * regions + r];
                        phi_new[i][j * regions + r] += 0.5 * h * ((1.0 - t0) * g0 + (1.0 - t1) * g1);
                        phi_new[i][(j + 1) * regions + r] += 0.5 * h * (t0 * g0 + t1 * g1);
                    }
                }
            }
        }
    }
    for (i, tg) in targets.iter().enumerate() {
        if let Some(tg) = tg {
            let g: Vec<f64> = (0..regions).map(|r| sc(i, n - 1, r) * psi[(n - 1) * regions + r]).collect();
            below_the_range(tg, inp.energies[0], &g, &mut phi_new[i], regions);
        }
    }
    for l in levels.iter_mut().take(start) {
        *l = 1.0;
    }
    PassOut { psi, levels, phi: phi_new }
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
    if inp.total.len() != n * regions || inp.scatter.len() != m * n * regions {
        return Err("total must be (n, R) and scatter (nuclides, n, R)".into());
    }
    if !single && inp.collide.len() != n * regions * regions {
        return Err("collide must be (n, R, R)".into());
    }
    if !inp.kernels.is_empty() && inp.kernels.len() != m {
        return Err("one kernel (or None) per nuclide".into());
    }
    let u_hist = (inp.energies[0] / inp.energies[start - 1]).ln();
    let widest = inp.alphas.iter().map(|a| (1.0 / a).ln()).fold(f64::MIN, f64::max);
    if u_hist < widest {
        return Err("the history must span the widest scattering window (ln 1/alpha)".into());
    }
    let mut targets: Vec<Option<Targets>> = Vec::with_capacity(m);
    for i in 0..m {
        targets.push(match inp.kernels.get(i) {
            Some(Some(k)) => Some(Targets::new(k, inp.energies, start)?),
            _ => None,
        });
    }
    let any = targets.iter().any(|t| t.is_some());
    let mut phi: Vec<Vec<f64>> = targets
        .iter()
        .map(|t| t.as_ref().map_or(Vec::new(), |t| vec![0.0; t.kernel.nodes * regions]))
        .collect();
    let first = march_pass(inp, &targets, &phi, true);
    let (mut psi, mut levels) = (first.psi, first.levels);
    phi = first.phi;
    if !any {
        return Ok((psi, levels));
    }
    for _ in 0..inp.max_passes {
        let next = march_pass(inp, &targets, &phi, false);
        let mut change: f64 = 0.0;
        for (new, old) in next.psi.iter().zip(psi.iter()) {
            change = change.max((new - old).abs() / new.abs().max(1e-300));
        }
        psi = next.psi;
        levels = next.levels;
        phi = next.phi;
        if change < inp.tol {
            return Ok((psi, levels));
        }
    }
    Err(format!("moving-target iterations: no convergence in {}", inp.max_passes))
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
            share: 0.0, outside: &[], kernels: Vec::new(), tol: 1e-10, max_passes: 10,
        })
        .unwrap();
        let worst = psi.iter().map(|p| (p - 1.0).abs()).fold(0.0, f64::max);
        assert!(worst < 1e-5, "worst {worst}");
    }
}
