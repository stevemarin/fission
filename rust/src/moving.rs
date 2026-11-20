//! Free-gas elastic kernel with the 0 K cross section (port of
//! `fission.pincell.moving`): erfcx, the hat weights W_j of the relative-speed
//! grid, and the normalized kernel P(E_m -> E_n) on the uniform x mesh, stored by
//! outgoing node. Threaded over incoming nodes (each thread keeps its own window
//! of per-node prefix/suffix sums); the result is identical to the serial one.

use std::collections::HashMap;
use std::f64::consts::{LN_2, PI};

use crate::doppler::{BOLTZMANN_EV_PER_K, NODES_4 as GAUSS_4_NODES, WEIGHTS_4 as GAUSS_4_WEIGHTS};

const Z_CUTOFF: f64 = 6.0;
const TAIL: f64 = 2.0 * Z_CUTOFF + 1.0;
const ERFCX_SWITCH: f64 = 2.0;
const SERIES_TERMS: usize = 60;
const FRACTION_TERMS: usize = 80;

/// exp(y^2) erfc(y) for y >= 0, the same fixed-depth series / continued
/// fraction as the Python reference.
pub fn erfcx(y: f64) -> f64 {
    let sqrt_pi = PI.sqrt();
    if y < ERFCX_SWITCH {
        let t2 = 2.0 * y * y;
        let mut term = y;
        let mut total = y;
        for n in 1..SERIES_TERMS {
            term = term * t2 / (2 * n + 1) as f64;
            total += term;
        }
        (y * y).exp() - 2.0 / sqrt_pi * total
    } else {
        let mut t = y;
        for k in (1..=FRACTION_TERMS).rev() {
            t = y + 0.5 * k as f64 / t;
        }
        1.0 / (sqrt_pi * t)
    }
}

/// numpy's logaddexp.
fn log_add(x: f64, y: f64) -> f64 {
    if x == y {
        return x + LN_2;
    }
    let tmp = x - y;
    if tmp > 0.0 {
        x + (-tmp).exp().ln_1p()
    } else if tmp <= 0.0 {
        y + tmp.exp().ln_1p()
    } else {
        tmp
    }
}

/// numpy's interp (clamped at the ends).
fn interp(x: f64, xp: &[f64], fp: &[f64]) -> f64 {
    let n = xp.len();
    if x <= xp[0] {
        return fp[0];
    }
    if x >= xp[n - 1] {
        return fp[n - 1];
    }
    let j = xp.partition_point(|&v| v <= x) - 1;
    let slope = (fp[j + 1] - fp[j]) / (xp[j + 1] - xp[j]);
    slope * (x - xp[j]) + fp[j]
}

/// W_j = Int hat_j(z) z sigma_0(z^2 / scale) dz on z_j = z0 + j dz.
pub fn hat_weights(z0: f64, dz: f64, count: usize, energies: &[f64], sigma: &[f64], scale: f64) -> Vec<f64> {
    let z_hi = z0 + (count - 1) as f64 * dz;
    let mut points: Vec<f64> = (0..count).map(|j| z0 + dz * j as f64).collect();
    points.extend(energies.iter().map(|&e| (scale * e).sqrt()).filter(|&z| z > z0 && z < z_hi));
    points.sort_by(|a, b| a.partial_cmp(b).unwrap());
    points.dedup();
    let pieces = points.len() - 1;
    let cells: Vec<usize> = (0..pieces)
        .map(|i| {
            let c = ((0.5 * (points[i] + points[i + 1]) - z0) / dz).floor() as i64;
            c.min(count as i64 - 2) as usize
        })
        .collect();
    let mut w = vec![0.0; count];
    let mut left = vec![0.0; count];
    let mut right = vec![0.0; count];
    for (&t, &wt) in GAUSS_4_NODES.iter().zip(GAUSS_4_WEIGHTS.iter()) {
        left.iter_mut().for_each(|v| *v = 0.0);
        right.iter_mut().for_each(|v| *v = 0.0);
        for i in 0..pieces {
            let (lo, hi) = (points[i], points[i + 1]);
            let z = 0.5 * (lo + hi) + 0.5 * (hi - lo) * t;
            let f = 0.5 * (hi - lo) * wt * z * interp(z * z / scale, energies, sigma);
            let c = cells[i];
            let r = (z - (z0 + c as f64 * dz)) / dz;
            left[c] += f * (1.0 - r);
            right[c + 1] += f * r;
        }
        for j in 0..count {
            w[j] += left[j];
        }
        for j in 0..count {
            w[j] += right[j];
        }
    }
    w
}

pub struct Mesh {
    pub scale: f64,
    pub base: f64,
    pub top: usize,
    pub nodes: usize,
    pub z0: f64,
    pub count: usize,
}

pub fn mesh(awr: f64, temperature: f64, e_low: f64, e_top: f64, dx: f64, dz: f64) -> Result<Mesh, String> {
    let scale = awr / (BOLTZMANN_EV_PER_K * temperature);
    let sqrt_alpha = (awr - 1.0) / (awr + 1.0);
    let (x_top, x_low) = ((scale * e_top).sqrt(), (scale * e_low).sqrt());
    let below = ((x_top - (sqrt_alpha * (x_low - TAIL) - TAIL)) / dx).ceil() as usize;
    let above = (TAIL / dx).ceil() as usize;
    let base = x_top - below as f64 * dx;
    if base < 10.0 {
        return Err(format!(
            "x = {base:.1} at the bottom of the mesh: the upper limit of the \
             centre-of-mass integral is not negligible (raise e_low)"
        ));
    }
    let nodes = below + above + 1;
    let x_max = base + (nodes - 1) as f64 * dx;
    let z0 = base - TAIL - Z_CUTOFF;
    let z1 = x_max / sqrt_alpha + TAIL + Z_CUTOFF;
    let count = ((z1 - z0) / dz).ceil() as usize + 1;
    Ok(Mesh { scale, base, top: below, nodes, z0, count })
}

/// Prefix and suffix log-sums of W D_s over node s's z range.
struct Span {
    i_lo: usize,
    prefix: Vec<f64>,
    suffix: Vec<f64>,
}

struct Assembler<'a> {
    awr: f64,
    a: f64,
    sqrt_alpha: f64,
    front: f64,
    base: f64,
    dx: f64,
    z0: f64,
    dz: f64,
    log_w: &'a [f64],
}

impl Assembler<'_> {
    fn span(&self, s: usize) -> Span {
        let xs = self.base + s as f64 * self.dx;
        let n_w = self.log_w.len();
        let i_lo = ((xs - TAIL - Z_CUTOFF - self.z0) / self.dz).floor().max(0.0) as usize;
        let i_hi = (((xs / self.sqrt_alpha + TAIL + Z_CUTOFF - self.z0) / self.dz).ceil() as usize).min(n_w - 1);
        let sqrt_a = self.a.sqrt();
        let log_d: Vec<f64> = (i_lo..=i_hi)
            .map(|i| {
                let z = self.z0 + self.dz * i as f64;
                let gap = (xs - self.a * z).abs();
                xs * xs / self.awr - z * z / (self.awr + 1.0) - gap * gap / self.a
                    + (self.front * erfcx(gap / sqrt_a)).ln()
                    + self.log_w[i]
            })
            .collect();
        let mut prefix = log_d.clone();
        for i in 1..prefix.len() {
            prefix[i] = log_add(prefix[i - 1], log_d[i]);
        }
        let mut suffix = log_d.clone();
        for i in (0..suffix.len() - 1).rev() {
            suffix[i] = log_add(suffix[i + 1], log_d[i]);
        }
        Span { i_lo, prefix, suffix }
    }
}

fn prefix(sp: &Span, j: i64) -> f64 {
    if j < sp.i_lo as i64 {
        return f64::NEG_INFINITY;
    }
    sp.prefix[((j - sp.i_lo as i64) as usize).min(sp.prefix.len() - 1)]
}

fn suffix(sp: &Span, j: i64) -> f64 {
    let idx = j - sp.i_lo as i64;
    if idx >= sp.suffix.len() as i64 {
        return f64::NEG_INFINITY;
    }
    sp.suffix[idx.max(0) as usize]
}

pub struct Kernel {
    pub scale: f64,
    pub base: f64,
    pub dx: f64,
    pub top: usize,
    pub nodes: usize,
    pub starts: Vec<i64>,
    pub rows: Vec<i64>,
    pub values: Vec<f64>,
    pub norm: Vec<f64>,
}

#[allow(clippy::too_many_arguments)]
pub fn free_gas_kernel(
    energies: &[f64],
    sigma0: &[f64],
    awr: f64,
    temperature: f64,
    e_low: f64,
    e_top: f64,
    dx: f64,
    dz: f64,
) -> Result<Kernel, String> {
    let ms = mesh(awr, temperature, e_low, e_top, dx, dz)?;
    let w = hat_weights(ms.z0, dz, ms.count, energies, sigma0, ms.scale);
    let log_w: Vec<f64> = w.iter().map(|&v| if v > 0.0 { v.ln() } else { f64::NEG_INFINITY }).collect();
    let a = awr / (awr + 1.0);
    let asm = Assembler {
        awr,
        a,
        sqrt_alpha: (awr - 1.0) / (awr + 1.0),
        front: (PI * a).sqrt() / 2.0,
        base: ms.base,
        dx,
        z0: ms.z0,
        dz,
        log_w: &log_w,
    };
    let (scale, base, nodes) = (ms.scale, ms.base, ms.nodes);
    let x_low = (scale * e_low).sqrt();
    let first = ((x_low - TAIL - base) / dx).floor().max(0.0) as usize;
    let incoming: Vec<usize> = (first..=ms.top).collect();
    let n_w = log_w.len() as i64;
    let threads = std::thread::available_parallelism().map_or(1, |n| n.get());
    let chunk = incoming.len().div_ceil(threads).max(1);
    // per incoming node: (m, lo, normalized row, norm)
    let rows_out: Vec<(usize, usize, Vec<f64>, f64)> = std::thread::scope(|s| {
        let handles: Vec<_> = incoming
            .chunks(chunk)
            .map(|block| {
                let asm = &asm;
                s.spawn(move || {
                    let mut cache: HashMap<usize, Span> = HashMap::new();
                    let mut out = Vec::with_capacity(block.len());
                    for &m in block {
                        let xm = base + m as f64 * dx;
                        let lo = ((asm.sqrt_alpha * xm - TAIL - base) / dx).ceil().max(0.0) as usize;
                        let hi = (((xm + TAIL - base) / dx).floor() as usize).min(nodes - 1);
                        let mut k = vec![0.0; hi - lo + 1];
                        for n in lo..=hi {
                            for s_ in [m, n] {
                                cache.entry(s_).or_insert_with(|| asm.span(s_));
                            }
                            let xn = base + n as f64 * dx;
                            let sm = &cache[&m];
                            let value = if n == m {
                                prefix(sm, n_w).exp()
                            } else {
                                let j = (((xm + xn) / (2.0 * a) - asm.z0) / dz).floor() as i64;
                                let balance = (xm * xm - xn * xn) / awr;
                                let sn = &cache[&n];
                                if n < m {
                                    prefix(sm, j).exp() + (balance + suffix(sn, j + 1)).exp()
                                } else {
                                    suffix(sm, j + 1).exp() + (balance + prefix(sn, j)).exp()
                                }
                            };
                            k[n - lo] = value;
                        }
                        let factor = (awr + 1.0).powi(2) / (2.0 * awr * awr * xm * xm / scale * PI.sqrt());
                        let mut total = 0.0;
                        for (idx, v) in k.iter_mut().enumerate() {
                            *v *= factor;
                            let xn = base + dx * (lo + idx) as f64;
                            total += *v * 2.0 * xn * dx / scale;
                        }
                        k.iter_mut().for_each(|v| *v /= total);
                        out.push((m, lo, k, total));
                        cache.retain(|&s_, _| s_ + 1 >= lo);
                    }
                    out
                })
            })
            .collect();
        handles.into_iter().flat_map(|h| h.join().expect("worker panicked")).collect()
    });
    let mut norm = vec![f64::NAN; nodes];
    let mut triples: Vec<(usize, usize, f64)> = Vec::new();
    for (m, lo, k, total) in rows_out {
        norm[m] = total;
        for (idx, v) in k.into_iter().enumerate() {
            triples.push((lo + idx, m, v));
        }
    }
    triples.sort_by(|p, q| (p.0, p.1).cmp(&(q.0, q.1)));
    let mut starts = vec![0i64; nodes + 1];
    for &(n, _, _) in &triples {
        starts[n + 1] += 1;
    }
    for i in 1..=nodes {
        starts[i] += starts[i - 1];
    }
    let rows = triples.iter().map(|t| t.1 as i64).collect();
    let values = triples.iter().map(|t| t.2).collect();
    Ok(Kernel { scale, base, dx, top: ms.top, nodes, starts, rows, values, norm })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn erfcx_limits() {
        assert!((erfcx(0.0) - 1.0).abs() < 1e-15);
        // asymptotic 1 / (sqrt(pi) y) (1 - 1/(2 y^2) + 3/(4 y^4))
        let y = 20.0f64;
        let asym = 1.0 / (PI.sqrt() * y) * (1.0 - 0.5 / (y * y) + 0.75 / y.powi(4));
        assert!((erfcx(y) / asym - 1.0).abs() < 1e-7);
    }
}
