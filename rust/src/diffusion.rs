//! Red-black SOR sweeps of the 7-point diffusion operator (port of
//! `fission.nodal.diffusion._sweep`, called `sweeps` times in a row for one
//! group): diag phi = b + sum over the six neighbours of coef * phi.
//!
//! Every node of one colour ((x + y + z) even or odd) depends only on nodes of
//! the other colour, so a half-sweep can update its nodes in any order, or in
//! parallel: the result is the same, bit for bit, and equal to the Python
//! reference, which evaluates the same expression for the whole array:
//!
//! ```text
//! s = b + (cxl phi_x- + cxu phi_x+) + (cyl phi_y- + cyu phi_y+)
//!       + (czl phi_z- + czu phi_z+),
//! phi += omega (s / diag - phi),
//! ```
//!
//! with a missing neighbour (outside the grid) counted as 0.

pub struct Grid {
    pub nz: usize,
    pub ny: usize,
    pub nx: usize,
}

pub struct Operator<'a> {
    pub diag: &'a [f64],
    /// (lower, upper) neighbour coefficients along x, y, z
    pub coef: [(&'a [f64], &'a [f64]); 3],
    pub inside: &'a [u8],
}

/// Nodes per thread below which a half-sweep stays on one thread (spawning
/// threads costs more than the work on a small mesh).
const PARALLEL_MIN: usize = 32_768;

#[allow(clippy::too_many_arguments)]
fn update(g: &Grid, op: &Operator, b: &[f64], phi: &[f64], c: usize, z: usize, y: usize,
          x: usize, omega: f64) -> f64 {
    let plane = g.ny * g.nx;
    let lo_x = if x > 0 { phi[c - 1] } else { 0.0 };
    let hi_x = if x + 1 < g.nx { phi[c + 1] } else { 0.0 };
    let lo_y = if y > 0 { phi[c - g.nx] } else { 0.0 };
    let hi_y = if y + 1 < g.ny { phi[c + g.nx] } else { 0.0 };
    let lo_z = if z > 0 { phi[c - plane] } else { 0.0 };
    let hi_z = if z + 1 < g.nz { phi[c + plane] } else { 0.0 };
    let [(xl, xu), (yl, yu), (zl, zu)] = op.coef;
    let mut s = b[c];
    s += xl[c] * lo_x + xu[c] * hi_x;
    s += yl[c] * lo_y + yu[c] * hi_y;
    s += zl[c] * lo_z + zu[c] * hi_z;
    let new = s / op.diag[c];
    phi[c] + omega * (new - phi[c])
}

/// New values of the colour's nodes in planes [z0, z1), written to `out`
/// (indexed from the start of plane z0); other entries are left alone.
#[allow(clippy::too_many_arguments)]
fn half(g: &Grid, op: &Operator, b: &[f64], phi: &[f64], colour: usize, z0: usize, z1: usize,
        omega: f64, out: &mut [f64]) {
    let plane = g.ny * g.nx;
    for z in z0..z1 {
        for y in 0..g.ny {
            let start = (z + y + colour) % 2; // first x of this colour on the row
            let row = z * plane + y * g.nx;
            for x in (start..g.nx).step_by(2) {
                let c = row + x;
                if op.inside[c] != 0 {
                    out[c - z0 * plane] = update(g, op, b, phi, c, z, y, x, omega);
                }
            }
        }
    }
}

pub fn sweeps(g: &Grid, op: &Operator, b: &[f64], phi: &mut [f64], count: usize, omega: f64) {
    let plane = g.ny * g.nx;
    let n = g.nz * plane;
    let threads = std::thread::available_parallelism()
        .map_or(1, |t| t.get())
        .min(g.nz)
        .min(n.div_ceil(PARALLEL_MIN))
        .max(1);
    let mut next = phi.to_vec();
    for _ in 0..count {
        for colour in 0..2 {
            if threads == 1 {
                half(g, op, b, phi, colour, 0, g.nz, omega, &mut next);
            } else {
                let per = g.nz.div_ceil(threads);
                let snapshot: &[f64] = phi;
                std::thread::scope(|s| {
                    for (t, chunk) in next.chunks_mut(per * plane).enumerate() {
                        let z0 = t * per;
                        let z1 = (z0 + per).min(g.nz);
                        s.spawn(move || half(g, op, b, snapshot, colour, z0, z1, omega, chunk));
                    }
                });
            }
            phi.copy_from_slice(&next);
        }
    }
}
