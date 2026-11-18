//! Reich-Moore cross sections. Rust port of `nodal/endf/reich_moore.py`.
//!
//! The physics and the order of operations follow the Python reference; the
//! one structural difference is that per-resonance quantities that do not
//! depend on the energy (penetrability at E_r, fission amplitudes) are
//! computed once in the constructor instead of on every call.

use crate::channel::{penetrability_shift_phase, spin_factor, wave_number};
use std::f64::consts::PI;
use std::ops::{Add, Div, Mul, Sub};

/// Minimal complex number (enough for a 3x3 R-matrix).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct C {
    pub re: f64,
    pub im: f64,
}

impl C {
    pub const ZERO: C = C { re: 0.0, im: 0.0 };
    pub const ONE: C = C { re: 1.0, im: 0.0 };
    pub fn new(re: f64, im: f64) -> C {
        C { re, im }
    }
    pub fn norm_sqr(self) -> f64 {
        self.re * self.re + self.im * self.im
    }
    pub fn abs(self) -> f64 {
        self.re.hypot(self.im)
    }
    pub fn scale(self, s: f64) -> C {
        C::new(self.re * s, self.im * s)
    }
}
impl Add for C {
    type Output = C;
    fn add(self, o: C) -> C {
        C::new(self.re + o.re, self.im + o.im)
    }
}
impl Sub for C {
    type Output = C;
    fn sub(self, o: C) -> C {
        C::new(self.re - o.re, self.im - o.im)
    }
}
impl Mul for C {
    type Output = C;
    fn mul(self, o: C) -> C {
        C::new(
            self.re * o.re - self.im * o.im,
            self.re * o.im + self.im * o.re,
        )
    }
}
impl Div for C {
    type Output = C;
    fn div(self, o: C) -> C {
        let d = o.norm_sqr();
        C::new(
            (self.re * o.re + self.im * o.im) / d,
            (self.im * o.re - self.re * o.im) / d,
        )
    }
}

// --------------------------------------------- inputs read from Python objects
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

// ------------------------------------------------------- prepared spin groups
struct Level {
    er: f64,
    gn: f64,       // neutron width at |E_r|
    gg: f64,       // gamma width
    pr: f64,       // penetrability at |E_r|
    fis: [f64; 2], // signed sqrt of the fission widths
}

struct SpinGroup {
    l: usize,
    g: f64,
    awri: f64,
    a: f64,
    ap: f64,
    nch: usize, // 1 (neutron) or 3 (neutron + two fission channels)
    levels: Vec<Level>,
}

pub struct ReichMoore {
    groups: Vec<SpinGroup>,
    potentials: Vec<PotentialTerm>,
}

/// One spin group at one energy (mirrors `SpinGroupSolution`).
pub struct Solution {
    pub u_nn: C,
    pub row: [C; 3],
    pub k: f64,
    pub gamma_flux: f64,
}

/// Solve a small complex system by Gaussian elimination with partial pivoting.
fn solve(mut m: [[C; 4]; 3], n: usize) -> [C; 3] {
    for col in 0..n {
        let mut pivot = col;
        for r in col + 1..n {
            if m[r][col].abs() > m[pivot][col].abs() {
                pivot = r;
            }
        }
        m.swap(col, pivot);
        for r in col + 1..n {
            let f = m[r][col] / m[col][col];
            for c in col..=n {
                m[r][c] = m[r][c] - f * m[col][c];
            }
        }
    }
    let mut x = [C::ZERO; 3];
    for r in (0..n).rev() {
        let mut s = C::ZERO;
        for c in r + 1..n {
            s = s + m[r][c] * x[c];
        }
        x[r] = (m[r][n] - s) / m[r][r];
    }
    x
}

impl ReichMoore {
    pub fn new(groups: Vec<SpinGroupIn>, potentials: Vec<PotentialTerm>, spi: f64) -> Self {
        let groups = groups
            .into_iter()
            .map(|sg| {
                let has_fission = sg.resonances.iter().any(|r| r.gfa != 0.0 || r.gfb != 0.0);
                let levels = sg
                    .resonances
                    .iter()
                    .map(|r| {
                        let kr = wave_number(r.energy, sg.awri);
                        let (pr, _, _) = penetrability_shift_phase(sg.l, kr * sg.a, kr * sg.ap);
                        let signed_sqrt = |f: f64| f.abs().sqrt().copysign(f);
                        Level {
                            er: r.energy,
                            gn: r.gn,
                            gg: r.gg,
                            pr,
                            fis: [signed_sqrt(r.gfa), signed_sqrt(r.gfb)],
                        }
                    })
                    .collect();
                SpinGroup {
                    l: sg.l,
                    g: spin_factor(sg.j, spi),
                    awri: sg.awri,
                    a: sg.a,
                    ap: sg.ap,
                    nch: if has_fission { 3 } else { 1 },
                    levels,
                }
            })
            .collect();
        ReichMoore { groups, potentials }
    }

    fn solve_group(sg: &SpinGroup, energy: f64) -> Solution {
        let k = wave_number(energy, sg.awri);
        let (p, _, phi) = penetrability_shift_phase(sg.l, k * sg.a, k * sg.ap);
        let n = sg.nch;
        let half_i = C::new(0.0, 0.5);
        let mut kmat = [[C::ZERO; 3]; 3];
        for lv in &sg.levels {
            let amp = [(lv.gn * p / lv.pr).sqrt(), lv.fis[0], lv.fis[1]];
            let scale = half_i / C::new(lv.er - energy, -0.5 * lv.gg);
            for c in 0..n {
                for c2 in 0..n {
                    kmat[c][c2] = kmat[c][c2] + scale.scale(amp[c]).scale(amp[c2]);
                }
            }
        }
        // augmented system (1 - K) x = e_n
        let mut m = [[C::ZERO; 4]; 3];
        for c in 0..n {
            for c2 in 0..n {
                let delta = if c == c2 { C::ONE } else { C::ZERO };
                m[c][c2] = delta - kmat[c][c2];
            }
            m[c][n] = if c == 0 { C::ONE } else { C::ZERO };
        }
        let row = solve(m, n);
        let phase = C::new((-2.0 * phi).cos(), (-2.0 * phi).sin());
        let u_nn = phase * (row[0].scale(2.0) - C::ONE);

        let mut gamma_flux = 0.0;
        for lv in &sg.levels {
            let amp = [(lv.gn * p / lv.pr).sqrt(), lv.fis[0], lv.fis[1]];
            let mut s = C::ZERO;
            for c in 0..n {
                s = s + row[c].scale(amp[c]);
            }
            let a_r = s / C::new(lv.er - energy, -0.5 * lv.gg);
            gamma_flux += lv.gg * a_r.abs().powi(2);
        }
        Solution {
            u_nn,
            row,
            k,
            gamma_flux,
        }
    }

    /// (elastic, capture, fission) in barns at one energy. Mirrors `reich_moore`.
    pub fn eval(&self, energy: f64) -> [f64; 3] {
        let (mut elastic, mut capture, mut fission) = (0.0, 0.0, 0.0);
        for sg in &self.groups {
            let sol = Self::solve_group(sg, energy);
            let pk = PI / (sol.k * sol.k);
            elastic += pk * sg.g * (C::ONE - sol.u_nn).abs().powi(2);
            if sg.nch == 3 {
                fission += 4.0 * pk * sg.g * (sol.row[1].abs().powi(2) + sol.row[2].abs().powi(2));
            }
            capture += pk * sg.g * sol.gamma_flux;
        }
        elastic += self.potentials.iter().map(|t| t.eval(energy)).sum::<f64>();
        [elastic, capture, fission]
    }
}
