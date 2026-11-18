//! Multi-Level Breit-Wigner cross sections. Rust port of `fission/endf/mlbw.py`
//! (together with the single-level pieces it uses from `slbw.single_resonance`).
//!
//! Capture and fission are sums of single-level terms; elastic scattering uses
//! the multi-level collision-matrix element
//!     U = exp(-2 i phi) [1 + i sum_r Gn_r / (E'_r - E - i G_r / 2)].
//! Energy-independent per-level quantities (penetrability and shift at E_r,
//! spin factor, fission width) are computed once in the constructor.

use crate::channel::{penetrability_shift_phase, spin_factor, wave_number};
use crate::complex::C;
use crate::spin_groups::{PotentialTerm, SpinGroupIn};
use std::f64::consts::PI;

struct Level {
    er: f64,
    gn: f64, // neutron width at |E_r|
    gg: f64, // gamma width
    gf: f64, // total fission width |GFA| + |GFB|
    pr: f64, // penetrability at |E_r|
    sr: f64, // shift factor at |E_r|
    g: f64,  // spin factor g_J
}

struct Group {
    l: usize,
    awri: f64,
    a: f64,
    ap: f64,
    levels: Vec<Level>,
}

pub struct Mlbw {
    groups: Vec<Group>,
    potentials: Vec<PotentialTerm>,
}

impl Mlbw {
    pub fn new(groups: Vec<SpinGroupIn>, potentials: Vec<PotentialTerm>, spi: f64) -> Self {
        let groups = groups
            .into_iter()
            .map(|sg| {
                let levels = sg
                    .resonances
                    .iter()
                    .map(|r| {
                        let kr = wave_number(r.energy, sg.awri);
                        let (pr, sr, _) = penetrability_shift_phase(sg.l, kr * sg.a, kr * sg.ap);
                        Level {
                            er: r.energy,
                            gn: r.gn,
                            gg: r.gg,
                            gf: r.gfa.abs() + r.gfb.abs(),
                            pr,
                            sr,
                            g: spin_factor(sg.j, spi),
                        }
                    })
                    .collect();
                Group {
                    l: sg.l,
                    awri: sg.awri,
                    a: sg.a,
                    ap: sg.ap,
                    levels,
                }
            })
            .collect();
        Mlbw { groups, potentials }
    }

    /// (elastic, capture, fission) in barns at one energy. Mirrors `mlbw`.
    pub fn eval(&self, energy: f64) -> [f64; 3] {
        let (mut elastic, mut capture, mut fission) = (0.0, 0.0, 0.0);
        for sg in &self.groups {
            let k = wave_number(energy, sg.awri);
            let (p, s, phi) = penetrability_shift_phase(sg.l, k * sg.a, k * sg.ap);
            let mut amplitude = C::ONE;
            let mut g = 0.0;
            for lv in &sg.levels {
                let gn = lv.gn * p / lv.pr;
                let shifted = lv.er + (lv.sr - s) * lv.gn / (2.0 * lv.pr);
                let width = gn + lv.gg + lv.gf;
                let de = energy - shifted;
                let denom = de * de + 0.25 * width * width;
                let scale = PI / (k * k) * lv.g / denom;
                capture += scale * gn * lv.gg;
                fission += scale * gn * lv.gf;
                amplitude = amplitude + C::new(0.0, gn) / C::new(shifted - energy, -0.5 * width);
                g = lv.g;
            }
            let phase = C::new((-2.0 * phi).cos(), (-2.0 * phi).sin());
            let u_nn = phase * amplitude;
            elastic += PI / (k * k) * g * (C::ONE - u_nn).abs().powi(2);
        }
        elastic += self.potentials.iter().map(|t| t.eval(energy)).sum::<f64>();
        [elastic, capture, fission]
    }
}
