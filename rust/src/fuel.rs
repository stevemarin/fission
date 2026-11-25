//! Steady-state fuel rod temperatures (port of `fission.feedback.fuel.rod`):
//! film, clad (conductivity at the mean clad temperature, iterated), gap,
//! and the pellet by Newton on int_{T_s}^{T_c} k dT = q' / (4 pi), with the
//! OECD/NEA MOX benchmark (2003) properties. Same operations in the same
//! order as the Python reference, so equal bit for bit (up to libm's log).

use std::f64::consts::PI;

fn k_fuel(t: f64) -> f64 {
    1.05 + 2150.0 / (t - 73.15)
}

fn k_fuel_integral(t: f64) -> f64 {
    1.05 * t + 2150.0 * (t - 73.15).ln()
}

fn k_clad(t: f64) -> f64 {
    7.51 + 2.09e-2 * t - 1.45e-5 * t.powf(2.0) + 7.67e-9 * t.powf(3.0)
}

pub struct Rod {
    pub pellet_radius: f64,
    pub clad_inner: f64,
    pub clad_outer: f64,
    pub gap_conductance: f64,
}

/// (centre, pellet surface, clad inner, clad outer) temperatures (K) for a
/// linear heat rate q' (W/m), coolant temperature (K) and film coefficient
/// (W/(m^2 K), None = clad surface at the coolant temperature).
pub fn rod(r: &Rod, q: f64, coolant: f64, film: Option<f64>) -> [f64; 4] {
    let mut t_co = coolant;
    if let Some(h) = film {
        t_co += q / (2.0 * PI * r.clad_outer * h);
    }
    let mut t_ci = t_co;
    for _ in 0..50 {
        let new = t_co + q * (r.clad_outer / r.clad_inner).ln()
            / (2.0 * PI * k_clad(0.5 * (t_ci + t_co)));
        if (new - t_ci).abs() < 1e-9 {
            break;
        }
        t_ci = new;
    }
    let t_s = t_ci + q / (2.0 * PI * r.pellet_radius * r.gap_conductance);
    let target = k_fuel_integral(t_s) + q / (4.0 * PI);
    let mut t_c = t_s + q / (4.0 * PI * k_fuel(t_s));
    for _ in 0..50 {
        let step = (k_fuel_integral(t_c) - target) / k_fuel(t_c);
        t_c -= step;
        if step.abs() < 1e-9 {
            break;
        }
    }
    [t_c, t_s, t_ci, t_co]
}
