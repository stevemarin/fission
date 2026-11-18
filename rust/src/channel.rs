//! Channel quantities. Rust port of `nodal/endf/channel.py`.

/// m_n c^2 (eV), CODATA 2018.
pub const NEUTRON_MASS_EV: f64 = 939.56542052e6;
/// hbar c (eV cm).
pub const HBAR_C_EV_CM: f64 = 197.3269804e6 * 1e-13;

/// k = K_CONSTANT * sqrt(E [eV]) in 1/(10^-12 cm).
pub fn k_constant() -> f64 {
    (2.0 * NEUTRON_MASS_EV).sqrt() / HBAR_C_EV_CM * 1e-12
}

/// Centre-of-mass neutron wave number k in 1/(10^-12 cm). Mirrors `wave_number`.
pub fn wave_number(energy: f64, awri: f64) -> f64 {
    k_constant() * awri / (awri + 1.0) * energy.abs().sqrt()
}

/// (P_l(rho), S_l(rho), phi_l(rho_hat)) for a hard sphere, l = 0..3.
/// Mirrors `penetrability_shift_phase`, including the order of operations.
pub fn penetrability_shift_phase(l: usize, rho: f64, rho_hat: f64) -> (f64, f64, f64) {
    let r2 = rho * rho;
    let h2 = rho_hat * rho_hat;
    match l {
        0 => (rho, 0.0, rho_hat),
        1 => (
            r2 * rho / (1.0 + r2),
            -1.0 / (1.0 + r2),
            rho_hat - rho_hat.atan(),
        ),
        2 => {
            let d = 9.0 + 3.0 * r2 + r2 * r2;
            (
                r2 * r2 * rho / d,
                -(18.0 + 3.0 * r2) / d,
                rho_hat - (3.0 * rho_hat).atan2(3.0 - h2),
            )
        }
        3 => {
            let d = 225.0 + 45.0 * r2 + 6.0 * r2 * r2 + r2 * r2 * r2;
            (
                r2 * r2 * r2 * rho / d,
                -(675.0 + 90.0 * r2 + 6.0 * r2 * r2) / d,
                rho_hat - (rho_hat * (15.0 - h2)).atan2(15.0 - 6.0 * h2),
            )
        }
        _ => panic!("l = {l} not implemented"),
    }
}

/// g_J = (2J+1) / (2 (2I+1)). Mirrors `spin_factor`.
pub fn spin_factor(j: f64, spi: f64) -> f64 {
    (2.0 * j.abs() + 1.0) / (2.0 * (2.0 * spi + 1.0))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn s_wave_penetrability_is_rho() {
        assert_eq!(penetrability_shift_phase(0, 0.3, 0.4), (0.3, 0.0, 0.4));
    }

    #[test]
    fn p_wave_suppressed_at_low_rho() {
        let (p, _, _) = penetrability_shift_phase(1, 1e-3, 1e-3);
        assert!((p / 1e-9 - 1.0).abs() < 1e-5);
    }
}
