//! Compressed liquid water (port of `fission.feedback.water` and
//! `fission.feedback.convection`): IAPWS-IF97 region 1, the IAPWS 2008
//! viscosity and 2011 thermal conductivity in their industrial form, and
//! the Dittus-Boelter film coefficient with Presser's bundle factor.
//!
//! IF97 region 1 gives the dimensionless Gibbs free energy
//! gamma(pi, tau) = sum_i n_i (7.1 - pi)^I_i (tau - 1.222)^J_i; every
//! property follows from its derivatives. The sums here run term by term,
//! numpy sums pairwise, so the two agree to rounding (~1e-15 relative), not
//! bit for bit.

pub const R: f64 = 0.461526e3; // J/(kg K)
const P_STAR: f64 = 16.53e6; // Pa
const T_STAR: f64 = 1386.0; // K

/// IF97 Table 2: I_i, J_i, n_i
const TABLE2: [(f64, f64, f64); 34] = [
    (0.0, -2.0, 0.14632971213167),
    (0.0, -1.0, -0.84548187169114),
    (0.0, 0.0, -0.37563603672040e1),
    (0.0, 1.0, 0.33855169168385e1),
    (0.0, 2.0, -0.95791963387872),
    (0.0, 3.0, 0.15772038513228),
    (0.0, 4.0, -0.16616417199501e-1),
    (0.0, 5.0, 0.81214629983568e-3),
    (1.0, -9.0, 0.28319080123804e-3),
    (1.0, -7.0, -0.60706301565874e-3),
    (1.0, -1.0, -0.18990068218419e-1),
    (1.0, 0.0, -0.32529748770505e-1),
    (1.0, 1.0, -0.21841717175414e-1),
    (1.0, 3.0, -0.52838357969930e-4),
    (2.0, -3.0, -0.47184321073267e-3),
    (2.0, 0.0, -0.30001780793026e-3),
    (2.0, 1.0, 0.47661393906987e-4),
    (2.0, 3.0, -0.44141845330846e-5),
    (2.0, 17.0, -0.72694996297594e-15),
    (3.0, -4.0, -0.31679644845054e-4),
    (3.0, 0.0, -0.28270797985312e-5),
    (3.0, 6.0, -0.85205128120103e-9),
    (4.0, -5.0, -0.22425281908000e-5),
    (4.0, -2.0, -0.65171222895601e-6),
    (4.0, 10.0, -0.14341729937924e-12),
    (5.0, -8.0, -0.40516996860117e-6),
    (8.0, -11.0, -0.12734301741641e-8),
    (8.0, -6.0, -0.17424871230634e-9),
    (21.0, -29.0, -0.68762131295531e-18),
    (23.0, -31.0, 0.14478307828521e-19),
    (29.0, -38.0, 0.26335781662795e-22),
    (30.0, -39.0, -0.11947622640071e-22),
    (31.0, -40.0, 0.18228094581404e-23),
    (32.0, -41.0, -0.93537087292458e-25),
];

/// The derivatives of gamma at (p, t) used by the properties.
struct Gamma {
    pi: f64,
    tau: f64,
    tautau: f64,
    pipi: f64,
    pitau: f64,
}

fn gamma(p: f64, t: f64) -> Gamma {
    let (a, b) = (7.1 - p / P_STAR, T_STAR / t - 1.222);
    let mut g = Gamma { pi: 0.0, tau: 0.0, tautau: 0.0, pipi: 0.0, pitau: 0.0 };
    for &(i, j, n) in &TABLE2 {
        g.pi += -n * i * a.powf(i - 1.0) * b.powf(j);
        g.tau += n * a.powf(i) * j * b.powf(j - 1.0);
        g.tautau += n * a.powf(i) * j * (j - 1.0) * b.powf(j - 2.0);
        g.pipi += n * i * (i - 1.0) * a.powf(i - 2.0) * b.powf(j);
        g.pitau += -n * i * a.powf(i - 1.0) * j * b.powf(j - 1.0);
    }
    g
}

/// kg/m^3 at p (Pa), t (K).
pub fn density(p: f64, t: f64) -> f64 {
    1.0 / ((p / P_STAR) * gamma(p, t).pi * R * t / p)
}

/// J/kg.
pub fn enthalpy(p: f64, t: f64) -> f64 {
    (T_STAR / t) * gamma(p, t).tau * R * t
}

/// Isobaric heat capacity, J/(kg K).
pub fn heat_capacity(p: f64, t: f64) -> f64 {
    let tau = T_STAR / t;
    -(tau * tau) * gamma(p, t).tautau * R
}

/// The temperature (K) of enthalpy h (J/kg) at p (Pa): Newton on h(T),
/// slope c_p. None if it does not converge in 50 steps.
pub fn temperature(p: f64, h: f64, guess: f64) -> Option<f64> {
    let mut t = guess;
    for _ in 0..50 {
        let step = (enthalpy(p, t) - h) / heat_capacity(p, t);
        t -= step;
        if step.abs() < 1e-10 {
            return Some(t);
        }
    }
    None
}

// ---------------------------------------------------------------- transport
const T_REF: f64 = 647.096; // K
const RHO_REF: f64 = 322.0; // kg/m^3
const P_REF: f64 = 22.064e6; // Pa
const R_COND: f64 = 461.51805; // J/(kg K)

const H0: [f64; 4] = [1.67752, 2.20462, 0.6366564, -0.241605];
/// R12-08 Table 2: (i, j, H_ij), in the order of the Python dict
const H1: [(f64, f64, f64); 21] = [
    (0.0, 0.0, 5.20094e-1), (1.0, 0.0, 8.50895e-2), (2.0, 0.0, -1.08374),
    (3.0, 0.0, -2.89555e-1), (0.0, 1.0, 2.22531e-1), (1.0, 1.0, 9.99115e-1),
    (2.0, 1.0, 1.88797), (3.0, 1.0, 1.26613), (5.0, 1.0, 1.20573e-1),
    (0.0, 2.0, -2.81378e-1), (1.0, 2.0, -9.06851e-1), (2.0, 2.0, -7.72479e-1),
    (3.0, 2.0, -4.89837e-1), (4.0, 2.0, -2.57040e-1), (0.0, 3.0, 1.61913e-1),
    (1.0, 3.0, 2.57399e-1), (0.0, 4.0, -3.25372e-2), (3.0, 4.0, 6.98452e-2),
    (4.0, 5.0, 8.72102e-3), (3.0, 6.0, -4.35673e-3), (5.0, 6.0, -5.93264e-4),
];
const L0: [f64; 5] = [2.443221e-3, 1.323095e-2, 6.770357e-3, -3.454586e-3, 4.096266e-4];
const L1: [[f64; 6]; 5] = [
    [1.60397357, -0.646013523, 0.111443906, 0.102997357, -0.0504123634, 0.00609859258],
    [2.33771842, -2.78843778, 1.53616167, -0.463045512, 0.0832827019, -0.00719201245],
    [2.19650529, -4.54580785, 3.55777244, -1.40944978, 0.275418278, -0.0205938816],
    [-1.21051378, 1.60812989, -0.621178141, 0.0716373224, 0.0, 0.0],
    [-2.7203370, 4.57586331, -3.18369245, 1.1168348, -0.19268305, 0.012913842],
];
const LAMBDA: f64 = 177.8514;
const QD: f64 = 1.0 / 0.40;
const NU: f64 = 0.630;
const GAMMA_C: f64 = 1.239;
const XI0: f64 = 0.13;
const GAMMA0: f64 = 0.06;
const TR: f64 = 1.5;
const A: [[f64; 5]; 6] = [
    [6.53786807199516, 6.52717759281799, 5.35500529896124, 1.55225959906681, 1.11999926419994],
    [-5.61149954923348, -6.30816983387575, -3.96415689925446, 0.464621290821181, 0.595748562571649],
    [3.39624167361325, 8.08379285492595, 8.91990208918795, 8.93237374861479, 9.88952565078920],
    [-2.27492629730878, -9.82240510197603, -12.0338729505790, -11.0321960061126, -10.3255051147040],
    [10.2631854662709, 12.1358413791395, 9.19494865194302, 6.16780999933360, 4.66861294457414],
    [1.97815050331519, -5.54349664571295, -2.16866274479712, -0.965458722086812, -0.503243546373828],
];
const A_EDGES: [f64; 4] = [0.310559006, 0.776397516, 1.242236025, 1.863354037];

/// Pa s, R12-08 Eq. (10) with mu_2 = 1.
pub fn viscosity(t: f64, rho: f64) -> f64 {
    let (tb, rb) = (t / T_REF, rho / RHO_REF);
    let mut d = 0.0;
    for (i, h) in H0.iter().enumerate() {
        d += h / tb.powf(i as f64);
    }
    let mu0 = 100.0 * tb.sqrt() / d;
    let mut s = 0.0;
    for &(i, j, h) in &H1 {
        s += h * (1.0 / tb - 1.0).powf(i) * (rb - 1.0).powf(j);
    }
    1e-6 * mu0 * (rb * s).exp()
}

fn zeta_ref(rb: f64) -> f64 {
    let j = A_EDGES.iter().filter(|&&edge| rb > edge).count();
    let mut s = 0.0;
    for (i, row) in A.iter().enumerate() {
        s += row[j] * rb.powf(i as f64);
    }
    1.0 / s
}

/// W/(m K), R15-11 industrial formulation (critical enhancement included).
pub fn conductivity(p: f64, t: f64) -> f64 {
    let g = gamma(p, t);
    let tau = T_STAR / t;
    let rho = 1.0 / ((p / P_STAR) * g.pi * R * t / p);
    let (tb, rb) = (t / T_REF, rho / RHO_REF);
    let mut d = 0.0;
    for (k, c) in L0.iter().enumerate() {
        d += c / tb.powf(k as f64);
    }
    let lam0 = tb.sqrt() / d;
    let mut s = 0.0;
    for (i, row) in L1.iter().enumerate() {
        for (j, l) in row.iter().enumerate() {
            s += l * (1.0 / tb - 1.0).powf(i as f64) * (rb - 1.0).powf(j as f64);
        }
    }
    let lam1 = (rb * s).exp();
    let cp = -(tau * tau) * g.tautau * R;
    let cv = R * (-(tau * tau) * g.tautau + (g.pi - tau * g.pitau).powi(2) / g.pipi);
    let drho_dp = -(rho * rho) * g.pipi * R * t / (P_STAR * P_STAR);
    let zeta = drho_dp * P_REF / RHO_REF;
    let dchi = (rb * (zeta - zeta_ref(rb) * TR / tb)).max(0.0);
    let xi = XI0 * (dchi / GAMMA0).powf(NU / GAMMA_C);
    let y = QD * xi;
    let mut lam2 = 0.0;
    if y >= 1.2e-7 {
        let kappa = cp / cv;
        let z = 2.0 / (std::f64::consts::PI * y)
            * (((1.0 - 1.0 / kappa) * y.atan() + y / kappa)
                - (1.0 - (-1.0 / (1.0 / y + y * y / (3.0 * rb * rb))).exp()));
        let mub = viscosity(t, rho) / 1e-6;
        lam2 = LAMBDA * rb * (cp / R_COND) * tb / mub * z;
    }
    1e-3 * (lam0 * lam1 + lam2)
}

/// Presser's bundle factor psi for a square array.
pub fn presser(x: f64) -> f64 {
    0.9217 + 0.1478 * x - 0.1130 * (-7.0 * (x - 1.0)).exp()
}

/// Dittus-Boelter h (W/(m^2 K)) of water at p (Pa), t (K), mass flux G
/// (kg/(m^2 s)), hydraulic diameter D_h (m); Presser's factor if P/D given.
pub fn film_coefficient(p: f64, t: f64, mass_flux: f64, diameter: f64,
                        pitch_over_diameter: Option<f64>) -> f64 {
    let rho = density(p, t);
    let mu = viscosity(t, rho);
    let k = conductivity(p, t);
    let re = mass_flux * diameter / mu;
    let pr = heat_capacity(p, t) * mu / k;
    let mut nu = 0.023 * re.powf(0.8) * pr.powf(0.4);
    if let Some(x) = pitch_over_diameter {
        nu *= presser(x);
    }
    nu * k / diameter
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn if97_table5_values() {
        // IF97 Table 5: T = 300 K, p = 3 MPa and T = 500 K, p = 3 MPa
        let v = 1.0 / density(3e6, 300.0);
        assert!((v / 0.100215168e-2 - 1.0).abs() < 1e-8);
        assert!((enthalpy(3e6, 300.0) / 0.115331273e6 - 1.0).abs() < 1e-8);
        assert!((heat_capacity(3e6, 500.0) / 0.465580682e4 - 1.0).abs() < 1e-8);
        let t = temperature(15.5e6, enthalpy(15.5e6, 580.0), 560.0).unwrap();
        assert!((t - 580.0).abs() < 1e-9);
    }
}
