//! Gauss-Legendre rules of any order, and tables of Legendre polynomials.

/// Nodes (ascending) and weights of the n-point Gauss-Legendre rule on
/// [-1, 1], by Newton's method on P_n (as numpy.polynomial.legendre.leggauss).
pub fn gauss_legendre(n: usize) -> (Vec<f64>, Vec<f64>) {
    let mut x = vec![0.0; n];
    let mut w = vec![0.0; n];
    for i in 0..n.div_ceil(2) {
        // Tricomi's initial guess for the i-th largest root
        let mut z = (std::f64::consts::PI * (i as f64 + 0.75) / (n as f64 + 0.5)).cos();
        let mut dp = 0.0;
        for _ in 0..100 {
            let (mut p0, mut p1) = (1.0, z);
            for k in 2..=n {
                let p2 = ((2 * k - 1) as f64 * z * p1 - (k - 1) as f64 * p0) / k as f64;
                p0 = p1;
                p1 = p2;
            }
            // P_n = p1, P_(n-1) = p0; derivative from the standard identity
            dp = n as f64 * (z * p1 - p0) / (z * z - 1.0);
            let dz = p1 / dp;
            z -= dz;
            if dz.abs() < 1e-16 {
                break;
            }
        }
        let weight = 2.0 / ((1.0 - z * z) * dp * dp);
        x[i] = -z;
        x[n - 1 - i] = z;
        w[i] = weight;
        w[n - 1 - i] = weight;
    }
    if n % 2 == 1 {
        x[n / 2] = 0.0;
    }
    (x, w)
}

/// P_0..P_lmax at mu (Bonnet's recurrence).
pub fn legendre(lmax: usize, mu: f64) -> Vec<f64> {
    let mut p = vec![1.0; lmax + 1];
    if lmax >= 1 {
        p[1] = mu;
    }
    for l in 1..lmax {
        p[l + 1] = ((2 * l + 1) as f64 * mu * p[l] - l as f64 * p[l - 1]) / (l + 1) as f64;
    }
    p
}

/// sum_l (2l+1)/2 a_l P_l(mu) by Clenshaw's recurrence (as numpy legval).
pub fn legendre_series(a: &[f64], mu: f64) -> f64 {
    // coefficients c_l = (2l+1)/2 a_l of the plain Legendre series
    let n = a.len();
    if n == 0 {
        return 0.0;
    }
    let c = |l: usize| (2 * l + 1) as f64 / 2.0 * a[l];
    if n == 1 {
        return c(0);
    }
    // numpy legval: c0 = c[-2], c1 = c[-1], then downward
    let (mut c0, mut c1) = (c(n - 2), c(n - 1));
    let mut nd = n as f64;
    for i in 3..=n {
        let tmp = c0;
        nd -= 1.0;
        c0 = c(n - i) - c1 * (nd - 1.0) / nd;
        c1 = tmp + c1 * mu * (2.0 * nd - 1.0) / nd;
    }
    c0 + c1 * mu
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rules_integrate_polynomials_exactly() {
        for n in [1, 2, 5, 16, 33] {
            let (x, w) = gauss_legendre(n);
            for degree in 0..(2 * n) {
                let integral: f64 = x
                    .iter()
                    .zip(&w)
                    .map(|(xi, wi)| wi * xi.powi(degree as i32))
                    .sum();
                let exact = if degree % 2 == 1 {
                    0.0
                } else {
                    2.0 / (degree + 1) as f64
                };
                assert!((integral - exact).abs() < 1e-13, "n={n} degree={degree}");
            }
        }
    }

    #[test]
    fn series_matches_direct_sum() {
        let a = [1.0, 0.3, -0.2, 0.05];
        let p = legendre(3, 0.4);
        let direct: f64 = (0..4).map(|l| (2 * l + 1) as f64 / 2.0 * a[l] * p[l]).sum();
        assert!((legendre_series(&a, 0.4) - direct).abs() < 1e-15);
    }
}
