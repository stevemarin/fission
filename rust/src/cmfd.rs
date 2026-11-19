//! Coarse-mesh (CMFD) eigenvalue solve (port of
//! `fission.transport.cmfd.solve_coarse`): power iteration over k, groups
//! swept Gauss-Seidel in energy, and red-black SOR passes per group on the
//! 5-point operator diag phi_c + E phi_E + W phi_W + N phi_N + S phi_S.

/// The coarse operator; every array is (groups, ny, nx) except `inscatter`,
/// which is (g', g, ny, nx) = V Sigma_s,g'->g (zero on the diagonal).
pub struct Operator<'a> {
    pub groups: usize,
    pub ny: usize,
    pub nx: usize,
    pub diag: &'a [f64],
    pub east: &'a [f64],
    pub west: &'a [f64],
    pub north: &'a [f64],
    pub south: &'a [f64],
    pub inscatter: &'a [f64],
    pub fission: &'a [f64],
    pub chi: &'a [f64],
}

pub struct Settings {
    pub tol: f64,
    pub max_outer: usize,
    pub inner: usize,
    pub omega: f64,
}

/// Returns (k, phi, outer iterations); phi keeps the starting production.
pub fn solve(op: &Operator, phi0: &[f64], k0: f64, s: &Settings) -> (f64, Vec<f64>, usize) {
    let (gn, ny, nx) = (op.groups, op.ny, op.nx);
    let cells = ny * nx;
    let mut phi = phi0.to_vec();
    let total = |phi: &[f64]| -> f64 { op.fission.iter().zip(phi).map(|(f, p)| f * p).sum() };
    let target = total(&phi);
    let mut production = target;
    let mut k = k0;
    let mut fission = vec![0.0; cells];
    let mut b = vec![0.0; cells];
    for outer in 1..=s.max_outer {
        for (c, f) in fission.iter_mut().enumerate() {
            let mut sum = 0.0;
            for g in 0..gn {
                sum += op.fission[g * cells + c] * phi[g * cells + c];
            }
            *f = sum / k;
        }
        let old = phi.clone();
        for g in 0..gn {
            let row = g * cells;
            for c in 0..cells {
                let mut sum = 0.0;
                for h in 0..gn {
                    sum += op.inscatter[(h * gn + g) * cells + c] * phi[h * cells + c];
                }
                b[c] = op.chi[row + c] * fission[c] + sum;
            }
            for _ in 0..s.inner {
                for colour in 0..2 {
                    for y in 0..ny {
                        for x in 0..nx {
                            if (x + y) % 2 != colour {
                                continue;
                            }
                            let c = y * nx + x;
                            let at = |xx: usize, yy: usize| phi[row + yy * nx + xx];
                            let e = if x + 1 < nx { at(x + 1, y) } else { 0.0 };
                            let w = if x > 0 { at(x - 1, y) } else { 0.0 };
                            let n = if y + 1 < ny { at(x, y + 1) } else { 0.0 };
                            let so = if y > 0 { at(x, y - 1) } else { 0.0 };
                            let neighbours = op.east[row + c] * e
                                + op.west[row + c] * w
                                + op.north[row + c] * n
                                + op.south[row + c] * so;
                            let update = (b[c] - neighbours) / op.diag[row + c];
                            phi[row + c] = (1.0 - s.omega) * phi[row + c] + s.omega * update;
                        }
                    }
                }
            }
        }
        let new_production = total(&phi);
        let k_new = k * new_production / production;
        let scale = target / new_production;
        for p in phi.iter_mut() {
            *p *= scale;
        }
        production = target;
        let biggest = phi.iter().fold(0.0f64, |m, p| m.max(p.abs()));
        let change = phi
            .iter()
            .zip(&old)
            .fold(0.0f64, |m, (p, o)| m.max((p - o).abs()))
            / biggest;
        let converged = (k_new - k).abs() < s.tol && change < 100.0 * s.tol;
        k = k_new;
        if converged {
            return (k, phi, outer);
        }
    }
    (k, phi, s.max_outer)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn one_cell_one_group_gives_nu_sigma_f_over_removal() {
        // infinite medium: k = nu Sigma_f / Sigma_a
        let op = Operator {
            groups: 1,
            ny: 1,
            nx: 1,
            diag: &[0.1],
            east: &[0.0],
            west: &[0.0],
            north: &[0.0],
            south: &[0.0],
            inscatter: &[0.0],
            fission: &[0.12],
            chi: &[1.0],
        };
        let s = Settings {
            tol: 1e-13,
            max_outer: 100,
            inner: 2,
            omega: 1.0,
        };
        let (k, phi, _) = solve(&op, &[1.0], 1.0, &s);
        assert!((k - 1.2).abs() < 1e-12);
        assert!((phi[0] - 1.0).abs() < 1e-12);
    }
}
