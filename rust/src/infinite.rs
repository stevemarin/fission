//! k-infinity of a homogeneous mixture by power iteration.
//! Port of `fission/pincell/infinite.py:solve` (which uses a dense eigen-solve).

/// LU factorization with partial pivoting of a dense n x n matrix (row-major).
pub struct Lu {
    n: usize,
    lu: Vec<f64>,
    pivot: Vec<usize>,
}

impl Lu {
    pub fn new(mut a: Vec<f64>, n: usize) -> Result<Self, String> {
        let mut pivot: Vec<usize> = (0..n).collect();
        for col in 0..n {
            let p = (col..n)
                .max_by(|&i, &j| a[i * n + col].abs().total_cmp(&a[j * n + col].abs()))
                .unwrap();
            if a[p * n + col] == 0.0 {
                return Err("singular matrix".into());
            }
            if p != col {
                for k in 0..n {
                    a.swap(p * n + k, col * n + k);
                }
                pivot.swap(p, col);
            }
            for row in col + 1..n {
                let f = a[row * n + col] / a[col * n + col];
                a[row * n + col] = f;
                for k in col + 1..n {
                    a[row * n + k] -= f * a[col * n + k];
                }
            }
        }
        Ok(Lu { n, lu: a, pivot })
    }

    pub fn solve(&self, b: &[f64]) -> Vec<f64> {
        let n = self.n;
        let mut x: Vec<f64> = self.pivot.iter().map(|&p| b[p]).collect();
        for i in 0..n {
            for k in 0..i {
                x[i] -= self.lu[i * n + k] * x[k];
            }
        }
        for i in (0..n).rev() {
            for k in i + 1..n {
                x[i] -= self.lu[i * n + k] * x[k];
            }
            x[i] /= self.lu[i * n + i];
        }
        x
    }
}

fn matvec(m: &[f64], v: &[f64], n: usize) -> Vec<f64> {
    (0..n)
        .map(|i| (0..n).map(|j| m[i * n + j] * v[j]).sum())
        .collect()
}

/// k and the flux for (diag(total) - S^T) phi = F phi / k.
/// `scatter[g * n + g']` is the transfer g -> g'; `fission[g * n + g']` = chi_g nuSigma_f,g'.
pub fn k_infinity(
    total: &[f64],
    scatter: &[f64],
    fission: &[f64],
    tol: f64,
) -> Result<(f64, Vec<f64>, usize), String> {
    let n = total.len();
    let mut a = vec![0.0; n * n];
    for g in 0..n {
        for gp in 0..n {
            a[g * n + gp] = -scatter[gp * n + g];
        }
        a[g * n + g] += total[g];
    }
    let lu = Lu::new(a, n)?;
    let mut phi = vec![1.0; n];
    let mut source = matvec(fission, &phi, n);
    let mut k = 1.0;
    for iteration in 1..=10_000 {
        let x = lu.solve(&source);
        let new_source = matvec(fission, &x, n);
        let k_new = new_source.iter().sum::<f64>() / source.iter().sum::<f64>();
        let converged = (k_new - k).abs() < tol * k_new;
        k = k_new;
        phi = x.iter().map(|v| v / k).collect();
        source = new_source.iter().map(|v| v / k).collect();
        if converged && iteration > 2 {
            return Ok((k, phi, iteration));
        }
    }
    Err("power iteration did not converge".into())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn two_group_formula() {
        let (st1, st2, s11, s12, s22, nf1, nf2) = (0.52, 1.30, 0.48, 0.02, 1.20, 0.008, 0.16);
        let total = [st1, st2];
        let scatter = [s11, s12, 0.0, s22];
        let fission = [nf1, nf2, 0.0, 0.0]; // chi = (1, 0)
        let (k, _, _) = k_infinity(&total, &scatter, &fission, 1e-14).unwrap();
        let expected = (nf1 + nf2 * s12 / (st2 - s22)) / (st1 - s11);
        assert!((k - expected).abs() < 1e-12 * expected, "{k} vs {expected}");
    }
}
