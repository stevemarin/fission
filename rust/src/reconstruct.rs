//! Adaptive pointwise reconstruction. Rust port of `nodal/endf/reconstruct.py`.

/// Midpoint test for all three reactions. Mirrors `_converged`.
fn converged(left: &[f64; 3], mid: &[f64; 3], right: &[f64; 3], tol: f64, atol: f64) -> bool {
    (0..3).all(|j| (mid[j] - 0.5 * (left[j] + right[j])).abs() <= tol * mid[j].abs() + atol)
}

/// Refine `seeds` until linear interpolation of `f` is within tolerance.
/// Mirrors `adaptive_grid`: depth-first bisection with a stack of right ends.
pub fn adaptive_grid<F: Fn(f64) -> [f64; 3]>(
    f: F,
    seeds: &[f64],
    tol: f64,
    atol: f64,
    min_rel_width: f64,
) -> (Vec<f64>, Vec<[f64; 3]>) {
    let mut energies = vec![seeds[0]];
    let mut values = vec![f(seeds[0])];
    for &right_seed in &seeds[1..] {
        let mut stack = vec![(right_seed, f(right_seed))];
        while let Some(&(e2, v2)) = stack.last() {
            let e1 = *energies.last().unwrap();
            let v1 = *values.last().unwrap();
            if e2 - e1 > min_rel_width * e2 {
                let mid = 0.5 * (e1 + e2);
                let vm = f(mid);
                if !converged(&v1, &vm, &v2, tol, atol) {
                    stack.push((mid, vm));
                    continue;
                }
            }
            energies.push(e2);
            values.push(v2);
            stack.pop();
        }
    }
    (energies, values)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn linear_function_keeps_only_seeds() {
        let seeds = [1.0, 2.0, 7.0];
        let (e, _) = adaptive_grid(|x| [2.0 * x + 1.0, 3.0, 0.0], &seeds, 1e-3, 1e-5, 1e-9);
        assert_eq!(e, seeds);
    }
}
