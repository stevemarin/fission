//! Multilinear interpolation of a cross-section table at many points (port
//! of `fission.tables.multilinear`, called once per point).
//!
//! The table holds, at every node of a d-dimensional grid, a vector of m
//! numbers (all the constants of a node, side by side). At a point x, axis
//! after axis: find the segment [x_i, x_i+1] that holds x_a (the end
//! segment outside the grid, which extrapolates), and replace the leading
//! axis by
//!
//! ```text
//! v <- (1 - t) v[i] + t v[i + 1],   t = (x_a - x_i) / (x_i+1 - x_i),
//! ```
//!
//! an axis with a single node is constant (v <- v[0]). Each step is the
//! same elementwise expression as the numpy reference, evaluated in the
//! same order, so the result is equal bit for bit.

/// `values` is C-ordered [n_0, ..., n_d-1, m]; `points` is [p, d].
/// Returns [p, m].
pub fn multilinear(nodes: &[&[f64]], values: &[f64], m: usize, points: &[f64]) -> Vec<f64> {
    let d = nodes.len();
    let p = if d == 0 { 0 } else { points.len() / d };
    let mut out = Vec::with_capacity(p * m);
    // the size of one leading-axis slice before each step
    let mut strides = vec![m; d];
    for a in (0..d.saturating_sub(1)).rev() {
        strides[a] = strides[a + 1] * nodes[a + 1].len();
    }
    let mut buffer: Vec<f64> = Vec::new();
    let mut next: Vec<f64> = Vec::new();
    for x in points.chunks_exact(d) {
        let mut v: &[f64] = values;
        for a in 0..d {
            let (axis, stride) = (nodes[a], strides[a]);
            let n = axis.len();
            if n == 1 {
                next.clear();
                next.extend_from_slice(&v[..stride]);
            } else {
                let xa = x[a];
                // numpy searchsorted (side = left), minus one, clipped
                let i = axis.partition_point(|&node| node < xa).saturating_sub(1).min(n - 2);
                let t = (xa - axis[i]) / (axis[i + 1] - axis[i]);
                let (lo, hi) = (&v[i * stride..(i + 1) * stride], &v[(i + 1) * stride..(i + 2) * stride]);
                next.clear();
                next.extend(lo.iter().zip(hi).map(|(&l, &h)| (1.0 - t) * l + t * h));
            }
            std::mem::swap(&mut buffer, &mut next);
            v = &buffer;
        }
        out.extend_from_slice(v);
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bilinear_function_is_reproduced() {
        // f(x, y) = 1 + 2x + 3y + 4xy on x in {0, 1, 3}, y in {0, 2}
        let (xs, ys) = ([0.0, 1.0, 3.0], [0.0, 2.0]);
        let f = |x: f64, y: f64| 1.0 + 2.0 * x + 3.0 * y + 4.0 * x * y;
        let mut values = vec![];
        for &x in &xs {
            for &y in &ys {
                values.extend([f(x, y), -f(x, y)]);
            }
        }
        let points = [0.5, 1.0, 2.5, -1.0, 4.0, 3.0];
        let out = multilinear(&[&xs, &ys], &values, 2, &points);
        for (k, pt) in points.chunks(2).enumerate() {
            let want = f(pt[0], pt[1]);
            assert!((out[2 * k] - want).abs() < 1e-12);
            assert!((out[2 * k + 1] + want).abs() < 1e-12);
        }
    }

    #[test]
    fn single_node_axis_is_constant() {
        let out = multilinear(&[&[5.0], &[0.0, 1.0]], &[1.0, 3.0], 1, &[7.0, 0.25]);
        assert_eq!(out, vec![1.5]);
    }
}
