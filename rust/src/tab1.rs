//! ENDF tabulated functions. Rust port of `nodal/endf/interpolation.py`.

/// A function y(x) on ENDF points and interpolation regions.
/// Extracted from a Python `Tabulated1D` by attribute name.
#[derive(Clone, Debug, pyo3::FromPyObject)]
pub struct Tab1 {
    pub x: Vec<f64>,
    pub y: Vec<f64>,
    pub breaks: Vec<usize>, // NBT: 1-based index of the last point of each region
    pub laws: Vec<u8>,      // INT: interpolation law of each region
}

/// Value at x between (x1, y1) and (x2, y2) using one ENDF law. Mirrors `interpolate`.
pub fn interpolate(law: u8, x: f64, x1: f64, x2: f64, y1: f64, y2: f64) -> f64 {
    if law == 1 || x2 == x1 {
        return y1;
    }
    assert!(
        (2..=5).contains(&law),
        "unsupported interpolation law {law}"
    );
    let log_x = (law == 3 || law == 5) && x1 > 0.0 && x > 0.0;
    let log_y = (law == 4 || law == 5) && y1 > 0.0 && y2 > 0.0;
    let t = if log_x {
        (x / x1).ln() / (x2 / x1).ln()
    } else {
        (x - x1) / (x2 - x1)
    };
    if log_y {
        y1 * (t * (y2 / y1).ln()).exp()
    } else {
        y1 + t * (y2 - y1)
    }
}

impl Tab1 {
    /// Law of the interval between points i and i+1 (0-based).
    fn law(&self, i: usize) -> u8 {
        let r = self.breaks.partition_point(|&b| b < i + 2); // bisect_left
        self.laws[r]
    }

    /// Mirrors `Tabulated1D.__call__`: 0 outside the table, upper side at jumps.
    pub fn eval(&self, x: f64) -> f64 {
        let xs = &self.x;
        let n = xs.len();
        if x < xs[0] || x > xs[n - 1] {
            return 0.0;
        }
        if x == xs[n - 1] {
            return self.y[n - 1];
        }
        let i = xs.partition_point(|&v| v <= x) - 1; // bisect_right - 1
        interpolate(self.law(i), x, xs[i], xs[i + 1], self.y[i], self.y[i + 1])
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn log_log_is_exact_for_one_over_v() {
        let f = |e: f64| 3.0 / e.sqrt();
        let v = interpolate(5, 2.0, 1.0, 4.0, f(1.0), f(4.0));
        assert!((v / f(2.0) - 1.0).abs() < 1e-14);
    }

    #[test]
    fn regions_and_outside() {
        let t = Tab1 {
            x: vec![1., 2., 3., 4.],
            y: vec![10., 20., 30., 40.],
            breaks: vec![2, 4],
            laws: vec![1, 2],
        };
        assert_eq!(t.eval(1.5), 10.0);
        assert!((t.eval(2.5) - 25.0).abs() < 1e-12);
        assert_eq!(t.eval(0.5), 0.0);
    }
}
