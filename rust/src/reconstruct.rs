//! Adaptive pointwise reconstruction. Rust port of `fission/endf/reconstruct.py`.

use std::collections::HashSet;

/// Intervals (as bit patterns of their end energies) that must not be bisected:
/// shaded pairs straddling a discontinuity. Mirrors `unsplittable`.
pub type Unsplittable = HashSet<(u64, u64)>;

/// Midpoint test for every reaction. Mirrors `_converged`.
fn converged(left: &[f64], mid: &[f64], right: &[f64], tol: f64, atol: f64) -> bool {
    left.iter()
        .zip(mid)
        .zip(right)
        .all(|((a, m), b)| (m - 0.5 * (a + b)).abs() <= tol * m.abs() + atol)
}

/// Refine `seeds` until linear interpolation of `f` is within tolerance.
/// Mirrors `adaptive_grid`: depth-first bisection with a stack of right ends.
/// `f` returns one value per reaction (a fixed array, or a Vec).
pub fn adaptive_grid<V: AsRef<[f64]>, F: Fn(f64) -> V>(
    f: F,
    seeds: &[f64],
    tol: f64,
    atol: f64,
    min_rel_width: f64,
    unsplittable: &Unsplittable,
) -> (Vec<f64>, Vec<V>) {
    let mut energies = vec![seeds[0]];
    let mut values = vec![f(seeds[0])];
    for &right_seed in &seeds[1..] {
        let mut stack = vec![(right_seed, f(right_seed))];
        while let Some((e2, v2)) = stack.last() {
            let e1 = *energies.last().unwrap();
            if *e2 - e1 > min_rel_width * *e2
                && !unsplittable.contains(&(e1.to_bits(), e2.to_bits()))
            {
                let mid = 0.5 * (e1 + *e2);
                let vm = f(mid);
                let v1 = values.last().unwrap();
                if !converged(v1.as_ref(), vm.as_ref(), v2.as_ref(), tol, atol) {
                    stack.push((mid, vm));
                    continue;
                }
            }
            let (e2, v2) = stack.pop().unwrap();
            energies.push(e2);
            values.push(v2);
        }
    }
    (energies, values)
}

/// A refined grid: energies and the values at each.
type Grid<V> = (Vec<f64>, Vec<V>);

/// `adaptive_grid` on all cores. Bisection between two seeds never looks
/// outside them, so chunks of seeds are refined independently and joined;
/// the result is identical to the serial one.
pub fn adaptive_grid_parallel<V, F>(
    f: F,
    seeds: &[f64],
    tol: f64,
    atol: f64,
    min_rel_width: f64,
    unsplittable: &Unsplittable,
) -> (Vec<f64>, Vec<V>)
where
    V: AsRef<[f64]> + Send,
    F: Fn(f64) -> V + Sync,
{
    use std::sync::Mutex;
    use std::sync::atomic::{AtomicUsize, Ordering};

    let threads = std::thread::available_parallelism().map_or(1, |n| n.get());
    let intervals = seeds.len().saturating_sub(1);
    if threads == 1 || intervals < 64 {
        return adaptive_grid(f, seeds, tol, atol, min_rel_width, unsplittable);
    }
    // many more chunks than threads, handed out on demand, for load balance
    let n_chunks = (threads * 16).min(intervals);
    let bounds: Vec<usize> = (0..=n_chunks).map(|c| c * intervals / n_chunks).collect();
    let next = AtomicUsize::new(0);
    let results: Mutex<Vec<Option<Grid<V>>>> = Mutex::new((0..n_chunks).map(|_| None).collect());
    std::thread::scope(|s| {
        for _ in 0..threads {
            s.spawn(|| {
                loop {
                    let c = next.fetch_add(1, Ordering::Relaxed);
                    if c >= n_chunks {
                        break;
                    }
                    let chunk = &seeds[bounds[c]..=bounds[c + 1]];
                    let r = adaptive_grid(&f, chunk, tol, atol, min_rel_width, unsplittable);
                    results.lock().unwrap()[c] = Some(r);
                }
            });
        }
    });
    let mut energies = Vec::new();
    let mut values = Vec::new();
    for (c, r) in results.into_inner().unwrap().into_iter().enumerate() {
        let (e, v) = r.expect("every chunk is processed");
        let skip = usize::from(c > 0); // first point repeats the previous chunk's last
        energies.extend(e.into_iter().skip(skip));
        values.extend(v.into_iter().skip(skip));
    }
    (energies, values)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn linear_function_keeps_only_seeds() {
        let seeds = [1.0, 2.0, 7.0];
        let (e, _) = adaptive_grid(
            |x| [2.0 * x + 1.0, 3.0, 0.0],
            &seeds,
            1e-3,
            1e-5,
            1e-9,
            &Unsplittable::new(),
        );
        assert_eq!(e, seeds);
    }

    #[test]
    fn parallel_grid_equals_serial() {
        let seeds: Vec<f64> = (0..500).map(|i| 1.0 + i as f64 * 0.01).collect();
        let f = |x: f64| [1.0 / (1.0 + 400.0 * (x - 3.3).powi(2)), x.sin()];
        let none = Unsplittable::new();
        let serial = adaptive_grid(f, &seeds, 1e-4, 1e-8, 1e-9, &none);
        let parallel = adaptive_grid_parallel(f, &seeds, 1e-4, 1e-8, 1e-9, &none);
        assert_eq!(serial.0, parallel.0);
        assert_eq!(serial.1, parallel.1);
    }

    #[test]
    fn unsplittable_interval_is_kept() {
        let seeds = [1.0, 2.0, 3.0];
        let step = |x: f64| [if x < 2.5 { 0.0 } else { 1.0 }];
        let keep: Unsplittable = [(2.0f64.to_bits(), 3.0f64.to_bits())].into();
        let (e, _) = adaptive_grid(step, &seeds, 1e-3, 1e-5, 1e-9, &keep);
        assert_eq!(e, seeds);
    }
}
