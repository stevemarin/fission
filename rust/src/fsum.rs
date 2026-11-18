//! Correctly rounded floating-point sum: a port of CPython's `math.fsum`
//! (Shewchuk's algorithm with the final half-even correction). Because the
//! result is the exactly rounded sum, it is bit-identical to Python's.

pub fn fsum(values: impl IntoIterator<Item = f64>) -> f64 {
    // partials: non-overlapping, increasing magnitude, summing exactly to the input so far
    let mut partials: Vec<f64> = Vec::new();
    for mut x in values {
        let mut i = 0;
        for j in 0..partials.len() {
            let mut y = partials[j];
            if x.abs() < y.abs() {
                std::mem::swap(&mut x, &mut y);
            }
            let hi = x + y;
            let lo = y - (hi - x);
            if lo != 0.0 {
                partials[i] = lo;
                i += 1;
            }
            x = hi;
        }
        partials.truncate(i);
        if x != 0.0 {
            partials.push(x);
        }
    }

    // add the partials from the top until the sum is inexact
    let mut n = partials.len();
    if n == 0 {
        return 0.0;
    }
    n -= 1;
    let mut hi = partials[n];
    let mut lo = 0.0;
    while n > 0 {
        let x = hi;
        n -= 1;
        let y = partials[n];
        hi = x + y;
        lo = y - (hi - x);
        if lo != 0.0 {
            break;
        }
    }
    // round half to even correctly when the remaining partials push past a tie
    if n > 0 && ((lo < 0.0 && partials[n - 1] < 0.0) || (lo > 0.0 && partials[n - 1] > 0.0)) {
        let y = lo * 2.0;
        let x = hi + y;
        if y == x - hi {
            hi = x;
        }
    }
    hi
}

#[cfg(test)]
mod tests {
    use super::fsum;

    #[test]
    fn exact_where_naive_sum_fails() {
        assert_eq!(fsum([1e100, 1.0, -1e100, 1e-100]), 1.0);
        assert_eq!(fsum([0.1; 10]), 1.0);
        assert_eq!(fsum(std::iter::empty()), 0.0);
    }
}
