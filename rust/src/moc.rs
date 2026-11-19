//! One MOC transport sweep (port of `fission.transport.moc.sweep`).
//!
//! Every track is swept forward and backward. Within one sweep the tracks are
//! independent (their incoming fluxes come from the previous iteration), so
//! they are split over threads; each thread keeps its own region tally and
//! the tallies are added at the end.

/// Track data shared by all sweeps.
pub struct Tracks<'a> {
    pub region: &'a [i64],
    pub length: &'a [f64],
    pub start: &'a [i64],
    pub next_fwd: &'a [i64],
    pub fwd_forward: &'a [i64],
    pub next_bwd: &'a [i64],
    pub bwd_forward: &'a [i64],
    /// W[t, p] = 2 pi w_a omega_p delta_a sin(theta_p), flattened (tracks x polar)
    pub weights: &'a [f64],
    pub sin: &'a [f64],
    pub groups: usize,
}

struct Partial {
    tally: Vec<f64>,
    handoff: Vec<(usize, Vec<f64>)>, // (destination slot t * 2 + d, psi)
    leakage: f64,
}

/// Sweep all tracks: returns (tally of W delta psi per region and group,
/// incoming fluxes for the next iteration, leakage).
pub fn sweep(
    tr: &Tracks,
    sigma_t: &[f64],
    source: &[f64],
    psi_in: &[f64],
) -> (Vec<f64>, Vec<f64>, f64) {
    let g_n = tr.groups;
    let p_n = tr.sin.len();
    let n_tracks = tr.start.len() - 1;
    let ratio: Vec<f64> = source.iter().zip(sigma_t).map(|(q, s)| q / s).collect();
    let inv_sin: Vec<f64> = tr.sin.iter().map(|s| 1.0 / s).collect();
    let slot = p_n * g_n;

    let run = |lo: usize, hi: usize| -> Partial {
        let mut tally = vec![0.0; sigma_t.len()];
        let mut handoff = Vec::with_capacity(2 * (hi - lo));
        let mut leakage = 0.0;
        let mut delta = vec![0.0; slot];
        for t in lo..hi {
            let (a, b) = (tr.start[t] as usize, tr.start[t + 1] as usize);
            let w = &tr.weights[t * p_n..(t + 1) * p_n];
            for d in 0..2 {
                let mut psi = psi_in[(t * 2 + d) * slot..(t * 2 + d + 1) * slot].to_vec();
                for step in 0..(b - a) {
                    let k = if d == 0 { a + step } else { b - 1 - step };
                    let r = tr.region[k] as usize;
                    let row = r * g_n;
                    for p in 0..p_n {
                        let l3 = tr.length[k] * inv_sin[p];
                        for g in 0..g_n {
                            let tau = sigma_t[row + g] * l3;
                            delta[p * g_n + g] =
                                (psi[p * g_n + g] - ratio[row + g]) * -(-tau).exp_m1();
                        }
                    }
                    for g in 0..g_n {
                        let mut sum = 0.0;
                        for p in 0..p_n {
                            sum += w[p] * delta[p * g_n + g];
                        }
                        tally[row + g] += sum;
                    }
                    for (x, dx) in psi.iter_mut().zip(&delta) {
                        *x -= dx;
                    }
                }
                let (next, forward) = if d == 0 {
                    (tr.next_fwd[t], tr.fwd_forward[t])
                } else {
                    (tr.next_bwd[t], tr.bwd_forward[t])
                };
                if next < 0 {
                    for p in 0..p_n {
                        for g in 0..g_n {
                            leakage += w[p] * psi[p * g_n + g];
                        }
                    }
                } else {
                    let dest = next as usize * 2 + if forward != 0 { 0 } else { 1 };
                    handoff.push((dest, psi));
                }
            }
        }
        Partial {
            tally,
            handoff,
            leakage,
        }
    };

    let threads = std::thread::available_parallelism().map_or(1, |n| n.get());
    let chunk = n_tracks.div_ceil(threads).max(1);
    let parts: Vec<Partial> = std::thread::scope(|s| {
        let handles: Vec<_> = (0..n_tracks)
            .step_by(chunk)
            .map(|lo| {
                let run = &run;
                s.spawn(move || run(lo, (lo + chunk).min(n_tracks)))
            })
            .collect();
        handles
            .into_iter()
            .map(|h| h.join().expect("worker panicked"))
            .collect()
    });

    let mut tally = vec![0.0; sigma_t.len()];
    let mut psi_next = vec![0.0; psi_in.len()];
    let mut leakage = 0.0;
    for part in parts {
        for (x, y) in tally.iter_mut().zip(&part.tally) {
            *x += y;
        }
        for (dest, psi) in part.handoff {
            psi_next[dest * slot..(dest + 1) * slot].copy_from_slice(&psi);
        }
        leakage += part.leakage;
    }
    (tally, psi_next, leakage)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_pure_absorber_with_a_flat_source_reaches_q_over_sigma() {
        // one long segment, one group: psi_out -> Q / Sigma_t from psi_in = 0
        let tr = Tracks {
            region: &[0],
            length: &[100.0],
            start: &[0, 1],
            next_fwd: &[-1],
            fwd_forward: &[1],
            next_bwd: &[-1],
            bwd_forward: &[1],
            weights: &[1.0],
            sin: &[1.0],
            groups: 1,
        };
        let (tally, next, leakage) = sweep(&tr, &[1.0], &[2.0], &[0.0, 0.0]);
        // each sweep: delta psi = -2 (1 - e^-100); leaks psi_out = 2 (1 - e^-100)
        assert!((tally[0] + 4.0).abs() < 1e-12);
        assert!((leakage - 4.0).abs() < 1e-12);
        assert_eq!(next, vec![0.0, 0.0]);
    }
}
