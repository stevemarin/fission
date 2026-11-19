//! Ray tracing of tracks through a lattice of ring/sector cells (port of
//! `fission.transport.trace.trace_track`; same operations in the same
//! order, so the segments agree with Python to rounding).

use std::f64::consts::PI;

const TWO_PI: f64 = 2.0 * PI;
const MIN_LENGTH: f64 = 1e-10;

/// A lattice of square cells, bottom-left at the origin, cells numbered
/// iy * nx + ix. Each cell has ring radii, a sector count and the index of
/// its first flat-source region.
pub struct Lattice {
    pub pitch: f64,
    pub nx: usize,
    pub ny: usize,
    pub radii: Vec<Vec<f64>>,
    pub sectors: Vec<usize>,
    pub offsets: Vec<usize>,
}

impl Lattice {
    fn cell_at(&self, x: f64, y: f64) -> (usize, usize) {
        let clamp = |v: f64, n: usize| -> usize {
            let i = (v / self.pitch) as i64; // truncation, as int() in Python
            i.clamp(0, n as i64 - 1) as usize
        };
        (clamp(x, self.nx), clamp(y, self.ny))
    }

    /// Local region of the point (dx, dy) from the centre of a cell.
    fn local_region(&self, cell: usize, dx: f64, dy: f64) -> usize {
        let r2 = dx * dx + dy * dy;
        let ring = self.radii[cell].iter().filter(|&&r| r * r <= r2).count();
        let ns = self.sectors[cell];
        if ns == 1 {
            return ring;
        }
        let mut theta = dy.atan2(dx);
        if theta < 0.0 {
            theta += TWO_PI;
        }
        let sector = ((theta / (TWO_PI / ns as f64)) as usize).min(ns - 1);
        ring * ns + sector
    }

    /// (region, length) of every segment from (x0, y0) to (x1, y1).
    pub fn trace(&self, x0: f64, y0: f64, x1: f64, y1: f64) -> Vec<(usize, f64)> {
        let p = self.pitch;
        let length = ((x1 - x0) * (x1 - x0) + (y1 - y0) * (y1 - y0)).sqrt();
        let (dx, dy) = ((x1 - x0) / length, (y1 - y0) / length);
        let mut cuts = vec![0.0, length];
        if dx != 0.0 {
            for i in 1..self.nx {
                let t = (i as f64 * p - x0) / dx;
                if 0.0 < t && t < length {
                    cuts.push(t);
                }
            }
        }
        if dy != 0.0 {
            for j in 1..self.ny {
                let t = (j as f64 * p - y0) / dy;
                if 0.0 < t && t < length {
                    cuts.push(t);
                }
            }
        }
        cuts.sort_by(f64::total_cmp);
        let mut out: Vec<(usize, f64)> = Vec::new();
        for w in cuts.windows(2) {
            let (a, b) = (w[0], w[1]);
            if b - a < MIN_LENGTH {
                continue;
            }
            let tm = 0.5 * (a + b);
            let (ix, iy) = self.cell_at(x0 + tm * dx, y0 + tm * dy);
            let cell = iy * self.nx + ix;
            let (cx, cy) = ((ix as f64 + 0.5) * p, (iy as f64 + 0.5) * p);
            let mut local = vec![a, b];
            let (ox, oy) = (x0 - cx, y0 - cy);
            let bb = dx * ox + dy * oy;
            for &r in &self.radii[cell] {
                let disc = bb * bb - (ox * ox + oy * oy - r * r);
                if disc > 0.0 {
                    let root = disc.sqrt();
                    for t in [-bb - root, -bb + root] {
                        if a < t && t < b {
                            local.push(t);
                        }
                    }
                }
            }
            let ns = self.sectors[cell];
            if ns > 1 {
                for k in 0..ns {
                    let theta = TWO_PI * k as f64 / ns as f64;
                    let (ux, uy) = (theta.cos(), theta.sin());
                    let den = dx * uy - dy * ux;
                    if den != 0.0 {
                        let t = ((cx - x0) * uy - (cy - y0) * ux) / den;
                        if a < t && t < b {
                            local.push(t);
                        }
                    }
                }
            }
            local.sort_by(f64::total_cmp);
            for v in local.windows(2) {
                let (s, e) = (v[0], v[1]);
                if e - s < MIN_LENGTH {
                    continue;
                }
                let tm = 0.5 * (s + e);
                let region = self.offsets[cell]
                    + self.local_region(cell, x0 + tm * dx - cx, y0 + tm * dy - cy);
                match out.last_mut() {
                    Some(last) if last.0 == region => last.1 += e - s,
                    _ => out.push((region, e - s)),
                }
            }
        }
        out
    }
}

/// Trace many tracks on all cores; returns flat (regions, lengths) and the
/// start index of every track's segments (len = tracks + 1).
pub fn trace_all(
    lattice: &Lattice,
    x0: &[f64],
    y0: &[f64],
    x1: &[f64],
    y1: &[f64],
) -> (Vec<i64>, Vec<f64>, Vec<i64>) {
    let n = x0.len();
    let threads = std::thread::available_parallelism().map_or(1, |t| t.get());
    let chunk = n.div_ceil(threads).max(1);
    let per_track: Vec<Vec<(usize, f64)>> = std::thread::scope(|s| {
        let handles: Vec<_> = (0..n)
            .step_by(chunk)
            .map(|lo| {
                let hi = (lo + chunk).min(n);
                s.spawn(move || {
                    (lo..hi)
                        .map(|t| lattice.trace(x0[t], y0[t], x1[t], y1[t]))
                        .collect::<Vec<_>>()
                })
            })
            .collect();
        handles
            .into_iter()
            .flat_map(|h| h.join().expect("worker panicked"))
            .collect()
    });
    let total: usize = per_track.iter().map(Vec::len).sum();
    let mut regions = Vec::with_capacity(total);
    let mut lengths = Vec::with_capacity(total);
    let mut start = Vec::with_capacity(n + 1);
    start.push(0);
    for segs in per_track {
        for (r, l) in segs {
            regions.push(r as i64);
            lengths.push(l);
        }
        start.push(regions.len() as i64);
    }
    (regions, lengths, start)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn pin() -> Lattice {
        Lattice {
            pitch: 1.26,
            nx: 1,
            ny: 1,
            radii: vec![vec![0.4096, 0.418, 0.475]],
            sectors: vec![8],
            offsets: vec![0],
        }
    }

    #[test]
    fn chords_of_an_off_centre_track() {
        // horizontal, 0.1 cm above the centre: each circle gives a chord
        // 2 sqrt(r^2 - 0.01); the sector lines at 45, 90 and 135 degrees cut
        // the fuel chord into 4 pieces, so 7 + 3 segments
        let y = 0.63 + 0.1;
        let segs = pin().trace(0.0, y, 1.26, y);
        assert_eq!(segs.len(), 10);
        let total: f64 = segs.iter().map(|s| s.1).sum();
        assert!((total - 1.26).abs() < 1e-14);
        let fuel: f64 = segs.iter().filter(|s| s.0 / 8 == 0).map(|s| s.1).sum();
        assert!((fuel - 2.0 * (0.4096f64 * 0.4096 - 0.01).sqrt()).abs() < 1e-14);
    }

    #[test]
    fn a_horizontal_track_crosses_cells_and_merges_nothing_across_them() {
        let mut lat = pin();
        lat.nx = 2;
        lat.radii.push(vec![0.3]);
        lat.sectors.push(1);
        lat.offsets.push(32);
        let segs = lat.trace(0.0, 0.2, 2.52, 0.2);
        let total: f64 = segs.iter().map(|s| s.1).sum();
        assert!((total - 2.52).abs() < 1e-14);
        assert!(segs.iter().any(|s| s.0 >= 32) && segs.iter().any(|s| s.0 < 32));
        for w in segs.windows(2) {
            assert_ne!(w[0].0, w[1].0);
        }
    }
}
