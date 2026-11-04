//! The boundary of a 2D signed distance field as a bag of oriented segments (marching squares on
//! a fine grid, every vertex projected exactly onto `φ = 0` by Newton steps along the exact
//! gradient). Used for area-uniform surface sampling of extrusions and revolutions and for the 2D
//! point-cloud generators (SPEC.md §11: "from 2D boundary samples, so no ray-casting bias").
//!
//! Accuracy: the polyline length and enclosed area are those of the true curve up to `O(h²κ²)`
//! relative on smooth parts and `O(h)` per sharp corner (`h` the cell size); sampled points are
//! projected back onto the curve, so they satisfy `|φ| < 1e-12`. Segments are oriented with the
//! inside on the left, so the shoelace sum is the (positive) area.

use rayon::prelude::*;

/// A traced boundary.
#[derive(Debug, Clone)]
pub struct Boundary2d {
    /// Oriented segments (inside on the left), vertices on `φ = 0`.
    pub segments: Vec<([f64; 2], [f64; 2])>,
    /// Total polyline length.
    pub length: f64,
    /// Enclosed area (shoelace sum of the segments).
    pub area: f64,
    /// The box that was traced.
    pub bbox: [[f64; 2]; 2],
    cum: Vec<f64>,
}

/// Project `p` onto `φ = 0` by `iters` Newton steps `p ← p − φ ∇φ` (exact for an exact SDF
/// away from the medial axis).
pub fn project(eval: &(dyn Fn([f64; 2]) -> (f64, [f64; 2]) + Sync), mut p: [f64; 2], iters: usize) -> [f64; 2] {
    for _ in 0..iters {
        let (d, g) = eval(p);
        if !d.is_finite() || d == 0.0 {
            break;
        }
        p = [p[0] - d * g[0], p[1] - d * g[1]];
    }
    p
}

impl Boundary2d {
    /// Trace `{φ = 0}` inside `bbox` (padded by 3 %) with `cells` cells along the longest side.
    pub fn trace(eval: &(dyn Fn([f64; 2]) -> (f64, [f64; 2]) + Sync), bbox: [[f64; 2]; 2], cells: usize) -> Boundary2d {
        let cells = cells.max(8);
        let ext = [bbox[1][0] - bbox[0][0], bbox[1][1] - bbox[0][1]];
        let pad = 0.03 * ext[0].max(ext[1]).max(1e-9);
        let lo = [bbox[0][0] - pad, bbox[0][1] - pad];
        let size = [ext[0] + 2.0 * pad, ext[1] + 2.0 * pad];
        let h = size[0].max(size[1]) / cells as f64;
        let n = [(size[0] / h).ceil() as usize, (size[1] / h).ceil() as usize];
        let (nx, ny) = (n[0] + 1, n[1] + 1);
        let at = |i: usize, j: usize| [lo[0] + h * i as f64, lo[1] + h * j as f64];
        let vals: Vec<f64> = (0..nx * ny).into_par_iter().map(|k| eval(at(k % nx, k / nx)).0).collect();
        let mut segments = Vec::new();
        for j in 0..n[1] {
            for i in 0..n[0] {
                let corners = [(i, j), (i + 1, j), (i + 1, j + 1), (i, j + 1)];
                let v: Vec<f64> = corners.iter().map(|&(a, b)| vals[b * nx + a]).collect();
                let ins: Vec<bool> = v.iter().map(|&x| x < 0.0).collect();
                if ins.iter().all(|&b| b) || ins.iter().all(|&b| !b) {
                    continue;
                }
                // Crossings in CCW edge order: (edge k, point, leaving?).
                let mut cr: Vec<(usize, [f64; 2], bool)> = Vec::new();
                for k in 0..4 {
                    let k1 = (k + 1) % 4;
                    if ins[k] != ins[k1] {
                        let (pa, pb) = (at(corners[k].0, corners[k].1), at(corners[k1].0, corners[k1].1));
                        let t = v[k] / (v[k] - v[k1]);
                        cr.push((k, [pa[0] + t * (pb[0] - pa[0]), pa[1] + t * (pb[1] - pa[1])], ins[k]));
                    }
                }
                let mut push = |a: &(usize, [f64; 2], bool), b: &(usize, [f64; 2], bool)| {
                    // From the leaving crossing to the entering one.
                    if a.2 { segments.push((a.1, b.1)) } else { segments.push((b.1, a.1)) }
                };
                if cr.len() == 2 {
                    push(&cr[0], &cr[1]);
                } else if cr.len() == 4 {
                    let centre = 0.25 * (v[0] + v[1] + v[2] + v[3]);
                    // Pair the crossings around the corners that the centre value separates.
                    let around_outside = centre < 0.0;
                    for k in 0..4 {
                        let corner = (k + 1) % 4; // corner between edge k and edge k+1
                        if ins[corner] != around_outside {
                            // edges k and k+1 both cross (alternating signs)
                            push(&cr[k], &cr[(k + 1) % 4]);
                        }
                    }
                }
            }
        }
        // Project the vertices onto the curve (shared vertices of neighbouring cells coincide
        // before projection, hence after it).
        let proj = |p: [f64; 2]| project(eval, p, 3);
        let segments: Vec<([f64; 2], [f64; 2])> =
            segments.into_par_iter().map(|(a, b)| (proj(a), proj(b))).filter(|(a, b)| a != b).collect();
        let mut cum = Vec::with_capacity(segments.len());
        let (mut length, mut area) = (0.0, 0.0);
        for (a, b) in &segments {
            length += (b[0] - a[0]).hypot(b[1] - a[1]);
            area += 0.5 * (a[0] * b[1] - b[0] * a[1]);
            cum.push(length);
        }
        Boundary2d { segments, length, area, bbox, cum }
    }

    /// The segment containing arclength `s ∈ [0, length)` and the parameter `t ∈ [0, 1]` in it.
    pub fn locate(&self, s: f64) -> Option<(usize, f64)> {
        if self.segments.is_empty() {
            return None;
        }
        let i = self.cum.partition_point(|&c| c <= s).min(self.segments.len() - 1);
        let prev = if i == 0 { 0.0 } else { self.cum[i - 1] };
        let len = self.cum[i] - prev;
        let t = if len > 0.0 { ((s - prev) / len).clamp(0.0, 1.0) } else { 0.5 };
        Some((i, t))
    }

    /// A point of the boundary, uniform in arclength from `u ∈ [0, 1)`, and its exact normal.
    pub fn sample(&self, eval: &(dyn Fn([f64; 2]) -> (f64, [f64; 2]) + Sync), u: f64) -> Option<([f64; 2], [f64; 2])> {
        let (i, t) = self.locate(u * self.length)?;
        let (a, b) = self.segments[i];
        let p = project(eval, [a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1])], 3);
        let (_, g) = eval(p);
        let l = g[0].hypot(g[1]);
        if l > 0.0 { Some((p, [g[0] / l, g[1] / l])) } else { None }
    }
}
