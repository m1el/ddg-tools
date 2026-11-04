//! 2D generators for the trained 2D network (SPEC.md §18 item 4; WP9b decides the training
//! format): shape instances from the Quilez SDFs, random simple polygons, 2D noise blobs, a
//! simulated 2D scanner, oriented 2D point clouds with exact normals, and the ground-truth 2D
//! signed distance.
//!
//! Everything is deterministic in its seed. A [`Curve2d`] is a bag of oriented segments with the
//! inside on the left (a closed polygon, or the traced zero level set of a noise implicit); its
//! signed distance is exact for the polyline (distance to the nearest segment, even-odd sign) and
//! costs `O(segments)` per query.

use super::perlin::Perlin;
use super::rng::Rng;
use crate::PatError;
use crate::analytic::boundary::Boundary2d;
use crate::analytic::sdf2d::{ALL_KINDS, Shape2d};
use rayon::prelude::*;
use std::f64::consts::{PI, TAU};

/// An oriented 2D point cloud.
#[derive(Debug, Clone, PartialEq)]
pub struct Cloud2d {
    /// Positions.
    pub positions: Vec<[f64; 2]>,
    /// Unit normals.
    pub normals: Vec<[f64; 2]>,
}

/// A closed curve as oriented segments (inside on the left).
#[derive(Debug, Clone)]
pub struct Curve2d {
    segments: Vec<([f64; 2], [f64; 2])>,
    cum: Vec<f64>,
    length: f64,
    bbox: [[f64; 2]; 2],
}

impl Curve2d {
    /// From segments (non-empty; zero-length ones are dropped).
    pub fn new(segments: Vec<([f64; 2], [f64; 2])>) -> Result<Curve2d, PatError> {
        let segments: Vec<_> = segments.into_iter().filter(|(a, b)| a != b && a.iter().chain(b).all(|c| c.is_finite())).collect();
        if segments.is_empty() {
            return Err(PatError::Format { reason: "curve has no segments".into() });
        }
        let mut cum = Vec::with_capacity(segments.len());
        let mut length = 0.0;
        let mut lo = [f64::INFINITY; 2];
        let mut hi = [f64::NEG_INFINITY; 2];
        for (a, b) in &segments {
            length += (b[0] - a[0]).hypot(b[1] - a[1]);
            cum.push(length);
            for p in [a, b] {
                for k in 0..2 {
                    lo[k] = lo[k].min(p[k]);
                    hi[k] = hi[k].max(p[k]);
                }
            }
        }
        Ok(Curve2d { segments, cum, length, bbox: [lo, hi] })
    }

    /// A closed polygon (any orientation; stored counter-clockwise). At least 3 vertices.
    pub fn polygon(vertices: &[[f64; 2]]) -> Result<Curve2d, PatError> {
        if vertices.len() < 3 {
            return Err(PatError::TooFewPoints { n: vertices.len(), need: 3 });
        }
        let mut v = vertices.to_vec();
        let area: f64 = (0..v.len()).map(|i| v[i][0] * v[(i + 1) % v.len()][1] - v[(i + 1) % v.len()][0] * v[i][1]).sum();
        if area < 0.0 {
            v.reverse();
        }
        Curve2d::new((0..v.len()).map(|i| (v[i], v[(i + 1) % v.len()])).collect())
    }

    /// The segments.
    pub fn segments(&self) -> &[([f64; 2], [f64; 2])] {
        &self.segments
    }

    /// Total length.
    pub fn length(&self) -> f64 {
        self.length
    }

    /// Bounding box of the vertices.
    pub fn bbox(&self) -> [[f64; 2]; 2] {
        self.bbox
    }

    /// Signed distance (negative inside) and the unit direction `∇φ` (the segment normal on a
    /// segment interior, `s (p − v)/|p − v|` at a vertex).
    pub fn signed_distance(&self, p: [f64; 2]) -> (f64, [f64; 2]) {
        let mut best = (f64::INFINITY, [0.0; 2], 0usize, 0.0);
        let mut inside = false;
        for (i, &(a, b)) in self.segments.iter().enumerate() {
            let e = [b[0] - a[0], b[1] - a[1]];
            let w = [p[0] - a[0], p[1] - a[1]];
            let t = ((w[0] * e[0] + w[1] * e[1]) / (e[0] * e[0] + e[1] * e[1])).clamp(0.0, 1.0);
            let q = [a[0] + t * e[0], a[1] + t * e[1]];
            let d2 = (p[0] - q[0]).powi(2) + (p[1] - q[1]).powi(2);
            if d2 < best.0 {
                best = (d2, q, i, t);
            }
            if (a[1] > p[1]) != (b[1] > p[1]) && p[0] < a[0] + (p[1] - a[1]) / (b[1] - a[1]) * (b[0] - a[0]) {
                inside = !inside;
            }
        }
        let (d2, q, i, t) = best;
        let d = d2.sqrt();
        let s = if inside { -1.0 } else { 1.0 };
        let (a, b) = self.segments[i];
        let l = (b[0] - a[0]).hypot(b[1] - a[1]);
        let n = [(b[1] - a[1]) / l, -(b[0] - a[0]) / l]; // right of the direction = outward
        let g = if t > 1e-12 && t < 1.0 - 1e-12 || d < 1e-13 { n } else { [s * (p[0] - q[0]) / d, s * (p[1] - q[1]) / d] };
        (s * d, g)
    }

    /// `n` points uniform in arclength with the outward normals of their segments.
    pub fn sample_uniform(&self, n: usize, rng: &mut Rng) -> Cloud2d {
        let mut positions = Vec::with_capacity(n);
        let mut normals = Vec::with_capacity(n);
        for _ in 0..n {
            let s = rng.f64() * self.length;
            let i = self.cum.partition_point(|&c| c <= s).min(self.segments.len() - 1);
            let prev = if i == 0 { 0.0 } else { self.cum[i - 1] };
            let (a, b) = self.segments[i];
            let l = self.cum[i] - prev;
            let t = ((s - prev) / l).clamp(0.0, 1.0);
            positions.push([a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1])]);
            normals.push([(b[1] - a[1]) / l, -(b[0] - a[0]) / l]);
        }
        Cloud2d { positions, normals }
    }

    /// First hit of the ray `o + t d` (`t > 0`): `(t, segment)`.
    pub fn intersect(&self, o: [f64; 2], d: [f64; 2]) -> Option<(f64, usize)> {
        let mut best: Option<(f64, usize)> = None;
        for (i, &(a, b)) in self.segments.iter().enumerate() {
            let e = [b[0] - a[0], b[1] - a[1]];
            let den = d[0] * e[1] - d[1] * e[0];
            if den.abs() < 1e-300 {
                continue;
            }
            let w = [a[0] - o[0], a[1] - o[1]];
            let t = (w[0] * e[1] - w[1] * e[0]) / den;
            let u = (w[0] * d[1] - w[1] * d[0]) / den;
            if t > 0.0 && (-1e-12..=1.0 + 1e-12).contains(&u) && best.is_none_or(|(bt, bi)| t < bt || (t == bt && i < bi)) {
                best = Some((t, i));
            }
        }
        best
    }

    /// Outward unit normal of segment `i`.
    pub fn segment_normal(&self, i: usize) -> [f64; 2] {
        let (a, b) = self.segments[i];
        let l = (b[0] - a[0]).hypot(b[1] - a[1]);
        [(b[1] - a[1]) / l, -(b[0] - a[0]) / l]
    }
}

/// A random simple (star-shaped) polygon around the origin with `3..=max_vertices` vertices, radii
/// in `[0.3, 1]`, centered in `[-1, 1]²`.
pub fn random_polygon(rng: &mut Rng, max_vertices: usize) -> Result<Curve2d, PatError> {
    let n = 3 + rng.below(max_vertices.saturating_sub(2).max(1));
    // sorted angles with a minimum gap so that the polygon is simple and not degenerate
    let mut ang: Vec<f64> = (0..n).map(|_| rng.uniform(0.0, TAU)).collect();
    ang.sort_by(f64::total_cmp);
    let gap_ok = (0..n).all(|i| (ang[(i + 1) % n] - ang[i]).rem_euclid(TAU) > 0.05 && (ang[(i + 1) % n] - ang[i]).rem_euclid(TAU) < PI - 0.05);
    if !gap_ok {
        // fall back to equally spaced angles with jitter
        for (i, a) in ang.iter_mut().enumerate() {
            *a = TAU * (i as f64 + rng.uniform(0.1, 0.9)) / n as f64;
        }
    }
    let v: Vec<[f64; 2]> = ang
        .iter()
        .map(|&a| {
            let r = rng.uniform(0.3, 1.0);
            [r * a.cos(), r * a.sin()]
        })
        .collect();
    Curve2d::polygon(&v)
}

/// Parameters of a 2D noise blob: `F(x) = ‖x‖² − (1 + 2α n(scale x))² = 0`.
#[derive(Debug, Clone, PartialEq)]
pub struct Blob2dParams {
    /// Amplitude `α ∈ (0.1, 0.7)`.
    pub alpha: f64,
    /// Noise scale `∈ (0.5, 2)`.
    pub scale: f64,
    /// Octaves `1..=4`.
    pub octaves: usize,
    /// Noise field number.
    pub noise_base: u64,
}

/// A 2D noise blob traced on a `cells`-cell grid (default 384) over `±1.2(1 + α)`, centered and
/// scaled to `[-1, 1]²`. Fails when the contour is empty or encloses no area.
pub fn noise_blob(params: &Blob2dParams, cells: usize) -> Result<Curve2d, PatError> {
    let noise = Perlin::new(params.noise_base);
    let (alpha, scale, oct) = (params.alpha, params.scale, params.octaves);
    let f = |p: [f64; 2]| {
        let r = 1.0 + 2.0 * alpha * noise.fbm2([p[0] * scale, p[1] * scale], oct, 0.5, 2.0);
        p[0] * p[0] + p[1] * p[1] - r * r
    };
    // F normalized by its gradient: a first-order distance, so Newton projection works.
    let eval = |p: [f64; 2]| {
        let h = 1e-5;
        let g = [(f([p[0] + h, p[1]]) - f([p[0] - h, p[1]])) / (2.0 * h), (f([p[0], p[1] + h]) - f([p[0], p[1] - h])) / (2.0 * h)];
        let l = g[0].hypot(g[1]).max(1e-12);
        (f(p) / l, [g[0] / l, g[1] / l])
    };
    let e = 1.2 * (1.0 + alpha);
    let b = Boundary2d::trace(&eval, [[-e, -e], [e, e]], cells.max(32));
    if b.area <= 1e-3 || b.segments.is_empty() {
        return Err(PatError::Format { reason: format!("noise blob {params:?} encloses no area") });
    }
    // center on the vertex mean, scale by the largest coordinate (as for meshes)
    let mut c = [0.0; 2];
    let mut cnt = 0.0;
    for (a, b2) in &b.segments {
        for p in [a, b2] {
            c[0] += p[0];
            c[1] += p[1];
            cnt += 1.0;
        }
    }
    c = [c[0] / cnt, c[1] / cnt];
    let m = b.segments.iter().flat_map(|(a, b2)| [a, b2]).map(|p| (p[0] - c[0]).abs().max((p[1] - c[1]).abs())).fold(0.0, f64::max);
    let s = if m > 1e-10 { 1.0 / m } else { 1.0 };
    Curve2d::new(b.segments.iter().map(|(a, b2)| ([(a[0] - c[0]) * s, (a[1] - c[1]) * s], [(b2[0] - c[0]) * s, (b2[1] - c[1]) * s])).collect())
}

/// A 2D training shape.
#[derive(Debug, Clone)]
pub enum Instance2d {
    /// One of the 16 Quilez SDFs (exact signed distance).
    Analytic(Shape2d),
    /// A random polygon.
    Polygon(Curve2d),
    /// A 2D noise blob.
    Blob(Curve2d, Blob2dParams),
}

impl Instance2d {
    /// Ground-truth signed distance and gradient.
    pub fn signed_distance(&self, p: [f64; 2]) -> (f64, [f64; 2]) {
        match self {
            Instance2d::Analytic(s) => s.eval(p),
            Instance2d::Polygon(c) | Instance2d::Blob(c, _) => c.signed_distance(p),
        }
    }

    /// Ground truth at many points, in parallel (output order = input order).
    pub fn signed_distances(&self, ps: &[[f64; 2]]) -> Vec<f64> {
        ps.par_iter().map(|&p| self.signed_distance(p).0).collect()
    }

    /// A box containing the shape.
    pub fn bbox(&self) -> [[f64; 2]; 2] {
        match self {
            Instance2d::Analytic(s) => s.bbox(),
            Instance2d::Polygon(c) | Instance2d::Blob(c, _) => c.bbox(),
        }
    }

    /// The boundary as a [`Curve2d`] (traced for analytic shapes at 2048 cells).
    pub fn curve(&self) -> Result<Curve2d, PatError> {
        match self {
            Instance2d::Analytic(s) => Curve2d::new(Boundary2d::trace(&|p| s.eval(p), s.bbox(), 2048).segments),
            Instance2d::Polygon(c) | Instance2d::Blob(c, _) => Ok(c.clone()),
        }
    }

    /// `n` boundary points uniform in arclength with exact normals (analytic shapes: the exact
    /// gradient at the projected point; curves: the segment normals).
    pub fn sample_uniform(&self, n: usize, rng: &mut Rng) -> Result<Cloud2d, PatError> {
        match self {
            Instance2d::Analytic(s) => {
                let b = Boundary2d::trace(&|p| s.eval(p), s.bbox(), 2048);
                let ev = |p: [f64; 2]| s.eval(p);
                let mut positions = Vec::with_capacity(n);
                let mut normals = Vec::with_capacity(n);
                for _ in 0..n {
                    let (p, g) = b.sample(&ev, rng.f64()).ok_or_else(|| PatError::Format { reason: "empty boundary".into() })?;
                    positions.push(p);
                    normals.push(g);
                }
                Ok(Cloud2d { positions, normals })
            }
            Instance2d::Polygon(c) | Instance2d::Blob(c, _) => Ok(c.sample_uniform(n, rng)),
        }
    }

    /// A simulated scan (see [`simulate_scan_2d`]).
    pub fn scan(&self, opts: &Scan2dOptions, seed: u64) -> Result<Cloud2d, PatError> {
        simulate_scan_2d(&self.curve()?, opts, seed)
    }
}

/// Draw instance number `index` of the run `seed`: 40 % analytic (uniform over the 16 kinds), 30 %
/// polygons, 30 % noise blobs.
pub fn random_instance_2d(seed: u64, index: u64) -> Result<Instance2d, PatError> {
    let mut last = None;
    for attempt in 0..20u64 {
        let mut rng = Rng::for_item(seed, index.wrapping_mul(32).wrapping_add(attempt));
        let r = rng.f64();
        let res = if r < 0.4 {
            Ok(Instance2d::Analytic(Shape2d::random(ALL_KINDS[rng.below(ALL_KINDS.len())], &mut rng)))
        } else if r < 0.7 {
            random_polygon(&mut rng, 12).map(Instance2d::Polygon)
        } else {
            let params = Blob2dParams { alpha: rng.uniform(0.1, 0.7), scale: rng.uniform(0.5, 2.0), octaves: 1 + rng.below(4), noise_base: index.wrapping_add(attempt) % 256 };
            noise_blob(&params, 384).map(|c| Instance2d::Blob(c, params))
        };
        match res {
            Ok(i) => return Ok(i),
            Err(e) => last = Some(e),
        }
    }
    Err(last.unwrap_or(PatError::Format { reason: "no attempt".into() }))
}

/// Settings of the 2D scanner.
#[derive(Debug, Clone, PartialEq)]
pub struct Scan2dOptions {
    /// Cameras on the circle (equally spaced under a random rotation); default 3.
    pub cameras: usize,
    /// Spacing of the parallel rays, drawn per camera; default `(0.003, 0.03)`.
    pub spacing: (f64, f64),
    /// Maximum points kept; default 512.
    pub max_points: usize,
    /// Probability of subsampling by farthest points rather than at random; default 0.5.
    pub fps_probability: f64,
}

impl Default for Scan2dOptions {
    fn default() -> Self {
        Scan2dOptions { cameras: 3, spacing: (0.003, 0.03), max_points: 512, fps_probability: 0.5 }
    }
}

/// Farthest-point sampling in 2D (start at `start`, ties to the smallest index).
pub fn farthest_point_sampling_2d(points: &[[f64; 2]], m: usize, start: usize) -> Vec<usize> {
    let n = points.len();
    if n == 0 || m >= n {
        return (0..n).collect();
    }
    let mut mind = vec![f64::INFINITY; n];
    let mut chosen = vec![start.min(n - 1)];
    while chosen.len() < m {
        let q = points[*chosen.last().unwrap_or(&0)];
        let mut best = (-1.0, 0usize);
        for (i, p) in points.iter().enumerate() {
            let d = (p[0] - q[0]).powi(2) + (p[1] - q[1]).powi(2);
            if d < mind[i] {
                mind[i] = d;
            }
            if mind[i] > best.0 {
                best = (mind[i], i);
            }
        }
        chosen.push(best.1);
        mind[best.1] = 0.0;
    }
    chosen
}

/// The 2D analog of the simulated scanner: cameras at distance `√2` from the centroid of the
/// curve's vertices, parallel rays across a sensor of side `2√2`, first hits, normals flipped
/// towards the ray, subsample to `max_points`.
pub fn simulate_scan_2d(curve: &Curve2d, opts: &Scan2dOptions, seed: u64) -> Result<Cloud2d, PatError> {
    if opts.cameras == 0 || !(opts.spacing.0 > 0.0 && opts.spacing.1 >= opts.spacing.0) || opts.max_points == 0 {
        return Err(PatError::InvalidOption { what: format!("2D scanner options {opts:?}") });
    }
    let mut rng = Rng::new(seed);
    let (mut c, mut cnt) = ([0.0; 2], 0.0);
    for (a, b) in curve.segments() {
        c[0] += a[0] + b[0];
        c[1] += a[1] + b[1];
        cnt += 2.0;
    }
    c = [c[0] / cnt, c[1] / cnt];
    let rot = rng.uniform(0.0, TAU);
    let radius = 2f64.sqrt();
    let mut hits: Vec<([f64; 2], [f64; 2])> = Vec::new();
    for k in 0..opts.cameras {
        let a = rot + TAU * k as f64 / opts.cameras as f64;
        let pos = [a.cos(), a.sin()];
        let dir = [-pos[0], -pos[1]];
        let side = [-dir[1], dir[0]];
        let spacing = rng.uniform(opts.spacing.0, opts.spacing.1);
        let width = 2.0 * radius;
        let per = ((width / spacing).ceil() as usize).max(1);
        let half = 0.5 * (per as f64 - 1.0) * spacing;
        let got: Vec<([f64; 2], [f64; 2])> = (0..per)
            .into_par_iter()
            .filter_map(|i| {
                let off = (i as f64 + 0.5) * spacing - half;
                let o = [c[0] + radius * pos[0] + off * side[0], c[1] + radius * pos[1] + off * side[1]];
                let (t, seg) = curve.intersect(o, dir)?;
                let n = curve.segment_normal(seg);
                let n = if n[0] * dir[0] + n[1] * dir[1] > 0.0 { [-n[0], -n[1]] } else { n };
                Some(([o[0] + t * dir[0], o[1] + t * dir[1]], n))
            })
            .collect();
        hits.extend(got);
    }
    let keep = hits.len().min(opts.max_points);
    let pts: Vec<[f64; 2]> = hits.iter().map(|h| h.0).collect();
    let idx: Vec<usize> = if hits.len() <= opts.max_points {
        (0..hits.len()).collect()
    } else if rng.f64() < opts.fps_probability {
        farthest_point_sampling_2d(&pts, keep, 0)
    } else {
        let mut all: Vec<usize> = (0..hits.len()).collect();
        rng.shuffle(&mut all);
        all.truncate(keep);
        all
    };
    Ok(Cloud2d { positions: idx.iter().map(|&i| hits[i].0).collect(), normals: idx.iter().map(|&i| hits[i].1).collect() })
}
