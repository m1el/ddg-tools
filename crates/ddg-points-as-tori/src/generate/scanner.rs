//! The simulated scanner (NOTES-REFERENCE §11, mode −1; `sample_even_raycasted_point_cloud`):
//! `cameras` directions by Fibonacci sampling on the unit sphere under a random global rotation,
//! each camera at distance `√3` from the mesh centroid looking at it, a square grid of parallel
//! rays of side `sensor_size · 2√3` with spacing `U(spacing_min, spacing_max)` per camera, the
//! first hit of every ray against a triangle BVH, face normals flipped towards the ray, and a
//! subsample to at most `max_points` by farthest-point sampling or at random (50/50).
//!
//! Per-camera spacing, the rotation, the sensor basis and the subsample draw come from the
//! seed's stream; ray casting is a parallel map whose output order is the ray order, so results
//! do not depend on the thread count.

use super::bvh::TriBvh;
use super::rng::Rng;
use super::sampling::farthest_point_sampling;
use crate::PatError;
use rayon::prelude::*;

/// Scanner settings.
#[derive(Debug, Clone, PartialEq)]
pub struct ScanOptions {
    /// Number of cameras (4).
    pub cameras: usize,
    /// Sensor side as a fraction of the bounding-sphere diameter `2√3` (1).
    pub sensor_size: f64,
    /// Range of the grid spacing (0.01, 0.1).
    pub spacing: (f64, f64),
    /// Maximum points kept (2048).
    pub max_points: usize,
    /// Probability of subsampling by farthest points rather than at random (0.5).
    pub fps_probability: f64,
}

impl Default for ScanOptions {
    fn default() -> Self {
        ScanOptions { cameras: 4, sensor_size: 1.0, spacing: (0.01, 0.1), max_points: 2048, fps_probability: 0.5 }
    }
}

/// A scan: hit points and unit normals facing the camera.
#[derive(Debug, Clone, PartialEq)]
pub struct Scan {
    /// Points on the mesh.
    pub positions: Vec<[f64; 3]>,
    /// Face normals flipped to face the incoming ray.
    pub normals: Vec<[f64; 3]>,
    /// Rays cast.
    pub rays: usize,
    /// Rays that hit.
    pub hits: usize,
}

/// `n` directions by Fibonacci sampling (the reference's `acos(1 − 2i/n)`, golden-angle azimuth).
pub fn fibonacci_sphere(n: usize) -> Vec<[f64; 3]> {
    let golden = 0.5 * (5f64.sqrt() + 1.0);
    (0..n)
        .map(|i| {
            let phi = (1.0 - 2.0 * i as f64 / n as f64).acos();
            let theta = 2.0 * i as f64 * std::f64::consts::PI / golden;
            [phi.sin() * theta.cos(), phi.sin() * theta.sin(), phi.cos()]
        })
        .collect()
}

fn mat_vec(m: &[[f64; 3]; 3], v: [f64; 3]) -> [f64; 3] {
    [0, 1, 2].map(|i| m[i][0] * v[0] + m[i][1] * v[1] + m[i][2] * v[2])
}

/// Simulate one scan of the mesh in `bvh` (stream `seed`).
pub fn simulate_scan(bvh: &TriBvh, opts: &ScanOptions, seed: u64) -> Result<Scan, PatError> {
    if opts.cameras == 0 || !(opts.sensor_size > 0.0) || !(opts.spacing.0 > 0.0 && opts.spacing.1 >= opts.spacing.0) || opts.max_points == 0 {
        return Err(PatError::InvalidOption { what: format!("scanner options {opts:?}") });
    }
    let mut rng = Rng::new(seed);
    let verts = bvh.vertices();
    let mut centroid = [0.0; 3];
    for v in verts {
        for a in 0..3 {
            centroid[a] += v[a] / verts.len() as f64;
        }
    }
    let rot = rng.rotation();
    let diameter = 2.0 * 3f64.sqrt();
    let radius = 0.5 * diameter;
    struct Cam {
        origin: [f64; 3],
        dir: [f64; 3],
        s: [f64; 3],
        t: [f64; 3],
        spacing: f64,
        per_row: usize,
    }
    let mut cams = Vec::new();
    for p in fibonacci_sphere(opts.cameras) {
        let pos = mat_vec(&rot, p);
        let dir = [-pos[0], -pos[1], -pos[2]];
        // random orthonormal basis of the sensor plane
        let r = rng.unit_sphere();
        let dd = r[0] * dir[0] + r[1] * dir[1] + r[2] * dir[2];
        let s = [r[0] - dd * dir[0], r[1] - dd * dir[1], r[2] - dd * dir[2]];
        let sl = (s[0] * s[0] + s[1] * s[1] + s[2] * s[2]).sqrt();
        let s = [s[0] / sl, s[1] / sl, s[2] / sl];
        let t = [dir[1] * s[2] - dir[2] * s[1], dir[2] * s[0] - dir[0] * s[2], dir[0] * s[1] - dir[1] * s[0]];
        let spacing = rng.uniform(opts.spacing.0, opts.spacing.1);
        let per_row = ((opts.sensor_size * diameter / spacing).ceil() as usize).max(1);
        cams.push(Cam { origin: [centroid[0] + radius * pos[0], centroid[1] + radius * pos[1], centroid[2] + radius * pos[2]], dir, s, t, spacing, per_row });
    }
    let mut rays = 0usize;
    let mut hits: Vec<([f64; 3], [f64; 3])> = Vec::new();
    for c in &cams {
        let total = c.per_row * c.per_row;
        rays += total;
        let half = 0.5 * (c.per_row as f64 - 1.0) * c.spacing;
        let got: Vec<([f64; 3], [f64; 3])> = (0..total)
            .into_par_iter()
            .filter_map(|k| {
                let (i, j) = (k / c.per_row, k % c.per_row);
                let (x, y) = ((i as f64 + 0.5) * c.spacing - half, (j as f64 + 0.5) * c.spacing - half);
                let o = [0, 1, 2].map(|a| c.origin[a] + x * c.s[a] + y * c.t[a]);
                let h = bvh.intersect(o, c.dir)?;
                let n = bvh.face_normal(h.face);
                let facing = n[0] * c.dir[0] + n[1] * c.dir[1] + n[2] * c.dir[2] > 0.0;
                let n = if facing { [-n[0], -n[1], -n[2]] } else { n };
                Some(([0, 1, 2].map(|a| o[a] + h.t * c.dir[a]), n))
            })
            .collect();
        hits.extend(got);
    }
    let n_hits = hits.len();
    let keep = n_hits.min(opts.max_points);
    let pts: Vec<[f64; 3]> = hits.iter().map(|h| h.0).collect();
    let idx: Vec<usize> = if n_hits <= opts.max_points {
        (0..n_hits).collect()
    } else if rng.f64() < opts.fps_probability {
        farthest_point_sampling(&pts, keep, 0)
    } else {
        let mut all: Vec<usize> = (0..n_hits).collect();
        rng.shuffle(&mut all);
        all.truncate(keep);
        all
    };
    Ok(Scan { positions: idx.iter().map(|&i| hits[i].0).collect(), normals: idx.iter().map(|&i| hits[i].1).collect(), rays, hits: n_hits })
}
