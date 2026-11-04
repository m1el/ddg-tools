//! Per-point circles: the predictor interface, the classical weighted least-squares fit, and the
//! conversion of quadratic coefficients into circles (SPEC.md §8.7).

use super::Cloud2;
use crate::spatial::Neighborhoods;
use crate::PatError;
use crate::pool::make_pool;
use rayon::prelude::*;

/// Coefficients of the local height function `Q(s) = a0 + a1 s + a2 s²` of one point, in world
/// units, in the frame of [`local_frame`] (`s` along the tangent `t̂`, `Q` along the input normal).
/// The 2D analogue of the six coefficients `(a00, a01, a10, a11, a02, a20)` of the network.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Quad2 {
    /// Offset along the input normal.
    pub a0: f64,
    /// Slope.
    pub a1: f64,
    /// Quadratic coefficient (`κ = 2 a2 / (1 + a1²)^{3/2}`).
    pub a2: f64,
}

/// The frame of a point: `(t̂, n)` with `t̂ = (−n_y, n_x)` (the normal rotated by +90 degrees).
pub fn local_frame(n: [f64; 2]) -> ([f64; 2], [f64; 2]) {
    ([-n[1], n[0]], n)
}

/// A per-point predictor of [`Quad2`] coefficients: the seam where the learned 2D network of WP9b
/// replaces the classical fit.
///
/// Contract:
/// - [`QuadraticPredictor::neighborhood_size`] is the number of tokens per point (the point itself
///   first, then its nearest neighbors by `(d², index)`); the caller builds
///   [`Neighborhoods`] of that size (fewer if the cloud is smaller).
/// - [`QuadraticPredictor::predict`] returns one [`Quad2`] per point in **world units** (a learned
///   predictor applies its own scale `σ_i` and output scaling, as `pat-output-scaling` does in 3D).
///   Non-finite coefficients are an error naming the rows.
/// - Implementations are `Sync` and must not depend on the thread count.
pub trait QuadraticPredictor: Sync {
    /// Tokens per point.
    fn neighborhood_size(&self) -> usize;

    /// One coefficient triple per point, and the per-point scale `σ_i` used (NaN when unused).
    fn predict(&self, cloud: &Cloud2, nbhds: &Neighborhoods, threads: Option<usize>) -> Result<Predicted2, PatError>;
}

/// What a [`QuadraticPredictor`] returns.
#[derive(Debug, Clone, PartialEq)]
pub struct Predicted2 {
    /// One coefficient triple per point (world units).
    pub coefficients: Vec<Quad2>,
    /// The scale `σ_i` the predictor used per point (NaN if it has none).
    pub sigma: Vec<f64>,
    /// Rows for which the predictor fell back to a simpler model (the classical fit: a line).
    pub fallbacks: usize,
}

/// The classical predictor: weighted least squares of `h` against `s` over the `k` nearest points,
/// with Gaussian weights of the distance to the point.
///
/// [formula: pat2d-circle-fit]
#[derive(Debug, Clone, PartialEq)]
pub struct ClassicalFit {
    /// Points in the fit, self included (default 16).
    pub k: usize,
    /// Gaussian weights `exp(−d²/(2 b²))` with `b = bandwidth · d_k` (`d_k` the distance of the
    /// farthest neighbor in the set). `f64::INFINITY` gives unweighted least squares. Default 0.6.
    pub bandwidth: f64,
}

impl Default for ClassicalFit {
    fn default() -> Self {
        ClassicalFit { k: 16, bandwidth: 0.6 }
    }
}

/// Weighted quadratic fit of one point. Returns coefficients and whether the fit degenerated to a
/// straight line (fewer than three distinct tangent coordinates or a singular system).
pub(crate) fn fit_quadratic(p: [f64; 2], n: [f64; 2], nb: &[[f64; 2]], bandwidth: f64) -> (Quad2, f64, bool) {
    let (t, nn) = local_frame(n);
    let mut dk: f64 = 0.0;
    let loc: Vec<(f64, f64)> = nb
        .iter()
        .map(|q| {
            let d = [q[0] - p[0], q[1] - p[1]];
            dk = dk.max((d[0] * d[0] + d[1] * d[1]).sqrt());
            (d[0] * t[0] + d[1] * t[1], d[0] * nn[0] + d[1] * nn[1])
        })
        .collect();
    let l = if dk > 0.0 { dk } else { 1.0 };
    // normal equations in u = s/L, v = h/L
    let mut m = [[0.0f64; 3]; 3];
    let mut r = [0.0f64; 3];
    for &(s, h) in &loc {
        let (u, v) = (s / l, h / l);
        let w = if bandwidth.is_finite() { (-(u * u + v * v) / (2.0 * bandwidth * bandwidth)).exp() } else { 1.0 };
        let b = [1.0, u, u * u];
        for a in 0..3 {
            for c in 0..3 {
                m[a][c] += w * b[a] * b[c];
            }
            r[a] += w * b[a] * v;
        }
    }
    match solve3(m, r) {
        Some(c) if c.iter().all(|x| x.is_finite()) => (Quad2 { a0: l * c[0], a1: c[1], a2: c[2] / l }, l, false),
        _ => {
            // straight line through the weighted mean (a1 from the 2x2 system, a2 = 0)
            let (s00, s01, s11) = (m[0][0], m[0][1], m[1][1]);
            let det = s00 * s11 - s01 * s01;
            let c = if det.abs() > 1e-14 * s00 * s11.max(1e-300) { [(r[0] * s11 - s01 * r[1]) / det, (s00 * r[1] - s01 * r[0]) / det] } else { [r[0] / s00.max(1e-300), 0.0] };
            (Quad2 { a0: l * c[0], a1: c[1], a2: 0.0 }, l, true)
        }
    }
}

/// Solve a symmetric positive semidefinite 3x3 system by Cholesky with a pivot check.
fn solve3(m: [[f64; 3]; 3], r: [f64; 3]) -> Option<[f64; 3]> {
    let mut l = [[0.0f64; 3]; 3];
    let scale = m[0][0].max(1e-300);
    for i in 0..3 {
        for j in 0..=i {
            let mut s = m[i][j];
            for k in 0..j {
                s -= l[i][k] * l[j][k];
            }
            if i == j {
                if s <= 1e-12 * scale {
                    return None;
                }
                l[i][i] = s.sqrt();
            } else {
                l[i][j] = s / l[j][j];
            }
        }
    }
    let mut y = [0.0; 3];
    for i in 0..3 {
        let mut s = r[i];
        for k in 0..i {
            s -= l[i][k] * y[k];
        }
        y[i] = s / l[i][i];
    }
    let mut x = [0.0; 3];
    for i in (0..3).rev() {
        let mut s = y[i];
        for k in i + 1..3 {
            s -= l[k][i] * x[k];
        }
        x[i] = s / l[i][i];
    }
    Some(x)
}

impl QuadraticPredictor for ClassicalFit {
    fn neighborhood_size(&self) -> usize {
        self.k
    }

    fn predict(&self, cloud: &Cloud2, nbhds: &Neighborhoods, threads: Option<usize>) -> Result<Predicted2, PatError> {
        if self.k < 3 || self.bandwidth.is_nan() || self.bandwidth <= 0.0 {
            return Err(PatError::InvalidOption { what: format!("ClassicalFit k = {}, bandwidth = {}", self.k, self.bandwidth) });
        }
        let pool = make_pool(threads)?;
        let pos = cloud.positions();
        let out: Vec<(Quad2, f64, bool)> = pool.install(|| {
            (0..cloud.len())
                .into_par_iter()
                .map(|i| {
                    let nb: Vec<[f64; 2]> = nbhds.row(i).iter().map(|&j| pos[j as usize]).collect();
                    fit_quadratic(pos[i], cloud.normals()[i], &nb, self.bandwidth)
                })
                .collect()
        });
        let bad: Vec<u32> = out.iter().enumerate().filter(|(_, o)| ![o.0.a0, o.0.a1, o.0.a2].iter().all(|x| x.is_finite())).map(|(i, _)| i as u32).collect();
        if !bad.is_empty() {
            return Err(PatError::NonFiniteOutput { indices: bad });
        }
        let sigma = out.iter().map(|o| o.1).collect();
        let fallbacks = out.iter().filter(|o| o.2).count();
        Ok(Predicted2 { coefficients: out.into_iter().map(|o| o.0).collect(), sigma, fallbacks })
    }
}

/// One primitive of the blend.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum CirclePrimitive {
    /// `g(x) = sign (‖x − c‖ − radius)`.
    Circle {
        /// Center.
        center: [f64; 2],
        /// Radius `1/|κ|`.
        radius: f64,
        /// `+1` iff `κ < 0` (the center is on the side opposite the normal).
        sign: f64,
    },
    /// A half-plane `g(x) = ⟨x − point, normal⟩` (curvature below the floor).
    Line {
        /// A point of the line (`o = p + a0 n`).
        point: [f64; 2],
        /// Unit normal `n*`.
        normal: [f64; 2],
    },
}

/// The circles of a cloud with the quantities they came from.
#[derive(Debug, Clone, PartialEq)]
pub struct Circles {
    /// One primitive per point.
    pub primitives: Vec<CirclePrimitive>,
    /// Signed curvature `κ` per point.
    pub kappa: Vec<f64>,
    /// The contact points `o = p + a0 n` expressed through the fitted graph, `Q*(0)`.
    pub contact: Vec<[f64; 2]>,
    /// The fitted unit normals `n*`.
    pub normal: Vec<[f64; 2]>,
    /// The coefficients the circles came from.
    pub coefficients: Vec<Quad2>,
    /// The scales `σ_i` of the predictor (NaN when unused).
    pub sigma: Vec<f64>,
}

/// Counters of [`fit_circles`] and of the predictor.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Fit2dReport {
    /// Points.
    pub n: usize,
    /// Primitives that are lines (`|κ| σ_i` below the threshold).
    pub lines: usize,
    /// Primitives that are circles.
    pub circles: usize,
    /// Fits that fell back to a straight line because the least-squares system was singular
    /// (classical predictor only; filled in by [`Pat2d::build`]).
    pub degenerate_fits: usize,
    /// Tokens per neighborhood.
    pub tokens: usize,
}

/// Options of [`fit_circles`].
#[derive(Debug, Clone, PartialEq)]
pub struct Fit2dOptions {
    /// A point is a line if `|κ| · s < line_threshold`, where `s` is the predictor's scale `σ_i`
    /// (the farthest-neighbor distance for the classical fit; if `σ_i` is NaN, the cloud's mean
    /// nearest-neighbor spacing is not available and the absolute `|κ| < line_threshold` is used).
    /// This is the SPEC's `κ_floor`, made dimensionless; default 1e-6.
    pub line_threshold: f64,
    /// Threads.
    pub threads: Option<usize>,
}

impl Default for Fit2dOptions {
    fn default() -> Self {
        Fit2dOptions { line_threshold: 1e-6, threads: None }
    }
}

/// Circles from coefficients: `κ = 2a2/(1 + a1²)^{3/2}` with respect to `n* = (−a1 t̂ + n)/√(1+a1²)`
/// (in the frame `(t̂, n)`), center `Q*(0) − (sign/|κ|) n* = o + n*/κ` (`o = p + a0 n`), and
/// `g = sign (‖x − c‖ − 1/|κ|)`, `sign = +1` iff `κ < 0`. Lines below the threshold.
///
/// [formula: pat2d-circle-fit]
/// [formula: pat2d-circle-sdf]
pub fn fit_circles(cloud: &Cloud2, coeffs: &[Quad2], sigma: &[f64], opts: &Fit2dOptions) -> Result<(Circles, Fit2dReport), PatError> {
    let n = cloud.len();
    if coeffs.len() != n || sigma.len() != n {
        return Err(PatError::SizeMismatch { what: format!("{} coefficient rows and {} scales for {n} points", coeffs.len(), sigma.len()) });
    }
    if opts.line_threshold.is_nan() || opts.line_threshold < 0.0 {
        return Err(PatError::InvalidOption { what: format!("Fit2dOptions::line_threshold = {}", opts.line_threshold) });
    }
    let pool = make_pool(opts.threads)?;
    let rows: Vec<(CirclePrimitive, f64, [f64; 2], [f64; 2])> = pool.install(|| {
        (0..n)
            .into_par_iter()
            .map(|i| {
                let (t, nn) = local_frame(cloud.normals()[i]);
                let p = cloud.positions()[i];
                let q = coeffs[i];
                let a = (1.0 + q.a1 * q.a1).sqrt();
                let kappa = 2.0 * q.a2 / (a * a * a);
                // n* = (-a1 t + n)/A, o = p + a0 n
                let ns = [(-q.a1 * t[0] + nn[0]) / a, (-q.a1 * t[1] + nn[1]) / a];
                let o = [p[0] + q.a0 * nn[0], p[1] + q.a0 * nn[1]];
                let scale = if sigma[i].is_finite() && sigma[i] > 0.0 { sigma[i] } else { 1.0 };
                let prim = if !(kappa.abs() * scale >= opts.line_threshold) || kappa == 0.0 {
                    CirclePrimitive::Line { point: o, normal: ns }
                } else {
                    let c = [o[0] + ns[0] / kappa, o[1] + ns[1] / kappa];
                    CirclePrimitive::Circle { center: c, radius: 1.0 / kappa.abs(), sign: if kappa < 0.0 { 1.0 } else { -1.0 } }
                };
                (prim, kappa, o, ns)
            })
            .collect()
    });
    let bad: Vec<u32> = rows
        .iter()
        .enumerate()
        .filter(|(_, r)| {
            !(r.1.is_finite()
                && match r.0 {
                    CirclePrimitive::Circle { center, radius, .. } => center.iter().all(|c| c.is_finite()) && radius.is_finite(),
                    CirclePrimitive::Line { point, normal } => point.iter().chain(&normal).all(|c| c.is_finite()),
                })
        })
        .map(|(i, _)| i as u32)
        .collect();
    if !bad.is_empty() {
        return Err(PatError::NonFiniteTorus { indices: bad });
    }
    let lines = rows.iter().filter(|r| matches!(r.0, CirclePrimitive::Line { .. })).count();
    let report = Fit2dReport { n, lines, circles: n - lines, degenerate_fits: 0, tokens: 0 };
    Ok((
        Circles {
            primitives: rows.iter().map(|r| r.0).collect(),
            kappa: rows.iter().map(|r| r.1).collect(),
            contact: rows.iter().map(|r| r.2).collect(),
            normal: rows.iter().map(|r| r.3).collect(),
            coefficients: coeffs.to_vec(),
            sigma: sigma.to_vec(),
        },
        report,
    ))
}

/// Signed distance of one primitive.
///
/// [formula: pat2d-circle-sdf]
#[inline]
pub fn circle_sdf(prim: &CirclePrimitive, x: [f64; 2]) -> f64 {
    match *prim {
        CirclePrimitive::Circle { center, radius, sign } => {
            let d = [x[0] - center[0], x[1] - center[1]];
            sign * ((d[0] * d[0] + d[1] * d[1]).sqrt() - radius)
        }
        CirclePrimitive::Line { point, normal } => (x[0] - point[0]) * normal[0] + (x[1] - point[1]) * normal[1],
    }
}

/// Gradient of one primitive's SDF; 0 at the center of a circle (counted by the caller via the
/// returned flag).
#[inline]
pub fn circle_sdf_gradient(prim: &CirclePrimitive, x: [f64; 2]) -> ([f64; 2], bool) {
    match *prim {
        CirclePrimitive::Circle { center, sign, .. } => {
            let d = [x[0] - center[0], x[1] - center[1]];
            let l = (d[0] * d[0] + d[1] * d[1]).sqrt();
            if l == 0.0 { ([0.0, 0.0], true) } else { ([sign * d[0] / l, sign * d[1] / l], false) }
        }
        CirclePrimitive::Line { normal, .. } => (normal, false),
    }
}
