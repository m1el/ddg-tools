//! Measures of a flattening: angle distortion, area distortion, flipped triangles and overlap
//! (Sec. 7.4 of the paper, `Distortion.cpp` of the reference implementation), plus the
//! similarity fit used to compare two layouts.
//!
//! All measures are intrinsic: they compare the edge lengths of a triangle on the surface with
//! its image in the plane, so no vertex positions are needed.

use std::collections::HashMap;

/// Minimum, maximum and mean of a per-face quantity. The mean is weighted by the face areas on
/// the surface.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Stats {
    /// Smallest value.
    pub min: f64,
    /// Largest value.
    pub max: f64,
    /// Area-weighted mean.
    pub mean: f64,
}

/// Twice the signed area of the plane triangle `(a, b, c)`: positive if the corners are in
/// counter-clockwise order.
///
/// [formula: signed-area-2d]
#[inline]
pub fn signed_double_area(a: [f64; 2], b: [f64; 2], c: [f64; 2]) -> f64 {
    (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])
}

/// The singular values `(σ_max, σ_min)` of the affine map that takes a triangle with side
/// lengths `lengths = [ℓ_01, ℓ_12, ℓ_20]` and area `area` to the plane triangle `uv`, and the
/// determinant of that map (negative if the triangle is flipped).
///
/// With the triangle laid out as `p_0 = 0`, `p_1 = (ℓ_01, 0)`, `p_2 = (x, y)`, the columns of
/// the Jacobian are `J_s = (q_1 − q_0)/ℓ_01` and `J_t = ((q_2 − q_0) − x J_s)/y`; with
/// `a = J_s·J_s`, `b = J_s·J_t`, `c = J_t·J_t` the singular values are
/// `σ² = ½ (a + c ± √((a − c)² + 4b²))` (Sander et al. 2001). `σ_min` is evaluated as
/// `|det J| / σ_max` to avoid cancellation.
///
/// [formula: quasi-conformal-error]
pub fn singular_values(lengths: [f64; 3], area: f64, uv: [[f64; 2]; 3]) -> (f64, f64, f64) {
    let [l01, l12, l20] = lengths;
    let x = (l01 * l01 + (l20 - l12) * (l20 + l12)) / (2.0 * l01);
    let y = 2.0 * area / l01;
    let js = [(uv[1][0] - uv[0][0]) / l01, (uv[1][1] - uv[0][1]) / l01];
    let jt = [
        ((uv[2][0] - uv[0][0]) - x * js[0]) / y,
        ((uv[2][1] - uv[0][1]) - x * js[1]) / y,
    ];
    let a = js[0] * js[0] + js[1] * js[1];
    let b = js[0] * jt[0] + js[1] * jt[1];
    let c = jt[0] * jt[0] + jt[1] * jt[1];
    let det = js[0] * jt[1] - js[1] * jt[0];
    let sigma_max = (0.5 * (a + c + ((a - c) * (a - c) + 4.0 * b * b).sqrt())).sqrt();
    let sigma_min = if sigma_max > 0.0 {
        det.abs() / sigma_max
    } else {
        0.0
    };
    (sigma_max, sigma_min, det)
}

/// The quasi-conformal error `Q = σ_max / σ_min ≥ 1` of a triangle map (paper Sec. 7.4): the
/// ratio of the singular values of its Jacobian. `Q = 1` for a similarity; `+∞` if the image
/// is degenerate.
///
/// [formula: quasi-conformal-error]
pub fn quasi_conformal_error(lengths: [f64; 3], area: f64, uv: [[f64; 2]; 3]) -> f64 {
    let (sigma_max, sigma_min, _) = singular_values(lengths, area, uv);
    if sigma_min > 0.0 {
        sigma_max / sigma_min
    } else {
        f64::INFINITY
    }
}

/// The log scale factor of a triangle map, `u = ½ ln(B / A)` for a triangle of area `A` on the
/// surface and `B` in the plane: the logarithm of the factor by which lengths grow on average.
/// (The reference implementation reports `ln(B / A) = 2u`.)
///
/// [formula: face-log-scale-factor]
#[inline]
pub fn log_scale_factor(surface_area: f64, plane_area: f64) -> f64 {
    0.5 * (plane_area / surface_area).ln()
}

/// Distortion of the faces of one chart.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DistortionReport {
    /// Number of faces.
    pub n_faces: usize,
    /// Total area on the surface.
    pub surface_area: f64,
    /// Total signed area in the plane (flipped triangles count negatively).
    pub plane_area: f64,
    /// Quasi-conformal error `Q` over the non-degenerate faces; 1 is perfectly conformal.
    pub quasi_conformal: Stats,
    /// Log scale factor `u = ½ ln(B/A)` over the non-degenerate faces.
    pub log_scale: Stats,
    /// Area-weighted root mean square of `u − mean(u)`: the area distortion that a global
    /// rescaling cannot remove.
    pub log_scale_deviation: f64,
    /// Number of faces whose image has negative orientation.
    pub flipped_triangles: usize,
    /// Total plane area (as a positive number) of the flipped faces.
    pub flipped_area: f64,
    /// Number of faces whose image has zero area.
    pub degenerate_triangles: usize,
}

/// Accumulates [`DistortionReport`]s face by face.
#[derive(Debug, Clone)]
pub(crate) struct DistortionAccumulator {
    report: DistortionReport,
    weight: f64,
    q_sum: f64,
    u_sum: f64,
    u_square_sum: f64,
}

impl DistortionAccumulator {
    pub(crate) fn new() -> Self {
        let empty = Stats {
            min: f64::INFINITY,
            max: f64::NEG_INFINITY,
            mean: f64::NAN,
        };
        DistortionAccumulator {
            report: DistortionReport {
                n_faces: 0,
                surface_area: 0.0,
                plane_area: 0.0,
                quasi_conformal: empty,
                log_scale: empty,
                log_scale_deviation: f64::NAN,
                flipped_triangles: 0,
                flipped_area: 0.0,
                degenerate_triangles: 0,
            },
            weight: 0.0,
            q_sum: 0.0,
            u_sum: 0.0,
            u_square_sum: 0.0,
        }
    }

    /// Adds a face with the given side lengths `[ℓ_01, ℓ_12, ℓ_20]`, area and image.
    pub(crate) fn add(&mut self, lengths: [f64; 3], area: f64, uv: [[f64; 2]; 3]) {
        let r = &mut self.report;
        r.n_faces += 1;
        r.surface_area += area;
        let plane = 0.5 * signed_double_area(uv[0], uv[1], uv[2]);
        r.plane_area += plane;
        if plane < 0.0 {
            r.flipped_triangles += 1;
            r.flipped_area -= plane;
        }
        let (sigma_max, sigma_min, _) = singular_values(lengths, area, uv);
        if !(plane != 0.0 && sigma_min > 0.0 && sigma_max.is_finite()) {
            r.degenerate_triangles += 1;
            return;
        }
        let q = sigma_max / sigma_min;
        let u = log_scale_factor(area, plane.abs());
        r.quasi_conformal.min = r.quasi_conformal.min.min(q);
        r.quasi_conformal.max = r.quasi_conformal.max.max(q);
        r.log_scale.min = r.log_scale.min.min(u);
        r.log_scale.max = r.log_scale.max.max(u);
        self.weight += area;
        self.q_sum += area * q;
        self.u_sum += area * u;
        self.u_square_sum += area * u * u;
    }

    pub(crate) fn finish(mut self) -> DistortionReport {
        if self.weight > 0.0 {
            let mean_u = self.u_sum / self.weight;
            self.report.quasi_conformal.mean = self.q_sum / self.weight;
            self.report.log_scale.mean = mean_u;
            self.report.log_scale_deviation =
                (self.u_square_sum / self.weight - mean_u * mean_u).max(0.0).sqrt();
        }
        self.report
    }
}

/// The interior angle at `apex` of the plane triangle `(apex, b, c)`, signed: positive if
/// `apex, b, c` are in counter-clockwise order. In `(−π, π]`.
#[inline]
pub(crate) fn signed_corner_angle(apex: [f64; 2], b: [f64; 2], c: [f64; 2]) -> f64 {
    let (ux, uy) = (b[0] - apex[0], b[1] - apex[1]);
    let (vx, vy) = (c[0] - apex[0], c[1] - apex[1]);
    (ux * vy - uy * vx).atan2(ux * vx + uy * vy)
}

/// Whether the open segments `a0 a1` and `b0 b1` cross in a single interior point, and if so
/// that point. Touching and collinear overlap do not count.
fn proper_crossing(a0: [f64; 2], a1: [f64; 2], b0: [f64; 2], b1: [f64; 2]) -> Option<[f64; 2]> {
    let d1 = signed_double_area(a0, a1, b0);
    let d2 = signed_double_area(a0, a1, b1);
    let d3 = signed_double_area(b0, b1, a0);
    let d4 = signed_double_area(b0, b1, a1);
    let opposite = |x: f64, y: f64| (x > 0.0 && y < 0.0) || (x < 0.0 && y > 0.0);
    if opposite(d1, d2) && opposite(d3, d4) {
        let t = d3 / (d3 - d4);
        Some([a0[0] + t * (a1[0] - a0[0]), a0[1] + t * (a1[1] - a0[1])])
    } else {
        None
    }
}

/// Counts the pairs of non-adjacent edges of a closed polygon that cross each other in a
/// single interior point. A simple polygon has none. Edges that only touch, or overlap along a
/// line, are not counted.
///
/// Uses a uniform grid; the expected cost is linear in the number of edges for polygons with
/// edges of comparable length.
pub fn polygon_self_crossings(polygon: &[[f64; 2]]) -> usize {
    let n = polygon.len();
    if n < 4 || polygon.iter().any(|p| !(p[0].is_finite() && p[1].is_finite())) {
        return 0;
    }
    let segment = |i: usize| (polygon[i], polygon[(i + 1) % n]);
    // Cell size: the median edge length, but not so small that the polygon needs more than
    // about 8 n cells along its length.
    let mut lengths: Vec<f64> = (0..n)
        .map(|i| {
            let (a, b) = segment(i);
            (b[0] - a[0]).hypot(b[1] - a[1])
        })
        .collect();
    let total: f64 = lengths.iter().sum();
    let (_, median, _) = lengths.select_nth_unstable_by(n / 2, |a, b| a.total_cmp(b));
    let h = median.max(total / (8.0 * n as f64));
    if !(h > 0.0 && h.is_finite()) {
        return 0;
    }
    let cell_of = |p: [f64; 2]| ((p[0] / h).floor() as i64, (p[1] / h).floor() as i64);

    // Register every edge in the cells it passes through (piece by piece, by bounding box).
    let mut cells: HashMap<(i64, i64), Vec<u32>> = HashMap::new();
    for i in 0..n {
        let (a, b) = segment(i);
        let length = (b[0] - a[0]).hypot(b[1] - a[1]);
        let pieces = (length / h).ceil().max(1.0) as usize;
        for k in 0..pieces {
            let t0 = k as f64 / pieces as f64;
            let t1 = (k + 1) as f64 / pieces as f64;
            let p = [a[0] + t0 * (b[0] - a[0]), a[1] + t0 * (b[1] - a[1])];
            let q = [a[0] + t1 * (b[0] - a[0]), a[1] + t1 * (b[1] - a[1])];
            let (lo, hi) = (
                cell_of([p[0].min(q[0]), p[1].min(q[1])]),
                cell_of([p[0].max(q[0]), p[1].max(q[1])]),
            );
            for cx in lo.0..=hi.0 {
                for cy in lo.1..=hi.1 {
                    let list = cells.entry((cx, cy)).or_default();
                    if list.last() != Some(&(i as u32)) {
                        list.push(i as u32);
                    }
                }
            }
        }
    }

    // A crossing is counted in the cell that contains the crossing point.
    let mut count = 0usize;
    for (&cell, list) in &cells {
        for (x, &i) in list.iter().enumerate() {
            for &j in &list[x + 1..] {
                let (i, j) = (i.min(j) as usize, i.max(j) as usize);
                if j == i + 1 || (i == 0 && j == n - 1) || i == j {
                    continue;
                }
                let (a0, a1) = segment(i);
                let (b0, b1) = segment(j);
                if let Some(point) = proper_crossing(a0, a1, b0, b1)
                    && cell_of(point) == cell
                {
                    count += 1;
                }
            }
        }
    }
    count
}

/// The orientation-preserving similarity `z ↦ a z + b` (complex `a`, `b`) that best maps one
/// point set onto another in the least-squares sense, and how well it fits.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Similarity {
    /// `|a|`.
    pub scale: f64,
    /// `arg a`, in radians.
    pub rotation: f64,
    /// `b`.
    pub translation: [f64; 2],
    /// Root mean square of `|a p_i + b − q_i|`.
    pub rms: f64,
    /// Largest `|a p_i + b − q_i|`.
    pub max: f64,
    /// Index of the point with the largest residual.
    pub worst: usize,
}

impl Similarity {
    /// Applies the similarity to a point.
    pub fn apply(&self, p: [f64; 2]) -> [f64; 2] {
        let (c, s) = (self.scale * self.rotation.cos(), self.scale * self.rotation.sin());
        [
            c * p[0] - s * p[1] + self.translation[0],
            s * p[0] + c * p[1] + self.translation[1],
        ]
    }
}

/// Fits the similarity `q ≈ a p + b` to corresponding points: with centroids `p̄`, `q̄`,
/// `a = Σ (q_i − q̄) conj(p_i − p̄) / Σ |p_i − p̄|²` and `b = q̄ − a p̄` (complex arithmetic).
/// Returns `None` if there are no points or all `p_i` coincide.
///
/// [formula: similarity-fit]
pub fn fit_similarity(p: &[[f64; 2]], q: &[[f64; 2]]) -> Option<Similarity> {
    let n = p.len().min(q.len());
    if n == 0 {
        return None;
    }
    let mean = |points: &[[f64; 2]]| {
        let mut m = [0.0, 0.0];
        for point in &points[..n] {
            m[0] += point[0];
            m[1] += point[1];
        }
        [m[0] / n as f64, m[1] / n as f64]
    };
    let (pm, qm) = (mean(p), mean(q));
    let (mut re, mut im, mut norm) = (0.0, 0.0, 0.0);
    for i in 0..n {
        let (px, py) = (p[i][0] - pm[0], p[i][1] - pm[1]);
        let (qx, qy) = (q[i][0] - qm[0], q[i][1] - qm[1]);
        re += qx * px + qy * py;
        im += qy * px - qx * py;
        norm += px * px + py * py;
    }
    if !(norm > 0.0) {
        return None;
    }
    let (ar, ai) = (re / norm, im / norm);
    let translation = [
        qm[0] - (ar * pm[0] - ai * pm[1]),
        qm[1] - (ai * pm[0] + ar * pm[1]),
    ];
    let mut fit = Similarity {
        scale: ar.hypot(ai),
        rotation: ai.atan2(ar),
        translation,
        rms: 0.0,
        max: 0.0,
        worst: 0,
    };
    let mut sum = 0.0;
    for i in 0..n {
        let dx = ar * p[i][0] - ai * p[i][1] + translation[0] - q[i][0];
        let dy = ai * p[i][0] + ar * p[i][1] + translation[1] - q[i][1];
        let d = dx.hypot(dy);
        sum += d * d;
        if d > fit.max {
            fit.max = d;
            fit.worst = i;
        }
    }
    fit.rms = (sum / n as f64).sqrt();
    Some(fit)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn similarity_of_a_triangle_has_q_one() {
        // 3-4-5 triangle, rotated by 30 degrees and scaled by 2.
        let lengths = [3.0, 5.0, 4.0];
        let p = [[0.0, 0.0], [3.0, 0.0], [0.0, 4.0]];
        let (c, s) = (3.0f64.sqrt(), 1.0); // 2 cos 30°, 2 sin 30°
        let uv = p.map(|[x, y]| [c * x - s * y + 7.0, s * x + c * y - 1.0]);
        let (smax, smin, det) = singular_values(lengths, 6.0, uv);
        assert!((smax - 2.0).abs() < 1e-14 && (smin - 2.0).abs() < 1e-14);
        assert!((det - 4.0).abs() < 1e-14);
        assert!((quasi_conformal_error(lengths, 6.0, uv) - 1.0).abs() < 1e-14);
        assert!((log_scale_factor(6.0, 24.0) - 2.0f64.ln()).abs() < 1e-15);
    }

    #[test]
    fn stretch_and_flip_are_measured() {
        let lengths = [1.0, 2.0f64.sqrt(), 1.0];
        // Stretch x by 3: Q = 3.
        let uv = [[0.0, 0.0], [3.0, 0.0], [0.0, 1.0]];
        assert!((quasi_conformal_error(lengths, 0.5, uv) - 3.0).abs() < 1e-14);
        // Mirror image: still Q = 3, determinant negative.
        let flipped = [[0.0, 0.0], [3.0, 0.0], [0.0, -1.0]];
        let (_, _, det) = singular_values(lengths, 0.5, flipped);
        assert!((det + 3.0).abs() < 1e-14);
        let mut acc = DistortionAccumulator::new();
        acc.add(lengths, 0.5, uv);
        acc.add(lengths, 0.5, flipped);
        acc.add(lengths, 0.5, [[0.0, 0.0], [1.0, 0.0], [2.0, 0.0]]);
        let r = acc.finish();
        assert_eq!((r.n_faces, r.flipped_triangles, r.degenerate_triangles), (3, 1, 1));
        assert!((r.flipped_area - 1.5).abs() < 1e-15);
        assert!((r.quasi_conformal.mean - 3.0).abs() < 1e-14);
        assert!((r.log_scale.mean - 0.5 * 3.0f64.ln()).abs() < 1e-14);
        assert!(r.log_scale_deviation < 1e-7);
    }

    #[test]
    fn crossings_of_polygons() {
        let square = [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
        assert_eq!(polygon_self_crossings(&square), 0);
        let bowtie = [[0.0, 0.0], [1.0, 1.0], [1.0, 0.0], [0.0, 1.0]];
        assert_eq!(polygon_self_crossings(&bowtie), 1);
        // A pentagram has five crossings.
        let star: Vec<[f64; 2]> = (0..5)
            .map(|k| {
                let a = 4.0 * std::f64::consts::PI * k as f64 / 5.0;
                [a.cos(), a.sin()]
            })
            .collect();
        assert_eq!(polygon_self_crossings(&star), 5);
        // A fine circle traversed twice has no proper crossing with itself (it overlaps), but a
        // long chord through a fine zig-zag crosses every tooth.
        let mut zigzag: Vec<[f64; 2]> = (0..200)
            .map(|i| [i as f64 * 0.01, if i % 2 == 0 { 0.1 } else { -0.1 }])
            .collect();
        zigzag.push([3.0, 0.0]);
        zigzag.push([-1.0, 0.0]);
        // The closing chord from (-1, 0) back to the first point does not cross anything, the
        // chord from (3, 0) to (-1, 0) crosses all 199 teeth.
        assert_eq!(polygon_self_crossings(&zigzag), 199);
    }

    #[test]
    fn similarity_fit_recovers_a_similarity() {
        let p: Vec<[f64; 2]> = (0..20)
            .map(|i| [(i as f64 * 0.7).sin(), (i as f64 * 1.3).cos() * 2.0])
            .collect();
        let truth = Similarity {
            scale: 1.7,
            rotation: -2.1,
            translation: [0.3, -4.0],
            rms: 0.0,
            max: 0.0,
            worst: 0,
        };
        let q: Vec<[f64; 2]> = p.iter().map(|&x| truth.apply(x)).collect();
        let fit = fit_similarity(&p, &q).unwrap();
        assert!((fit.scale - 1.7).abs() < 1e-13);
        assert!((fit.rotation + 2.1).abs() < 1e-13);
        assert!(fit.max < 1e-13);
        assert!(fit_similarity(&[[1.0, 1.0]; 3], &q[..3]).is_none());
    }
}
