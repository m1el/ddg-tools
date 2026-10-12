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
/// The singular values are computed for the Jacobian divided by its largest entry and scaled
/// back, so they are finite whenever they are representable (a layout can be `10¹⁰⁰` times
/// larger or smaller than the surface: flattening a maze to a disk does that). The
/// determinant, a product of two of them, may still overflow to infinity or underflow to
/// zero. A degenerate triangle on the surface (`ℓ_01 = 0` or `area = 0`) gives values that
/// are not finite.
///
/// [formula: quasi-conformal-error]
pub fn singular_values(lengths: [f64; 3], area: f64, uv: [[f64; 2]; 3]) -> (f64, f64, f64) {
    let (sigma_max, sigma_min, det, scale) = scaled_singular_values(lengths, area, uv);
    (sigma_max * scale, sigma_min * scale, det * scale * scale)
}

/// [`singular_values`] of the Jacobian divided by `scale`, its largest entry in absolute
/// value: `(σ_max / scale, σ_min / scale, det / scale², scale)`. All zero if the Jacobian is.
fn scaled_singular_values(lengths: [f64; 3], area: f64, uv: [[f64; 2]; 3]) -> (f64, f64, f64, f64) {
    let [l01, l12, l20] = lengths;
    let x = (l01 * l01 + (l20 - l12) * (l20 + l12)) / (2.0 * l01);
    let y = 2.0 * area / l01;
    let mut js = [(uv[1][0] - uv[0][0]) / l01, (uv[1][1] - uv[0][1]) / l01];
    let mut jt = [
        ((uv[2][0] - uv[0][0]) - x * js[0]) / y,
        ((uv[2][1] - uv[0][1]) - x * js[1]) / y,
    ];
    let scale = js[0]
        .abs()
        .max(js[1].abs())
        .max(jt[0].abs())
        .max(jt[1].abs());
    if scale == 0.0 {
        return (0.0, 0.0, 0.0, 0.0);
    }
    js = js.map(|v| v / scale);
    jt = jt.map(|v| v / scale);
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
    (sigma_max, sigma_min, det, scale)
}

/// The quasi-conformal error `Q = σ_max / σ_min ≥ 1` of a triangle map (paper Sec. 7.4): the
/// ratio of the singular values of its Jacobian. `Q = 1` for a similarity; `+∞` if the image
/// is degenerate.
///
/// [formula: quasi-conformal-error]
pub fn quasi_conformal_error(lengths: [f64; 3], area: f64, uv: [[f64; 2]; 3]) -> f64 {
    // A ratio: the scale of the Jacobian cancels.
    let (sigma_max, sigma_min, _, _) = scaled_singular_values(lengths, area, uv);
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
    /// Number of faces whose image has zero area, or whose distortion cannot be represented
    /// (an area in the plane that overflows). They are left out of the statistics above.
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
        let q = quasi_conformal_error(lengths, area, uv);
        let u = log_scale_factor(area, plane.abs());
        if !(plane != 0.0 && q.is_finite() && u.is_finite()) {
            r.degenerate_triangles += 1;
            return;
        }
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
            self.report.log_scale_deviation = (self.u_square_sum / self.weight - mean_u * mean_u)
                .max(0.0)
                .sqrt();
        }
        self.report
    }
}

/// The interior angle at `apex` of the plane triangle `(apex, b, c)`, signed: positive if
/// `apex, b, c` are in counter-clockwise order. In `(−π, π]`: `atan2(u × v, u · v)` with
/// `u = b − apex`, `v = c − apex`.
///
/// [formula: signed-angle-2d]
#[inline]
pub(crate) fn signed_corner_angle(apex: [f64; 2], b: [f64; 2], c: [f64; 2]) -> f64 {
    let (ux, uy) = (b[0] - apex[0], b[1] - apex[1]);
    let (vx, vy) = (c[0] - apex[0], c[1] - apex[1]);
    (ux * vy - uy * vx).atan2(ux * vx + uy * vy)
}

/// An end point of one segment that is closer to the line of the other than this fraction of
/// the longer of the two segments counts as lying on it.
const TOUCH_TOLERANCE: f64 = 1e-9;

/// Whether the open segments `a0 a1` and `b0 b1` cross in a single interior point, and if so
/// that point. Touching and collinear overlap do not count, and neither do segments that
/// coincide up to rounding (see [`TOUCH_TOLERANCE`]): the two sides of a cut through a flat
/// region are laid on top of each other, which is contact, not a crossing.
///
/// The segments cross properly iff `b0`, `b1` lie strictly on opposite sides of the line
/// `a0 a1` and `a0`, `a1` strictly on opposite sides of the line `b0 b1`.
///
/// [formula: segment-proper-crossing]
fn proper_crossing(a0: [f64; 2], a1: [f64; 2], b0: [f64; 2], b1: [f64; 2]) -> Option<[f64; 2]> {
    let d1 = signed_double_area(a0, a1, b0);
    let d2 = signed_double_area(a0, a1, b1);
    let d3 = signed_double_area(b0, b1, a0);
    let d4 = signed_double_area(b0, b1, a1);
    let length_a = (a1[0] - a0[0]).hypot(a1[1] - a0[1]);
    let length_b = (b1[0] - b0[0]).hypot(b1[1] - b0[1]);
    let reach = TOUCH_TOLERANCE * length_a.max(length_b);
    // d / length is the signed distance from the line.
    let opposite = |x: f64, y: f64, length: f64| {
        let margin = reach * length;
        (x > margin && y < -margin) || (x < -margin && y > margin)
    };
    if opposite(d1, d2, length_a) && opposite(d3, d4, length_b) {
        let t = d3 / (d3 - d4);
        Some([a0[0] + t * (a1[0] - a0[0]), a0[1] + t * (a1[1] - a0[1])])
    } else {
        None
    }
}

/// Counts the pairs of non-adjacent edges of a closed polygon that cross each other in a
/// single interior point. A simple polygon has none. Edges that only touch, or overlap along a
/// line, are not counted, also when they do so only up to rounding (an end point within
/// `1e-9` of the length of the longer edge from the other edge's line counts as on it).
///
/// Uses a uniform grid; the expected cost is linear in the number of edges for polygons with
/// edges of comparable length.
pub fn polygon_self_crossings(polygon: &[[f64; 2]]) -> usize {
    let n = polygon.len();
    if n < 4
        || polygon
            .iter()
            .any(|p| !(p[0].is_finite() && p[1].is_finite()))
    {
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

/// How much of the plane a layout covers, and how much of that it covers more than once:
/// an estimate by point sampling on a regular grid.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Coverage {
    /// Area covered by at least one triangle.
    pub covered_area: f64,
    /// Area covered by two or more triangles (counted once, however many lie on top of each
    /// other). Zero for an injective layout. Flipped triangles lie on top of their neighbours
    /// and are included.
    pub overlapped_area: f64,
    /// Distance between neighbouring sample points. Features of the layout smaller than this
    /// are not resolved; the two areas are multiples of its square.
    pub sample_spacing: f64,
    /// Number of sample points covered by at least one triangle.
    pub n_samples: usize,
}

impl Coverage {
    /// `overlapped_area / covered_area`: the fraction of the footprint of the layout in which
    /// two or more points of the surface share a place in the plane. Zero if nothing is
    /// covered.
    pub fn overlapped_fraction(&self) -> f64 {
        if self.covered_area > 0.0 {
            self.overlapped_area / self.covered_area
        } else {
            0.0
        }
    }
}

/// Largest number of grid cells [`coverage`] uses.
const MAX_COVERAGE_CELLS: f64 = 16_777_216.0;

/// Estimates [`Coverage`] of a set of plane triangles by counting, for every point of a
/// regular grid, the triangles that contain it.
///
/// The grid spacing is chosen such that a triangle of average area contains about
/// `samples_per_triangle` points, but not finer than a grid of 2²⁴ cells over the bounding box.
/// A point on an edge shared by two triangles is counted for one of them (top-left rule), so
/// a triangulation without overlap has `overlapped_area == 0` exactly.
///
/// Returns `None` if there is nothing to measure: no triangles, `samples_per_triangle` not
/// positive, a coordinate that is not finite, or zero total area.
pub fn coverage(triangles: &[[[f64; 2]; 3]], samples_per_triangle: f64) -> Option<Coverage> {
    if triangles.is_empty() || samples_per_triangle.is_nan() || samples_per_triangle <= 0.0 {
        return None;
    }
    let (mut lo, mut hi) = ([f64::INFINITY; 2], [f64::NEG_INFINITY; 2]);
    let mut total = 0.0;
    for t in triangles {
        for p in t {
            if !(p[0].is_finite() && p[1].is_finite()) {
                return None;
            }
            for k in 0..2 {
                lo[k] = lo[k].min(p[k]);
                hi[k] = hi[k].max(p[k]);
            }
        }
        total += 0.5 * signed_double_area(t[0], t[1], t[2]).abs();
    }
    let (width, height) = (hi[0] - lo[0], hi[1] - lo[1]);
    if !(total > 0.0 && total.is_finite() && width > 0.0 && height > 0.0) {
        return None;
    }
    let mut h = (total / (samples_per_triangle * triangles.len() as f64)).sqrt();
    h = h.max((width * height / MAX_COVERAGE_CELLS).sqrt());
    if !(h > 0.0 && h.is_finite()) {
        return None;
    }
    // Sample (i, j) sits at lo + (i + ½, j + ½) h.
    let nx = (width / h).ceil().max(1.0) as usize;
    let ny = (height / h).ceil().max(1.0) as usize;
    let mut count = vec![0u8; nx * ny];
    // Whether an edge with direction (dx, dy) of a counter-clockwise triangle owns the points
    // on it: left edges and horizontal top edges do.
    let owns = |dx: f64, dy: f64| dy < 0.0 || (dy == 0.0 && dx < 0.0);
    for t in triangles {
        let double_area = signed_double_area(t[0], t[1], t[2]);
        if double_area == 0.0 {
            continue;
        }
        // Counter-clockwise order.
        let [a, b, c] = if double_area > 0.0 {
            *t
        } else {
            [t[0], t[2], t[1]]
        };
        // The samples lo + (i + ½) h with first <= i < end lie in [min, max].
        let first = |min: f64, origin: f64, n: usize| {
            (((min - origin) / h - 0.5).ceil().max(0.0) as usize).min(n)
        };
        let end = |max: f64, origin: f64, n: usize| {
            ((((max - origin) / h - 0.5).floor() + 1.0).max(0.0) as usize).min(n)
        };
        let (x0, x1) = (
            first(a[0].min(b[0]).min(c[0]), lo[0], nx),
            end(a[0].max(b[0]).max(c[0]), lo[0], nx),
        );
        let (y0, y1) = (
            first(a[1].min(b[1]).min(c[1]), lo[1], ny),
            end(a[1].max(b[1]).max(c[1]), lo[1], ny),
        );
        let edges = [(a, b), (b, c), (c, a)];
        for j in y0..y1 {
            let y = lo[1] + (j as f64 + 0.5) * h;
            for i in x0..x1 {
                let p = [lo[0] + (i as f64 + 0.5) * h, y];
                let inside = edges.iter().all(|&(u, v)| {
                    let side = signed_double_area(u, v, p);
                    side > 0.0 || (side == 0.0 && owns(v[0] - u[0], v[1] - u[1]))
                });
                if inside {
                    let slot = &mut count[j * nx + i];
                    *slot = slot.saturating_add(1);
                }
            }
        }
    }
    let covered = count.iter().filter(|&&c| c >= 1).count();
    let overlapped = count.iter().filter(|&&c| c >= 2).count();
    Some(Coverage {
        covered_area: covered as f64 * h * h,
        overlapped_area: overlapped as f64 * h * h,
        sample_spacing: h,
        n_samples: covered,
    })
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
        let (c, s) = (
            self.scale * self.rotation.cos(),
            self.scale * self.rotation.sin(),
        );
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
    if norm.is_nan() || norm <= 0.0 {
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
        // The same at sizes whose squares are not representable.
        for factor in [1e-170, 1e120, 1e170] {
            let big = uv.map(|p| p.map(|x| x * factor));
            let (smax, smin, _) = singular_values(lengths, 6.0, big);
            assert!(
                (smax / (2.0 * factor) - 1.0).abs() < 1e-13,
                "{factor}: {smax}"
            );
            assert!(
                (smin / (2.0 * factor) - 1.0).abs() < 1e-13,
                "{factor}: {smin}"
            );
            assert!((quasi_conformal_error(lengths, 6.0, big) - 1.0).abs() < 1e-13);
        }
        assert_eq!(
            singular_values(lengths, 6.0, [[1.0, 1.0]; 3]),
            (0.0, 0.0, 0.0)
        );
        assert_eq!(
            quasi_conformal_error(lengths, 6.0, [[1.0, 1.0]; 3]),
            f64::INFINITY
        );
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
        assert_eq!(
            (r.n_faces, r.flipped_triangles, r.degenerate_triangles),
            (3, 1, 1)
        );
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
        // A slit: a path walked out and back, the way back off by rounding to either side.
        // Contact, not a crossing.
        let mut slit: Vec<[f64; 2]> = vec![[0.0, 0.0], [4.0, 0.0], [4.0, 3.0], [2.0, 3.0]];
        let out = [[2.0, 2.0], [2.3, 1.4], [1.9, 1.0]];
        slit.extend(out);
        slit.push([2.1, 0.5]);
        slit.extend(
            out.iter()
                .rev()
                .enumerate()
                .map(|(k, p)| [p[0] + 1e-15 * (k as f64 - 1.0), p[1] - 2e-16]),
        );
        slit.extend([[2.0 - 1e-16, 3.0], [0.0, 3.0]]);
        assert_eq!(polygon_self_crossings(&slit), 0);
    }

    #[test]
    fn coverage_measures_overlap() {
        // Two triangles that make the unit square: no overlap, area 1, also when the sample
        // points fall on the shared diagonal and on the outer edges.
        let square = [
            [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0]],
            [[0.0, 0.0], [1.0, 1.0], [0.0, 1.0]],
        ];
        for samples in [2.0, 8.0, 50.0, 4096.0] {
            let c = coverage(&square, samples).unwrap();
            assert_eq!(c.overlapped_area, 0.0, "{samples}: {c:?}");
            // The covered area is right up to the samples along the perimeter.
            assert!(
                (c.covered_area - 1.0).abs() <= 4.0 * c.sample_spacing,
                "{samples}: {c:?}"
            );
            assert_eq!(c.overlapped_fraction(), 0.0);
        }
        // The second triangle folded over onto the first (a flip): everything is covered
        // twice, and the footprint is half the square.
        let folded = [
            [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0]],
            [[0.0, 0.0], [1.0, 1.0], [1.0, 0.0]],
        ];
        let c = coverage(&folded, 4096.0).unwrap();
        assert!((c.covered_area - 0.5).abs() < 0.02, "{c:?}");
        assert!((c.overlapped_fraction() - 1.0).abs() < 1e-12, "{c:?}");
        // Two unit squares that share a quarter of their area: 0.25 of 1.75.
        let shifted = |d: f64| {
            [
                [[d, d], [1.0 + d, d], [1.0 + d, 1.0 + d]],
                [[d, d], [1.0 + d, 1.0 + d], [d, 1.0 + d]],
            ]
        };
        let mut both = square.to_vec();
        both.extend(shifted(0.5));
        let c = coverage(&both, 10_000.0).unwrap();
        assert!((c.covered_area - 1.75).abs() < 0.01, "{c:?}");
        assert!((c.overlapped_area - 0.25).abs() < 0.01, "{c:?}");
        assert!((c.overlapped_fraction() - 0.25 / 1.75).abs() < 0.01);
        assert!(c.sample_spacing > 0.0 && c.n_samples > 10_000);
        // Nothing to measure.
        assert!(coverage(&[], 16.0).is_none());
        assert!(coverage(&square, 0.0).is_none());
        assert!(coverage(&[[[0.0, 0.0], [1.0, 0.0], [2.0, 0.0]]], 16.0).is_none());
        assert!(coverage(&[[[0.0, f64::NAN], [1.0, 0.0], [2.0, 1.0]]], 16.0).is_none());
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
