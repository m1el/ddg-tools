//! The area distortion of a UV map, by the measure the authors of the optimal-cones dataset
//! define in its README.
//!
//! For every triangle `ijk` with area `A` on the surface and `B` in the UV plane the log scale
//! factor is `u_ijk = ln √(B / A)`, and the distortion is `Σ A_ijk u_ijk²`.
//!
//! Two things have to be fixed to make that a number one can compare between maps:
//!
//! * **The unit of area.** Here `A` is divided by the total area of the measured triangles, so
//!   the weights sum to one (the area normalization of the repository's `CLAUDE.md`).
//! * **The scale of the UV map.** Scaling a map by `s` adds `ln s` to every `u`. The published
//!   maps are scaled to fit a texture, and so are packed maps, so the sum as written
//!   ([`AreaDistortion::raw`]) says nothing about them. [`AreaDistortion::scale_free`] is its
//!   minimum over all rescalings of the whole map, `Σ A (u − ū)²` with the area-weighted mean
//!   `ū`: the variance of the log scale factor. It is the figure to compare.
//!
//! The cone solver minimizes `½ Σ_i m_i u_i²` with `u = 0` on the boundary, without removing a
//! mean. Its distortion therefore corresponds to **half** of `raw` on a map at the scale the
//! flattener produced (`e^u`, unpacked), up to the difference between vertex and triangle
//! values; `scale_free ≤ raw`.

use ddg_bff::quality::signed_double_area;

/// The area distortion of a set of triangles.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct AreaDistortion {
    /// `Σ a u²` with `a = A / Σ A` and `u = ½ ln(B / A)`: the README's sum for the map as it is
    /// scaled. Zero if nothing was measured.
    pub raw: f64,
    /// The area-weighted mean `ū = Σ a u`.
    pub mean_log_scale: f64,
    /// `Σ a (u − ū)²`: the README's sum minimized over the scale of the map.
    pub scale_free: f64,
    /// Smallest and largest `u − ū` over the measured triangles (0 if none).
    pub min_log_scale: f64,
    /// See `min_log_scale`.
    pub max_log_scale: f64,
    /// Number of triangles that were measured: those with area on the surface and in the
    /// plane.
    pub n_measured: usize,
    /// Number of triangles without area on the surface (left out).
    pub n_degenerate: usize,
    /// Number of triangles with area on the surface but none in the plane (left out).
    pub n_collapsed: usize,
    /// Number of measured triangles whose image is mirrored (measured with `B = |area|`).
    pub n_flipped: usize,
    /// Total area on the surface of the measured triangles.
    pub surface_area: f64,
    /// Total unsigned area in the plane of the measured triangles.
    pub plane_area: f64,
}

/// Accumulates [`AreaDistortion`] triangle by triangle.
#[derive(Debug, Clone, Default)]
pub(crate) struct Accumulator {
    // (surface area, log scale)
    samples: Vec<(f64, f64)>,
    out: AreaDistortion,
}

impl Accumulator {
    /// Adds a triangle with area `surface` on the surface and signed area `plane` in the
    /// plane.
    pub(crate) fn add(&mut self, surface: f64, plane: f64) {
        if !(surface.is_finite() && surface > 0.0) {
            self.out.n_degenerate += 1;
            return;
        }
        let unsigned = plane.abs();
        if !(unsigned.is_finite() && unsigned > 0.0) {
            self.out.n_collapsed += 1;
            return;
        }
        let u = 0.5 * (unsigned / surface).ln();
        if !u.is_finite() {
            self.out.n_collapsed += 1;
            return;
        }
        if plane < 0.0 {
            self.out.n_flipped += 1;
        }
        self.out.surface_area += surface;
        self.out.plane_area += unsigned;
        self.samples.push((surface, u));
    }

    pub(crate) fn finish(mut self) -> AreaDistortion {
        self.out.n_measured = self.samples.len();
        let total = self.out.surface_area;
        if self.samples.is_empty() || total <= 0.0 {
            return self.out;
        }
        let mean = self.samples.iter().map(|&(a, u)| a * u).sum::<f64>() / total;
        let (mut raw, mut free) = (0.0, 0.0);
        let (mut lo, mut hi) = (f64::INFINITY, f64::NEG_INFINITY);
        for &(a, u) in &self.samples {
            raw += a * u * u;
            free += a * (u - mean) * (u - mean);
            lo = lo.min(u - mean);
            hi = hi.max(u - mean);
        }
        self.out.raw = raw / total;
        self.out.scale_free = free / total;
        self.out.mean_log_scale = mean;
        self.out.min_log_scale = lo;
        self.out.max_log_scale = hi;
        self.out
    }
}

/// The area of the triangle with corners `a`, `b`, `c` in space.
pub(crate) fn triangle_area_3d(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> f64 {
    let u = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
    let v = [c[0] - a[0], c[1] - a[1], c[2] - a[2]];
    let n = [
        u[1] * v[2] - u[2] * v[1],
        u[2] * v[0] - u[0] * v[2],
        u[0] * v[1] - u[1] * v[0],
    ];
    0.5 * (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt()
}

/// Measures a UV map given per polygon corner, on the given triangulation: `triangles` are
/// triples of corner indices, `corner_vertices[c]` is the vertex of corner `c` and `uvs[c]` its
/// UV. Indices that are out of range make the triangle count as degenerate.
///
/// See the [module documentation](self) for the definition.
pub fn area_distortion(
    positions: &[[f64; 3]],
    corner_vertices: &[u32],
    triangles: &[[u32; 3]],
    uvs: &[[f64; 2]],
) -> AreaDistortion {
    let mut accumulator = Accumulator::default();
    for triangle in triangles {
        let corner = |k: usize| -> Option<([f64; 3], [f64; 2])> {
            let c = triangle[k] as usize;
            let v = *corner_vertices.get(c)? as usize;
            Some((*positions.get(v)?, *uvs.get(c)?))
        };
        match (corner(0), corner(1), corner(2)) {
            (Some((p0, q0)), Some((p1, q1)), Some((p2, q2))) => accumulator.add(
                triangle_area_3d(p0, p1, p2),
                0.5 * signed_double_area(q0, q1, q2),
            ),
            _ => accumulator.add(0.0, 0.0),
        }
    }
    accumulator.finish()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_similarity_has_no_scale_free_distortion() {
        let positions = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
        ];
        let corner_vertices = [0, 1, 2, 0, 2, 3];
        let triangles = [[0, 1, 2], [3, 4, 5]];
        // Scaled by 3 and rotated by a quarter turn.
        let uv_of = |v: u32| {
            let p = positions[v as usize];
            [-3.0 * p[1], 3.0 * p[0]]
        };
        let uvs: Vec<[f64; 2]> = corner_vertices.iter().map(|&v| uv_of(v)).collect();
        let d = area_distortion(&positions, &corner_vertices, &triangles, &uvs);
        assert_eq!(d.n_measured, 2);
        assert_eq!(d.n_flipped, 0);
        assert!((d.mean_log_scale - 3f64.ln()).abs() < 1e-14);
        assert!((d.raw - 3f64.ln().powi(2)).abs() < 1e-14);
        assert!(d.scale_free.abs() < 1e-28);
    }

    #[test]
    fn unequal_scales_are_measured_and_flips_counted() {
        let positions = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
        ];
        let corner_vertices = [0, 1, 2, 0, 2, 3];
        let triangles = [[0, 1, 2], [3, 4, 5]];
        // First triangle at scale 1, second at scale e (lengths) and mirrored.
        let e = std::f64::consts::E;
        let uvs = [
            [0.0, 0.0],
            [1.0, 0.0],
            [1.0, 1.0],
            [0.0, 0.0],
            [-e, e],
            [0.0, e],
        ];
        let d = area_distortion(&positions, &corner_vertices, &triangles, &uvs);
        assert_eq!((d.n_measured, d.n_flipped), (2, 1));
        // u = 0 and u = 1 with equal weights: mean ½, raw ½, variance ¼.
        assert!((d.mean_log_scale - 0.5).abs() < 1e-14);
        assert!((d.raw - 0.5).abs() < 1e-14);
        assert!((d.scale_free - 0.25).abs() < 1e-14);
        assert!((d.min_log_scale + 0.5).abs() < 1e-14 && (d.max_log_scale - 0.5).abs() < 1e-14);
    }

    #[test]
    fn triangles_without_area_are_left_out() {
        let positions = [[0.0; 3], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let corner_vertices = [0, 1, 2, 0, 1, 3];
        let triangles = [[0, 1, 2], [3, 4, 5], [3, 4, 9]];
        let uvs = [
            [0.0, 0.0],
            [1.0, 0.0],
            [2.0, 0.0],
            [0.0, 0.0],
            [1.0, 0.0],
            [2.0, 0.0],
        ];
        let d = area_distortion(&positions, &corner_vertices, &triangles, &uvs);
        assert_eq!((d.n_measured, d.n_degenerate, d.n_collapsed), (0, 2, 1));
        assert_eq!(d.raw, 0.0);
    }
}
