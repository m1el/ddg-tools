//! Distortion of a UV map, measured on the result.
//!
//! These figures are computed by the binding from the final per-corner UVs and the input mesh,
//! whatever pipeline produced the UVs. They are a check on the whole chain (pipeline, corner
//! mapping, array hand-over) and the source of the numbers the add-on shows; they are not part
//! of the algorithm.
//!
//! Every triangle of the triangulation is compared with its image in the UV plane, using
//! ddg-core's edge length, triangle area and corner angle functions (formula ids
//! `edge-length-euclidean`, `triangle-area`, `corner-angle`):
//!
//! * **angle error**: `|θ_uv − θ_3d|` per triangle corner. Zero for a conformal map.
//! * **log scale**: `u = ½ ln(A_uv / A_3d)` per triangle, the logarithm of the factor by which
//!   lengths grow (the per-triangle analogue of the log conformal factor `u` of `CLAUDE.md`).
//!   Since a UV map is defined up to a uniform scale, `u` is reported relative to its
//!   area-weighted mean: standard deviation, minimum and maximum. Zero for an isometry.
//!
//! Means are weighted by 3D triangle area. Triangles without area in 3D are skipped and counted;
//! triangles that have area in 3D but none in UV are counted as collapsed and left out of the
//! statistics; triangles whose UV image is wound against the polygon are counted as flipped.

use ddg_core::geometry::{corner_angle, edge_length, triangle_area};

use crate::contract::{Report, UnwrapInput};

/// Measured distortion of a UV map. See the [module documentation](self).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Distortion {
    /// Triangles that entered the statistics.
    pub n_measured: usize,
    /// Triangles with zero area in 3D (skipped).
    pub n_degenerate_triangles: usize,
    /// Triangles with area in 3D whose UV image has zero area (skipped).
    pub n_collapsed_uv_triangles: usize,
    /// Triangles whose UV image has the opposite orientation.
    pub n_flipped_uv_triangles: usize,
    /// Largest corner angle error, radians.
    pub angle_error_max: f64,
    /// Area-weighted mean corner angle error, radians.
    pub angle_error_mean: f64,
    /// Area-weighted standard deviation of the log scale.
    pub log_scale_std: f64,
    /// Smallest log scale minus the mean (≤ 0).
    pub log_scale_min: f64,
    /// Largest log scale minus the mean (≥ 0).
    pub log_scale_max: f64,
}

/// Measures the distortion of `uvs` (one per polygon corner of `input`).
///
/// # Panics
/// If `uvs.len() != input.n_corners()`; callers check that first
/// ([`crate::contract::check_output`]).
pub fn measure_distortion(input: &UnwrapInput, uvs: &[[f64; 2]]) -> Distortion {
    assert_eq!(uvs.len(), input.n_corners(), "one UV per polygon corner");
    let positions = input.positions();
    let corner_vertices = input.corner_vertices();

    let mut n_measured = 0usize;
    let mut n_degenerate = 0usize;
    let mut n_collapsed = 0usize;
    let mut n_flipped = 0usize;
    let mut angle_max = 0.0f64;
    let (mut weight, mut angle_sum, mut u_sum, mut u_sq_sum) = (0.0f64, 0.0f64, 0.0f64, 0.0f64);
    let (mut u_min, mut u_max) = (f64::INFINITY, f64::NEG_INFINITY);

    for tri in input.triangles() {
        let p = tri.map(|c| positions[corner_vertices[c as usize] as usize]);
        let q = tri.map(|c| [uvs[c as usize][0], uvs[c as usize][1], 0.0]);
        // Side k is opposite corner k.
        let sides = |x: &[[f64; 3]; 3]| {
            [
                edge_length(x[1], x[2]),
                edge_length(x[2], x[0]),
                edge_length(x[0], x[1]),
            ]
        };
        let (l3, l2) = (sides(&p), sides(&q));
        let Some(area3) = triangle_area(l3[0], l3[1], l3[2]) else {
            n_degenerate += 1;
            continue;
        };
        let signed =
            (q[1][0] - q[0][0]) * (q[2][1] - q[0][1]) - (q[1][1] - q[0][1]) * (q[2][0] - q[0][0]);
        if signed < 0.0 {
            n_flipped += 1;
        }
        let Some(area2) = triangle_area(l2[0], l2[1], l2[2]) else {
            n_collapsed += 1;
            continue;
        };
        let mut angle_error = 0.0;
        let mut angles_ok = true;
        for k in 0..3 {
            let (i, j) = ((k + 1) % 3, (k + 2) % 3);
            match (
                corner_angle(l3[k], l3[i], l3[j]),
                corner_angle(l2[k], l2[i], l2[j]),
            ) {
                (Some(a3), Some(a2)) => {
                    let e = (a2 - a3).abs();
                    angle_max = angle_max.max(e);
                    angle_error += e / 3.0;
                }
                _ => angles_ok = false,
            }
        }
        if !angles_ok {
            n_collapsed += 1;
            continue;
        }
        let u = 0.5 * (area2 / area3).ln();
        n_measured += 1;
        weight += area3;
        angle_sum += area3 * angle_error;
        u_sum += area3 * u;
        u_sq_sum += area3 * u * u;
        u_min = u_min.min(u);
        u_max = u_max.max(u);
    }

    if n_measured == 0 {
        return Distortion {
            n_measured,
            n_degenerate_triangles: n_degenerate,
            n_collapsed_uv_triangles: n_collapsed,
            n_flipped_uv_triangles: n_flipped,
            angle_error_max: 0.0,
            angle_error_mean: 0.0,
            log_scale_std: 0.0,
            log_scale_min: 0.0,
            log_scale_max: 0.0,
        };
    }
    let mean = u_sum / weight;
    // The variance is a difference of nearly equal numbers when the map is close to an
    // isometry; rounding can make it slightly negative.
    let variance = (u_sq_sum / weight - mean * mean).max(0.0);
    Distortion {
        n_measured,
        n_degenerate_triangles: n_degenerate,
        n_collapsed_uv_triangles: n_collapsed,
        n_flipped_uv_triangles: n_flipped,
        angle_error_max: angle_max,
        angle_error_mean: angle_sum / weight,
        log_scale_std: variance.sqrt(),
        log_scale_min: u_min - mean,
        log_scale_max: u_max - mean,
    }
}

impl Distortion {
    /// Writes the figures into a report under the keys documented in the crate README
    /// (`measured_*`, `n_*_triangles`).
    pub fn write_to(&self, report: &mut Report) {
        report.set("n_degenerate_triangles", self.n_degenerate_triangles);
        report.set("n_collapsed_uv_triangles", self.n_collapsed_uv_triangles);
        report.set("n_flipped_uv_triangles", self.n_flipped_uv_triangles);
        report.set("measured_triangles", self.n_measured);
        report.set("measured_angle_error_max", self.angle_error_max);
        report.set("measured_angle_error_mean", self.angle_error_mean);
        report.set("measured_log_scale_std", self.log_scale_std);
        report.set("measured_log_scale_min", self.log_scale_min);
        report.set("measured_log_scale_max", self.log_scale_max);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::contract::test_meshes::*;

    fn identity_uvs(input: &UnwrapInput) -> Vec<[f64; 2]> {
        input
            .corner_vertices()
            .iter()
            .map(|&v| {
                let p = input.positions()[v as usize];
                [p[0], p[1]]
            })
            .collect()
    }

    #[test]
    fn isometry_has_no_distortion() {
        let input = UnwrapInput::from_buffers(grid(4, 3)).unwrap();
        let d = measure_distortion(&input, &identity_uvs(&input));
        assert_eq!(d.n_measured, 24);
        assert_eq!(
            (
                d.n_degenerate_triangles,
                d.n_collapsed_uv_triangles,
                d.n_flipped_uv_triangles
            ),
            (0, 0, 0)
        );
        assert!(d.angle_error_max < 1e-14);
        assert!(d.log_scale_std < 1e-14);
        assert!(d.log_scale_min.abs() < 1e-14 && d.log_scale_max.abs() < 1e-14);
    }

    #[test]
    fn similarity_is_free_and_a_mirror_image_is_flipped() {
        let input = UnwrapInput::from_buffers(grid(4, 3)).unwrap();
        // Rotate by 90 degrees and scale by 7: still conformal, constant scale.
        let uvs: Vec<[f64; 2]> = identity_uvs(&input)
            .iter()
            .map(|uv| [-7.0 * uv[1], 7.0 * uv[0]])
            .collect();
        let d = measure_distortion(&input, &uvs);
        assert!(d.angle_error_max < 1e-13 && d.log_scale_std < 1e-13);
        assert_eq!(d.n_flipped_uv_triangles, 0);
        let mirrored: Vec<[f64; 2]> = identity_uvs(&input)
            .iter()
            .map(|uv| [-uv[0], uv[1]])
            .collect();
        let d = measure_distortion(&input, &mirrored);
        assert_eq!(d.n_flipped_uv_triangles, 24);
        assert!(d.angle_error_max < 1e-13);
    }

    #[test]
    fn stretching_one_axis_is_measured() {
        // u = 2x, v = y: every unit right triangle of the grid becomes a 2 x 1 right triangle.
        let input = UnwrapInput::from_buffers(grid(2, 2)).unwrap();
        let uvs: Vec<[f64; 2]> = identity_uvs(&input)
            .iter()
            .map(|uv| [2.0 * uv[0], uv[1]])
            .collect();
        let d = measure_distortion(&input, &uvs);
        // The 45 degree corners become atan(1/2) and atan(2).
        let expected = std::f64::consts::FRAC_PI_4 - 0.5f64.atan();
        assert!((d.angle_error_max - expected).abs() < 1e-14);
        // Every triangle doubles its area: constant log scale ln(2)/2, so no deviation.
        assert!(d.log_scale_std < 1e-14);
    }

    #[test]
    fn scale_variation_is_measured() {
        // Two unit quads; the second is mapped at twice the size.
        let input = UnwrapInput::from_buffers(grid(2, 1)).unwrap();
        let mut uvs = identity_uvs(&input);
        for uv in &mut uvs[4..8] {
            *uv = [10.0 + 2.0 * uv[0], 2.0 * uv[1]];
        }
        let d = measure_distortion(&input, &uvs);
        // u = 0 on half the area and ln 2 on the other half.
        let half = 0.5 * std::f64::consts::LN_2;
        assert!((d.log_scale_std - half).abs() < 1e-14);
        assert!((d.log_scale_min + half).abs() < 1e-14);
        assert!((d.log_scale_max - half).abs() < 1e-14);
        assert!(d.angle_error_max < 1e-14);
    }

    #[test]
    fn degenerate_and_collapsed_triangles_are_counted() {
        // A quad with two coincident... a quad whose fourth vertex lies on the edge between the
        // first and third: its second fan triangle (0, 2, 3) has no area.
        let positions = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [2.0, 2.0, 0.0],
            [1.0, 1.0, 0.0],
        ];
        let input = UnwrapInput::from_buffers(buffers(&positions, &[&[0, 1, 2, 3]])).unwrap();
        let d = measure_distortion(&input, &identity_uvs(&input));
        assert_eq!((d.n_measured, d.n_degenerate_triangles), (1, 1));
        // Everything mapped to one point: the remaining triangle collapses.
        let d = measure_distortion(&input, &[[0.5, 0.5]; 4]);
        assert_eq!(
            (
                d.n_measured,
                d.n_degenerate_triangles,
                d.n_collapsed_uv_triangles
            ),
            (0, 1, 1)
        );
        assert_eq!(d.angle_error_max, 0.0);
    }
}
