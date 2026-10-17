//! `map_to_sphere`: a conformal map of a closed genus-0 triangle mesh to the unit sphere.
//!
//! The binding of [`ddg_bff::map_to_sphere`] (Boundary First Flattening, Sawhney & Crane 2017,
//! Sec. 7.5: remove one vertex star, flatten the rest to a round disk, project
//! stereographically, then Möbius-center). The Python glue decodes the buffers into a
//! [`SphereInput`]; [`map_to_sphere`] validates it (every message names the parameter), builds
//! the mesh and the geometry, runs the map and returns a [`SphereOutput`] of plain vectors and a
//! [`SphereReport`]. No Python is involved, so `cargo test -p ddg-py` covers all of it.

use ddg_bff::{BffError, CenteringStatus, PoleChoice, SphereOptions};
use ddg_core::{Geometry, TriMesh};

use crate::error::{InputError, UnwrapError};

/// The decoded arguments of `map_to_sphere`.
#[derive(Debug, Clone)]
pub struct SphereInput {
    /// 3 reals per vertex.
    pub positions: Vec<f64>,
    /// 3 vertex indices per triangle.
    pub triangles: Vec<u32>,
    /// The vertex to remove and send to the pole; `None` picks the vertex with the most regular
    /// star.
    pub pole: Option<u32>,
    /// Whether to Möbius-center the map (on by default).
    pub center: bool,
}

/// What the map reports, in plain values.
#[derive(Debug, Clone, PartialEq)]
pub struct SphereReport {
    /// The vertex that was removed and sent to the pole.
    pub pole: u32,
    /// No flipped or degenerate triangle and covering degree 1.
    pub is_injective: bool,
    /// Triangles whose image is oriented clockwise seen from outside.
    pub flipped_triangles: usize,
    /// Triangles whose image has zero area.
    pub degenerate_triangles: usize,
    /// Total signed area of the image (geodesic triangles) divided by `4π`.
    pub covering_degree: f64,
    /// Quasi-conformal error `σ_max / σ_min` per triangle: area-weighted mean.
    pub quasi_conformal_mean: f64,
    /// Quasi-conformal error: largest.
    pub quasi_conformal_max: f64,
    /// Area distortion: area-weighted rms of the log scale factor about its mean (`0` for a map
    /// that scales all areas alike).
    pub area_distortion: f64,
    /// Smallest log scale factor `½ ln((a/Σa) / (A/ΣA))` over the triangles.
    pub log_scale_min: f64,
    /// Largest log scale factor.
    pub log_scale_max: f64,
    /// Largest `| |x| − 1 |` over the output points.
    pub max_radius_error: f64,
    /// `"converged"`, `"max_iterations"`, `"stalled"`, `"singular_jacobian"`, or `None` when
    /// centering was off.
    pub centering_status: Option<&'static str>,
    /// Number of inversions the centering applied (0 when off).
    pub centering_iterations: usize,
    /// `|µ|`, the distance of the area center of mass from the origin, at the end.
    pub center_of_mass: Option<f64>,
    /// Wall time of the map, seconds.
    pub seconds: f64,
}

/// The result of [`map_to_sphere`].
#[derive(Debug, Clone, PartialEq)]
pub struct SphereOutput {
    /// 3 reals per vertex: a point of the unit sphere.
    pub positions: Vec<f64>,
    /// The report.
    pub report: SphereReport,
}

fn centering_name(status: CenteringStatus) -> &'static str {
    match status {
        CenteringStatus::Converged => "converged",
        CenteringStatus::MaxIterations => "max_iterations",
        CenteringStatus::Stalled => "stalled",
        CenteringStatus::SingularJacobian => "singular_jacobian",
    }
}

/// Maps a map error to the binding's errors.
fn sphere_error(e: BffError) -> UnwrapError {
    match &e {
        BffError::PoleOutOfRange { .. } => UnwrapError::Input(InputError::Rejected {
            name: "pole",
            reason: e.to_string(),
        }),
        BffError::NotASphere { .. } | BffError::Geometry(_) => {
            UnwrapError::Mesh(format!("triangles: {e}"))
        }
        _ => UnwrapError::Solve(e.to_string()),
    }
}

/// Maps the mesh to the unit sphere. See the [module documentation](self).
///
/// # Errors
/// [`UnwrapError::Input`] for malformed arguments (the message starts with the parameter name)
/// and a pole that is not a vertex; [`UnwrapError::Mesh`] for a mesh that is not a valid
/// manifold triangle mesh, has degenerate triangles, or is not a connected closed surface of
/// genus 0; [`UnwrapError::Solve`] for a failure of the computation itself.
pub fn map_to_sphere(input: &SphereInput) -> Result<SphereOutput, UnwrapError> {
    if !input.positions.len().is_multiple_of(3) {
        return Err(InputError::LengthNotMultiple {
            name: "positions",
            len: input.positions.len(),
            width: 3,
            row: "vertex",
        }
        .into());
    }
    if let Some(position) = input.positions.iter().position(|x| !x.is_finite()) {
        return Err(InputError::NonFinite {
            name: "positions",
            position,
            value: input.positions[position],
        }
        .into());
    }
    let n_vertices = input.positions.len() / 3;
    if !input.triangles.len().is_multiple_of(3) {
        return Err(InputError::LengthNotMultiple {
            name: "triangles",
            len: input.triangles.len(),
            width: 3,
            row: "triangle",
        }
        .into());
    }
    if let Some(position) = input
        .triangles
        .iter()
        .position(|&v| v as usize >= n_vertices)
    {
        return Err(InputError::IndexOutOfRange {
            name: "triangles",
            position,
            value: input.triangles[position],
            count: n_vertices,
            of: "vertices",
        }
        .into());
    }
    let faces: Vec<[u32; 3]> = input
        .triangles
        .chunks_exact(3)
        .map(|t| [t[0], t[1], t[2]])
        .collect();
    let positions: Vec<[f64; 3]> = input
        .positions
        .chunks_exact(3)
        .map(|p| [p[0], p[1], p[2]])
        .collect();
    let mesh = TriMesh::new(n_vertices, faces)
        .map_err(|e| UnwrapError::Mesh(format!("triangles: not a manifold mesh: {e}")))?;
    let geometry = Geometry::from_positions(&mesh, &positions)
        .map_err(|e| UnwrapError::Mesh(format!("positions: {e}")))?;
    let defaults = SphereOptions::default();
    let options = SphereOptions {
        pole: input
            .pole
            .map_or(PoleChoice::MostRegularStar, PoleChoice::Vertex),
        centering: if input.center {
            defaults.centering
        } else {
            None
        },
        ..defaults
    };
    let map = ddg_bff::map_to_sphere(&mesh, &geometry, &options).map_err(sphere_error)?;
    let r = &map.report;
    let d = &r.distortion.faces;
    let report = SphereReport {
        pole: r.pole,
        is_injective: r.is_injective(),
        flipped_triangles: d.flipped_triangles,
        degenerate_triangles: d.degenerate_triangles,
        covering_degree: r.distortion.covering_degree,
        quasi_conformal_mean: d.quasi_conformal.mean,
        quasi_conformal_max: d.quasi_conformal.max,
        area_distortion: d.log_scale_deviation,
        log_scale_min: d.log_scale.min,
        log_scale_max: d.log_scale.max,
        max_radius_error: r.distortion.max_radius_error,
        centering_status: r.centering.map(|c| centering_name(c.status)),
        centering_iterations: r.centering.map_or(0, |c| c.iterations),
        center_of_mass: r.centering.map(|c| c.final_center),
        seconds: r.timings.total,
    };
    Ok(SphereOutput {
        positions: map.positions.iter().flatten().copied().collect(),
        report,
    })
}

#[cfg(test)]
mod tests {
    use ddg_core::generators::{grid_disk, icosphere, torus};

    use super::*;

    fn input_of(mesh: &ddg_core::generators::GeneratedMesh) -> SphereInput {
        SphereInput {
            positions: mesh.positions.iter().flatten().copied().collect(),
            triangles: mesh.faces.iter().flatten().copied().collect(),
            pole: None,
            center: true,
        }
    }

    #[test]
    fn icosphere_maps_to_the_unit_sphere() {
        let mesh = icosphere(2).unwrap();
        let out = map_to_sphere(&input_of(&mesh)).unwrap();
        assert_eq!(out.positions.len(), mesh.positions.len() * 3);
        for p in out.positions.chunks_exact(3) {
            let r = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt();
            assert!((r - 1.0).abs() < 1e-12);
        }
        let r = &out.report;
        assert!(r.is_injective);
        assert_eq!(r.flipped_triangles, 0);
        assert_eq!(r.centering_status, Some("converged"));
        assert!(
            r.quasi_conformal_mean < 1.05 && r.area_distortion < 0.05,
            "{r:?}"
        );
        assert!((r.covering_degree - 1.0).abs() < 1e-9);

        let mut input = input_of(&mesh);
        input.pole = Some(5);
        input.center = false;
        let out = map_to_sphere(&input).unwrap();
        assert_eq!(out.report.pole, 5);
        assert_eq!(out.report.centering_status, None);
        assert_eq!(&out.positions[15..18], &[0.0, 0.0, -1.0]);
    }

    #[test]
    fn errors_name_their_cause() {
        let disk = grid_disk(3, 3).unwrap();
        match map_to_sphere(&input_of(&disk)) {
            Err(UnwrapError::Mesh(m)) => assert!(m.contains("genus 0"), "{m}"),
            other => panic!("{other:?}"),
        }
        let t = torus(8, 6, 2.0, 0.5).unwrap();
        match map_to_sphere(&input_of(&t)) {
            Err(UnwrapError::Mesh(m)) => assert!(m.contains("genus 1"), "{m}"),
            other => panic!("{other:?}"),
        }
        let mesh = icosphere(1).unwrap();
        let mut input = input_of(&mesh);
        input.pole = Some(1000);
        match map_to_sphere(&input) {
            Err(UnwrapError::Input(e)) => assert!(e.to_string().starts_with("pole"), "{e}"),
            other => panic!("{other:?}"),
        }
        let mut input = input_of(&mesh);
        input.triangles[4] = 999;
        assert!(matches!(
            map_to_sphere(&input),
            Err(UnwrapError::Input(InputError::IndexOutOfRange { .. }))
        ));
    }
}
