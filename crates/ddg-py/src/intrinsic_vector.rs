//! `smoothest_direction_field` and `adaptive_solve`: the tangent fields (paper Sec. 5.5) and the
//! adaptive refinement (Sec. 5.4) of `ddg-intrinsic-apps`, for a triangle mesh given as arrays.
//!
//! Plain Rust: the PyO3 glue (`python/intrinsic_vector.rs`) decodes the buffers and calls
//! [`direction_field`] and [`adaptive`]. Results are per input vertex; the intrinsic
//! triangulation the work ran on stays inside.

use std::time::Instant;

use ddg_core::{Geometry, TriMesh};
use ddg_intrinsic::{IntrinsicTriangulation, RefineOptions, SurfacePoint};
use ddg_intrinsic_apps::amr::{self, AmrOptions, AmrStatus, Problem};
use ddg_intrinsic_apps::vector::{
    self, ConnectionLaplacian, DirectionFieldOptions, EigenStatus,
};

use crate::contract::ReportValue;
use crate::error::{InputError, UnwrapError};

/// A triangle mesh as flat arrays.
#[derive(Debug, Clone, PartialEq)]
pub struct TriangleArrays {
    /// 3 reals per vertex.
    pub positions: Vec<f64>,
    /// 3 vertex indices per triangle.
    pub triangles: Vec<u32>,
}

impl TriangleArrays {
    fn build(&self) -> Result<(TriMesh, Geometry, Vec<[f64; 3]>), UnwrapError> {
        if !self.positions.len().is_multiple_of(3) {
            return Err(InputError::LengthNotMultiple {
                name: "positions",
                len: self.positions.len(),
                width: 3,
                row: "vertex",
            }
            .into());
        }
        if let Some(position) = self.positions.iter().position(|x| !x.is_finite()) {
            return Err(InputError::NonFinite {
                name: "positions",
                position,
                value: self.positions[position],
            }
            .into());
        }
        if !self.triangles.len().is_multiple_of(3) {
            return Err(InputError::LengthNotMultiple {
                name: "triangles",
                len: self.triangles.len(),
                width: 3,
                row: "triangle",
            }
            .into());
        }
        let n = self.positions.len() / 3;
        if let Some(position) = self.triangles.iter().position(|&v| v as usize >= n) {
            return Err(InputError::IndexOutOfRange {
                name: "triangles",
                position,
                value: self.triangles[position],
                count: n,
                of: "vertices",
            }
            .into());
        }
        let positions: Vec<[f64; 3]> = self
            .positions
            .chunks_exact(3)
            .map(|p| [p[0], p[1], p[2]])
            .collect();
        let faces: Vec<[u32; 3]> = self
            .triangles
            .chunks_exact(3)
            .map(|t| [t[0], t[1], t[2]])
            .collect();
        let mesh = TriMesh::new(n, faces).map_err(|e| UnwrapError::Mesh(e.to_string()))?;
        let geometry = Geometry::from_positions(&mesh, &positions)
            .map_err(|e| UnwrapError::Mesh(e.to_string()))?;
        Ok((mesh, geometry, positions))
    }
}

fn bad_option(name: &'static str, value: String, requirement: &'static str) -> UnwrapError {
    InputError::BadOption {
        name,
        value,
        requirement,
    }
    .into()
}

/// Which triangulation a computation runs on.
fn triangulation(
    mesh: &TriMesh,
    geometry: &Geometry,
    which: &str,
    min_angle_degrees: f64,
) -> Result<IntrinsicTriangulation, UnwrapError> {
    let mut tri =
        IntrinsicTriangulation::new(mesh, geometry).map_err(|e| UnwrapError::Mesh(e.to_string()))?;
    let solve = |e: ddg_intrinsic::IntrinsicError| UnwrapError::Solve(e.to_string());
    match which {
        "input" => {}
        "delaunay" => {
            tri.flip_to_delaunay(&Default::default()).map_err(solve)?;
        }
        "refined" => {
            tri.refine(&refine_options(min_angle_degrees)?).map_err(solve)?;
        }
        other => {
            return Err(bad_option(
                "triangulation",
                other.to_string(),
                "must be \"input\", \"delaunay\" or \"refined\"",
            ));
        }
    }
    Ok(tri)
}

fn refine_options(min_angle_degrees: f64) -> Result<RefineOptions, UnwrapError> {
    if !(0.0..=35.0).contains(&min_angle_degrees) {
        return Err(bad_option(
            "min_angle_degrees",
            min_angle_degrees.to_string(),
            "must be between 0 and 35",
        ));
    }
    Ok(RefineOptions {
        min_angle: min_angle_degrees.to_radians(),
        ..Default::default()
    })
}

/// The position in space of a point of the input mesh.
fn position_of(mesh: &TriMesh, positions: &[[f64; 3]], point: SurfacePoint) -> [f64; 3] {
    match point {
        SurfacePoint::Vertex(v) => positions[v as usize],
        SurfacePoint::Edge { edge, t } => {
            let [a, b] = mesh.edge_vertices(edge as usize);
            [0, 1, 2].map(|k| (1.0 - t) * positions[a][k] + t * positions[b][k])
        }
        SurfacePoint::Face { face, bary } => {
            let f = mesh.faces()[face as usize];
            [0, 1, 2].map(|k| (0..3).map(|c| bary[c] * positions[f[c] as usize][k]).sum())
        }
    }
}

/// Options of [`direction_field`].
#[derive(Debug, Clone, PartialEq)]
pub struct FieldOptions {
    /// The symmetry: 1 (vectors), 2 (lines), 4 (crosses), any positive integer.
    pub n: u32,
    /// "input", "delaunay" (default) or "refined".
    pub triangulation: String,
    /// "fem" (Knöppel et al., default) or "fd" (cotan weights and transport).
    pub laplacian: String,
    /// The angle bound of "refined", in degrees.
    pub min_angle_degrees: f64,
}

/// The smoothest field.
#[derive(Debug, Clone, PartialEq)]
pub struct FieldOutput {
    /// 3 per input vertex: one of the field's `n` directions, a unit vector in the plane of an
    /// input face at the vertex (zero where the field vanishes).
    pub vectors: Vec<f64>,
    /// 3 per singularity: its place on the surface (the centre of the intrinsic face).
    pub singularity_points: Vec<f64>,
    /// Per singularity: its index in units of `1/n`.
    pub singularity_indices: Vec<i32>,
    /// Plain values.
    pub report: Vec<(&'static str, ReportValue)>,
}

/// The smoothest `n`-direction field of a triangle mesh (paper Sec. 5.5).
///
/// # Errors
/// [`UnwrapError::Input`] for malformed arrays and options, [`UnwrapError::Mesh`] for a mesh
/// that is not a manifold triangle mesh, [`UnwrapError::Solve`] if the eigenproblem fails.
pub fn direction_field(mesh: &TriangleArrays, options: &FieldOptions) -> Result<FieldOutput, UnwrapError> {
    let start = Instant::now();
    if options.n == 0 {
        return Err(bad_option("n", "0".to_string(), "must be a positive integer"));
    }
    let laplacian = match options.laplacian.as_str() {
        "fem" => ConnectionLaplacian::FiniteElement,
        "fd" => ConnectionLaplacian::FiniteDifference,
        other => {
            return Err(bad_option("laplacian", other.to_string(), "must be \"fem\" or \"fd\""));
        }
    };
    let (tri_mesh, geometry, positions) = mesh.build()?;
    let tri = triangulation(&tri_mesh, &geometry, &options.triangulation, options.min_angle_degrees)?;
    let field_options = DirectionFieldOptions {
        n_sym: options.n,
        laplacian,
        ..Default::default()
    };
    let (field, report) = vector::smoothest_field(&tri, &field_options)
        .map_err(|e| UnwrapError::Solve(e.to_string()))?;
    let normalized = vector::TangentField {
        n_sym: field.n_sym,
        values: field
            .values
            .iter()
            .map(|u| if u.abs() > 0.0 { u.scale(1.0 / u.abs()) } else { *u })
            .collect(),
    };
    let vectors = normalized
        .input_vectors_3d(&tri, &positions)
        .map_err(|e| UnwrapError::Solve(e.to_string()))?;
    let sing = vector::singularities(&tri, &field).map_err(|e| UnwrapError::Solve(e.to_string()))?;
    let mut singularity_points = Vec::with_capacity(3 * sing.faces.len());
    let mut singularity_indices = Vec::with_capacity(sing.faces.len());
    for s in &sing.faces {
        let place = tri
            .point_on_input(SurfacePoint::Face {
                face: s.face,
                bary: [1.0 / 3.0; 3],
            })
            .map_err(|e| UnwrapError::Solve(e.to_string()))?;
        singularity_points.extend(position_of(&tri_mesh, &positions, place));
        singularity_indices.push(s.index);
    }
    let report = vec![
        (
            "status",
            ReportValue::Text(
                match report.status {
                    EigenStatus::Converged => "converged",
                    EigenStatus::IterationLimit => "iteration_limit",
                }
                .to_string(),
            ),
        ),
        ("iterations", ReportValue::Int(report.iterations as i64)),
        ("eigenvalue", ReportValue::Float(report.eigenvalue)),
        ("energy", ReportValue::Float(report.energy)),
        ("residual", ReportValue::Float(report.residual)),
        ("n_triangulation_vertices", ReportValue::Int(tri.n_vertices() as i64)),
        ("index_sum", ReportValue::Int(sing.index_sum)),
        ("euler_characteristic", ReportValue::Int(tri.euler_characteristic())),
        ("seconds", ReportValue::Float(start.elapsed().as_secs_f64())),
    ];
    Ok(FieldOutput {
        vectors: vectors.into_iter().flatten().collect(),
        singularity_points,
        singularity_indices,
        report,
    })
}

/// Options of [`adaptive`].
#[derive(Debug, Clone, PartialEq)]
pub struct AdaptiveOptions {
    /// "green", "heat" or "disk".
    pub problem: String,
    /// The source vertex (green, heat).
    pub source: u32,
    /// The time step of the heat kernel; `None`: the squared mean edge length.
    pub time: Option<f64>,
    /// The vertex budget.
    pub max_vertices: usize,
    /// The bulk fraction of the marking.
    pub bulk: f64,
    /// The angle bound of the refinement, degrees.
    pub min_angle_degrees: f64,
}

/// The solution of [`adaptive`].
#[derive(Debug, Clone, PartialEq)]
pub struct AdaptiveOutput {
    /// `n_components` per input vertex, interleaved.
    pub values: Vec<f64>,
    /// 1 (green, heat) or 2 (disk: the coordinates in the unit disk).
    pub n_components: usize,
    /// Plain values.
    pub report: Vec<(&'static str, ReportValue)>,
}

/// Solves a model problem with adaptive refinement (paper Sec. 5.4) and returns the solution at
/// the input vertices.
///
/// # Errors
/// [`UnwrapError::Input`] for malformed arrays and options (also a source vertex out of range),
/// [`UnwrapError::Mesh`] for a mesh that is not a manifold triangle mesh (or not a disk for
/// "disk"), [`UnwrapError::Solve`] if refinement or a solve fails.
pub fn adaptive(mesh: &TriangleArrays, options: &AdaptiveOptions) -> Result<AdaptiveOutput, UnwrapError> {
    let start = Instant::now();
    let (tri_mesh, geometry, _) = mesh.build()?;
    if options.source as usize >= tri_mesh.n_vertices() && options.problem != "disk" {
        return Err(InputError::IndexOutOfRange {
            name: "source",
            position: 0,
            value: options.source,
            count: tri_mesh.n_vertices(),
            of: "vertices",
        }
        .into());
    }
    let mean_edge =
        geometry.edge_lengths().iter().sum::<f64>() / geometry.edge_lengths().len().max(1) as f64;
    let problem = match options.problem.as_str() {
        "green" => Problem::GreenFunction {
            source: options.source,
        },
        "heat" => {
            let time = options.time.unwrap_or(mean_edge * mean_edge);
            if !(time.is_finite() && time > 0.0) {
                return Err(bad_option("time", time.to_string(), "must be positive and finite"));
            }
            Problem::HeatKernel {
                source: options.source,
                time,
            }
        }
        "disk" => Problem::HarmonicDisk,
        other => {
            return Err(bad_option(
                "problem",
                other.to_string(),
                "must be \"green\", \"heat\" or \"disk\"",
            ));
        }
    };
    if !(options.bulk > 0.0 && options.bulk <= 1.0) {
        return Err(bad_option("bulk", options.bulk.to_string(), "must be in (0, 1]"));
    }
    let amr_options = AmrOptions {
        refine: refine_options(options.min_angle_degrees)?,
        bulk: options.bulk,
        max_vertices: options.max_vertices,
        ..Default::default()
    };
    let mut tri =
        IntrinsicTriangulation::new(&tri_mesh, &geometry).map_err(|e| UnwrapError::Mesh(e.to_string()))?;
    let (solution, report) = amr::adaptive_solve(&mut tri, &problem, &amr_options).map_err(|e| match e {
        amr::AmrError::InvalidInput { reason } => UnwrapError::Mesh(reason),
        other => UnwrapError::Solve(other.to_string()),
    })?;
    let n = tri_mesh.n_vertices();
    let k = solution.components.len();
    let mut values = Vec::with_capacity(k * n);
    for v in 0..n {
        for c in &solution.components {
            values.push(c[v]);
        }
    }
    let last = report.iterations.last().copied();
    let report = vec![
        (
            "status",
            ReportValue::Text(
                match report.status {
                    AmrStatus::Tolerance => "tolerance",
                    AmrStatus::VertexLimit => "vertex_limit",
                    AmrStatus::IterationLimit => "iteration_limit",
                    AmrStatus::Stalled => "stalled",
                }
                .to_string(),
            ),
        ),
        ("iterations", ReportValue::Int(report.iterations.len() as i64)),
        ("n_triangulation_vertices", ReportValue::Int(tri.n_vertices() as i64)),
        ("estimate", ReportValue::Float(last.map_or(f64::NAN, |i| i.estimate))),
        ("oscillation", ReportValue::Float(last.map_or(f64::NAN, |i| i.oscillation))),
        ("seconds", ReportValue::Float(start.elapsed().as_secs_f64())),
    ];
    Ok(AdaptiveOutput {
        values,
        n_components: k,
        report,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn octahedron() -> TriangleArrays {
        TriangleArrays {
            positions: vec![
                1.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 1.0, 0.0,
                0.0, -1.0,
            ],
            triangles: vec![
                0, 2, 4, 2, 1, 4, 1, 3, 4, 3, 0, 4, 2, 0, 5, 1, 2, 5, 3, 1, 5, 0, 3, 5,
            ],
        }
    }

    #[test]
    fn a_field_on_the_octahedron_has_index_sum_two_n() {
        for n in [1u32, 2, 4] {
            let out = direction_field(
                &octahedron(),
                &FieldOptions {
                    n,
                    triangulation: "refined".to_string(),
                    laplacian: "fem".to_string(),
                    min_angle_degrees: 25.0,
                },
            )
            .unwrap();
            assert_eq!(out.vectors.len(), 18);
            let sum = out
                .report
                .iter()
                .find(|(k, _)| *k == "index_sum")
                .map(|(_, v)| v.clone());
            assert_eq!(sum, Some(ReportValue::Int(2 * n as i64)));
            assert_eq!(
                out.singularity_indices.iter().map(|&i| i as i64).sum::<i64>(),
                2 * n as i64
            );
            for v in out.vectors.chunks_exact(3) {
                let len = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
                assert!((len - 1.0).abs() < 1e-9);
            }
        }
    }

    #[test]
    fn bad_options_name_the_parameter() {
        let opts = FieldOptions {
            n: 1,
            triangulation: "nope".to_string(),
            laplacian: "fem".to_string(),
            min_angle_degrees: 25.0,
        };
        let e = direction_field(&octahedron(), &opts).unwrap_err();
        assert!(e.to_string().contains("triangulation"), "{e}");
        let opts = AdaptiveOptions {
            problem: "heat".to_string(),
            source: 99,
            time: None,
            max_vertices: 500,
            bulk: 0.5,
            min_angle_degrees: 25.0,
        };
        let e = adaptive(&octahedron(), &opts).unwrap_err();
        assert!(e.to_string().contains("source"), "{e}");
    }

    #[test]
    fn the_heat_kernel_is_solved_adaptively() {
        let opts = AdaptiveOptions {
            problem: "heat".to_string(),
            source: 4,
            time: Some(0.05),
            max_vertices: 800,
            bulk: 0.5,
            min_angle_degrees: 25.0,
        };
        let out = adaptive(&octahedron(), &opts).unwrap();
        assert_eq!(out.n_components, 1);
        assert_eq!(out.values.len(), 6);
        // The source is the largest value; the opposite pole the smallest.
        assert!(out.values.iter().all(|&x| x <= out.values[4]));
        assert!(out.values.iter().all(|&x| x >= out.values[5]));
    }
}
