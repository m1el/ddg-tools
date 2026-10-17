//! `reoptimize_cones`: the optimal cone angles for given cone positions.
//!
//! The binding of [`ddg_optimal_cones::reoptimize_cones`]. The Python glue decodes the buffers
//! into a [`ConesInput`]; [`reoptimize`] validates it (every message names the parameter), builds
//! the mesh and the geometry, runs the solver and returns a [`ConesOutput`] of plain vectors and
//! a [`ConesReport`]. No Python is involved, so `cargo test -p ddg-py` covers all of it.
//!
//! The problem is the canonical problem of the repository with the cone cost `0` at the given
//! vertices and `+∞` at every other interior vertex: the least-distortion angles for cones
//! whose positions are fixed. The solver solves the linear optimality system directly, or the
//! reduced system of the free sites when there are few cones, and the generic solver when a
//! bound is active.
//!
//! # Boundary treatments
//!
//! `"fixed"` (`u = 0` on the boundary), `"free"`, `"convex"` and `"penalized"`; see
//! [`ddg_optimal_cones::BoundaryTreatment`]. `"penalized"` charges the boundary curvature
//! `λ·|ν|` with the trade-off `lambda` ([`DEFAULT_LAMBDA`], the default of `unwrap`); the
//! other treatments do not depend on it.

use ddg_core::{Geometry, TriMesh};
use ddg_optimal_cones::{
    BoundaryTreatment, ConeError, ConeProblem, ReoptimizeMethod, ReoptimizeMethodUsed,
    ReoptimizeOptions, Status, reoptimize_cones,
};

use crate::error::{InputError, UnwrapError};

/// The default `λ` for `boundary="penalized"`: the cost of a
/// unit of boundary curvature relative to the distortion. It is `unwrap`'s default.
pub const DEFAULT_LAMBDA: f64 = 0.1;

/// The decoded arguments of `reoptimize_cones`.
#[derive(Debug, Clone, Default)]
pub struct ConesInput {
    /// 3 reals per vertex.
    pub positions: Vec<f64>,
    /// 3 vertex indices per triangle.
    pub triangles: Vec<u32>,
    /// The vertices that carry a cone whose angle is optimized.
    pub cone_vertices: Vec<u32>,
    /// `"fixed"`, `"free"`, `"convex"` or `"penalized"`.
    pub boundary: String,
    /// `"auto"`, `"direct"` or `"solve"`.
    pub method: String,
    /// Cost of the boundary curvature with `"penalized"`: finite, `>= 0`.
    pub lambda: f64,
    /// One weight per vertex (`>= 0`), if given.
    pub distortion_weights: Option<Vec<f64>>,
    /// One lower bound per vertex (`<= 0`), if given.
    pub lower_bounds: Option<Vec<f64>>,
    /// One upper bound per vertex (`>= 0`), if given.
    pub upper_bounds: Option<Vec<f64>>,
    /// Vertices with a fixed cone angle.
    pub prescribed_vertices: Vec<u32>,
    /// The angles of those cones, radians.
    pub prescribed_angles: Vec<f64>,
}

/// What the solver reports, in plain values.
#[derive(Debug, Clone, PartialEq)]
pub struct ConesReport {
    /// `"converged"`, `"not_converged"` or `"linear_solver_failed"`.
    pub status: &'static str,
    /// `"direct_kkt"` or `"solve"`.
    pub method: &'static str,
    /// The distortion `½ Σ wE m u²` (masses divided by the total area).
    pub distortion: f64,
    /// The objective (the distortion, plus the boundary cost when penalized).
    pub objective: f64,
    /// `Σ |angle|` over the given cones.
    pub total_cone_angle: f64,
    /// Number of connected components.
    pub n_components: usize,
    /// Whether the optimality certificate holds at the solver's tolerance (`1e-9`).
    pub certificate_holds: bool,
    /// Largest violation of the state equation.
    pub state_residual: f64,
    /// Primal minus dual objective.
    pub duality_gap: f64,
    /// Solver wall time, seconds.
    pub seconds: f64,
    /// What did not go to plan, if anything.
    pub message: Option<String>,
}

/// The result of [`reoptimize`].
#[derive(Debug, Clone, PartialEq)]
pub struct ConesOutput {
    /// The given cone vertices, in the order given.
    pub cone_vertices: Vec<u32>,
    /// Their optimal angles (radians; may be 0), parallel to `cone_vertices`.
    pub cone_angles: Vec<f64>,
    /// The log scale factor `u`, one per vertex.
    pub log_scale: Vec<f64>,
    /// The boundary curvature `ν` per vertex (0 away from an optimized boundary).
    pub boundary_curvature: Vec<f64>,
    /// The report.
    pub report: ConesReport,
}

fn bad_option(name: &'static str, value: &str, requirement: &'static str) -> UnwrapError {
    InputError::BadOption {
        name,
        value: format!("{value:?}"),
        requirement,
    }
    .into()
}

fn parse_boundary(value: &str) -> Result<BoundaryTreatment, UnwrapError> {
    match value {
        "fixed" => Ok(BoundaryTreatment::Fixed),
        "free" => Ok(BoundaryTreatment::Free),
        "convex" => Ok(BoundaryTreatment::Convex),
        "penalized" => Ok(BoundaryTreatment::Penalized),
        _ => Err(bad_option(
            "boundary",
            value,
            "expected \"fixed\", \"free\", \"convex\" or \"penalized\"",
        )),
    }
}

fn parse_method(value: &str) -> Result<ReoptimizeMethod, UnwrapError> {
    match value {
        "auto" => Ok(ReoptimizeMethod::Auto),
        "direct" => Ok(ReoptimizeMethod::DirectKkt),
        "solve" => Ok(ReoptimizeMethod::Solve),
        _ => Err(bad_option(
            "method",
            value,
            "expected \"auto\", \"direct\" or \"solve\"",
        )),
    }
}

fn check_indices(
    name: &'static str,
    values: &[u32],
    count: usize,
    of: &'static str,
) -> Result<(), InputError> {
    match values.iter().position(|&v| v as usize >= count) {
        Some(position) => Err(InputError::IndexOutOfRange {
            name,
            position,
            value: values[position],
            count,
            of,
        }),
        None => Ok(()),
    }
}

/// An optional per-vertex array must have one entry per vertex; entries must pass `ok`
/// (written so that a NaN fails), otherwise the error says `requirement`.
fn check_per_vertex(
    name: &'static str,
    values: &Option<Vec<f64>>,
    n_vertices: usize,
    requirement: &'static str,
    ok: impl Fn(f64) -> bool,
) -> Result<(), InputError> {
    let Some(values) = values else {
        return Ok(());
    };
    if values.len() != n_vertices {
        return Err(InputError::LengthMismatch {
            name,
            got: values.len(),
            expected: n_vertices,
            reason: "one value per vertex of positions",
        });
    }
    match values.iter().position(|&x| !ok(x)) {
        Some(position) => Err(InputError::BadOption {
            name,
            value: format!("[{position}] = {}", values[position]),
            requirement,
        }),
        None => Ok(()),
    }
}

/// Maps a solver error to the binding's errors, naming the parameter it is about.
fn cone_error(e: ConeError) -> UnwrapError {
    let rejected = |name: &'static str, reason: String| {
        UnwrapError::Input(InputError::Rejected { name, reason })
    };
    match &e {
        ConeError::Geometry(_) => UnwrapError::Mesh(e.to_string()),
        ConeError::InvalidCone { index, .. } => {
            rejected("cone_vertices", format!("entry {index}: {e}"))
        }
        ConeError::InvalidPrescribedCone { index, .. } => {
            rejected("prescribed_vertices", format!("entry {index}: {e}"))
        }
        ConeError::InvalidDistortionWeight { .. } | ConeError::NoDistortionWeight { .. } => {
            rejected("distortion_weights", e.to_string())
        }
        ConeError::InvalidBounds { .. } => rejected("lower_bounds/upper_bounds", e.to_string()),
        ConeError::Infeasible { .. } => rejected(
            "cone_vertices",
            format!(
                "the given cones cannot carry the total angle a closed component needs \
                 (Gauss-Bonnet): {e}"
            ),
        ),
        _ => UnwrapError::Solve(e.to_string()),
    }
}

fn status_name(status: Status) -> &'static str {
    match status {
        Status::Converged => "converged",
        Status::NotConverged => "not_converged",
        Status::LinearSolverFailed => "linear_solver_failed",
    }
}

fn method_name(method: ReoptimizeMethodUsed) -> &'static str {
    match method {
        ReoptimizeMethodUsed::DirectKkt => "direct_kkt",
        ReoptimizeMethodUsed::Solve => "solve",
    }
}

/// Finds the optimal angles of the cones at `input.cone_vertices`. See the
/// [module documentation](self).
///
/// # Errors
/// [`UnwrapError::Input`] for a malformed argument (the message starts with the parameter name),
/// including a cone vertex that is out of range, listed twice, on the boundary or prescribed, and
/// a closed component whose cones cannot carry `2πχ`; [`UnwrapError::Mesh`] for a mesh that
/// `TriMesh::new` or `Geometry::from_positions` rejects; [`UnwrapError::Solve`] for a failure
/// inside the solver that is not about the arguments.
pub fn reoptimize(input: &ConesInput) -> Result<ConesOutput, UnwrapError> {
    let boundary = parse_boundary(&input.boundary)?;
    let method = parse_method(&input.method)?;
    if !(input.lambda.is_finite() && input.lambda >= 0.0) {
        return Err(InputError::BadOption {
            name: "lambda_",
            value: input.lambda.to_string(),
            requirement: "must be finite and >= 0",
        }
        .into());
    }

    // Positions: 3 per vertex, finite (the buffer decoder checked the multiple).
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
    check_indices("triangles", &input.triangles, n_vertices, "vertices")?;
    check_indices(
        "cone_vertices",
        &input.cone_vertices,
        n_vertices,
        "vertices",
    )?;
    check_indices(
        "prescribed_vertices",
        &input.prescribed_vertices,
        n_vertices,
        "vertices",
    )?;
    if input.prescribed_vertices.len() != input.prescribed_angles.len() {
        return Err(InputError::PrescribedLengthMismatch {
            vertices: input.prescribed_vertices.len(),
            angles: input.prescribed_angles.len(),
        }
        .into());
    }
    if let Some(position) = input.prescribed_angles.iter().position(|x| !x.is_finite()) {
        return Err(InputError::NonFinite {
            name: "prescribed_angles",
            position,
            value: input.prescribed_angles[position],
        }
        .into());
    }
    check_per_vertex(
        "distortion_weights",
        &input.distortion_weights,
        n_vertices,
        "weights must be finite and >= 0",
        |x| x.is_finite() && x >= 0.0,
    )?;
    check_per_vertex(
        "lower_bounds",
        &input.lower_bounds,
        n_vertices,
        "lower bounds must be <= 0 (-inf allowed)",
        |x| x <= 0.0,
    )?;
    check_per_vertex(
        "upper_bounds",
        &input.upper_bounds,
        n_vertices,
        "upper bounds must be >= 0 (+inf allowed)",
        |x| x >= 0.0,
    )?;

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

    let mut problem = ConeProblem::new(&mesh, &geometry, input.lambda);
    problem.distortion_weights = input.distortion_weights.clone();
    problem.lower_bounds = input.lower_bounds.clone();
    problem.upper_bounds = input.upper_bounds.clone();
    problem.prescribed_cones = input
        .prescribed_vertices
        .iter()
        .zip(&input.prescribed_angles)
        .map(|(&vertex, &angle)| ddg_optimal_cones::Cone { vertex, angle })
        .collect();
    let options = ReoptimizeOptions {
        boundary,
        method,
        ..ReoptimizeOptions::default()
    };
    let tolerance = options.solve.certificate_tolerance;
    let result = reoptimize_cones(&problem, &input.cone_vertices, &options).map_err(cone_error)?;

    let report = ConesReport {
        status: status_name(result.report.status),
        method: method_name(result.report.method),
        distortion: result.distortion,
        objective: result.objective,
        total_cone_angle: result.total_cone_angle,
        n_components: result.report.components.len(),
        certificate_holds: result.report.certificate.holds(tolerance),
        state_residual: result.report.certificate.state_residual,
        duality_gap: result.report.certificate.duality_gap,
        seconds: result.report.seconds,
        message: result.report.message,
    };
    Ok(ConesOutput {
        cone_vertices: result.cones.iter().map(|c| c.vertex).collect(),
        cone_angles: result.cones.iter().map(|c| c.angle).collect(),
        log_scale: result.log_scale,
        boundary_curvature: result.boundary_curvature,
        report,
    })
}

#[cfg(test)]
mod tests {
    use std::f64::consts::PI;

    use ddg_core::generators::{grid_disk, icosphere};

    use super::*;

    fn input_of(mesh: &ddg_core::generators::GeneratedMesh, cones: &[u32]) -> ConesInput {
        ConesInput {
            positions: mesh.positions.iter().flatten().copied().collect(),
            triangles: mesh.faces.iter().flatten().copied().collect(),
            cone_vertices: cones.to_vec(),
            boundary: "fixed".to_owned(),
            method: "auto".to_owned(),
            lambda: DEFAULT_LAMBDA,
            ..ConesInput::default()
        }
    }

    fn input_error_text(input: &ConesInput) -> String {
        match reoptimize(input) {
            Err(UnwrapError::Input(e)) => e.to_string(),
            other => panic!("expected an InputError, got {other:?}"),
        }
    }

    #[test]
    fn sphere_with_three_cones_sums_to_four_pi() {
        let mesh = icosphere(1).unwrap();
        let out = reoptimize(&input_of(&mesh, &[0, 5, 17])).unwrap();
        assert_eq!(out.cone_vertices, vec![0, 5, 17]);
        assert_eq!(out.cone_angles.len(), 3);
        let sum: f64 = out.cone_angles.iter().sum();
        assert!((sum - 4.0 * PI).abs() < 1e-9, "sum {sum}");
        assert_eq!(out.log_scale.len(), mesh.positions.len());
        assert_eq!(out.boundary_curvature.len(), mesh.positions.len());
        assert_eq!(out.report.status, "converged");
        assert!(out.report.certificate_holds);
        assert_eq!(out.report.n_components, 1);
        assert!(out.report.distortion >= 0.0);
        // The same answer from both methods.
        for method in ["direct", "solve"] {
            let mut input = input_of(&mesh, &[0, 5, 17]);
            input.method = method.to_owned();
            let other = reoptimize(&input).unwrap();
            for (a, b) in out.cone_angles.iter().zip(&other.cone_angles) {
                assert!((a - b).abs() < 1e-7, "{method}: {a} vs {b}");
            }
        }
    }

    #[test]
    fn a_single_cone_carries_everything() {
        let mesh = icosphere(1).unwrap();
        let out = reoptimize(&input_of(&mesh, &[3])).unwrap();
        assert!((out.cone_angles[0] - 4.0 * PI).abs() < 1e-9);
    }

    #[test]
    fn disk_with_fixed_boundary() {
        let mesh = grid_disk(6, 6).unwrap();
        // Vertex (3, 3) is interior.
        let center = 3 * 7 + 3;
        let out = reoptimize(&input_of(&mesh, &[center])).unwrap();
        // A flat disk has no curvature to compensate: the angle is 0 and u = 0.
        assert!(out.cone_angles[0].abs() < 1e-9);
        assert!(out.log_scale.iter().all(|u| u.abs() < 1e-9));
        assert_eq!(out.report.status, "converged");
        // Free boundary and prescribed cone: the log scale is no longer trivial.
        let mut input = input_of(&mesh, &[center]);
        input.prescribed_vertices = vec![2 * 7 + 2];
        input.prescribed_angles = vec![0.5];
        let out = reoptimize(&input).unwrap();
        assert!(out.log_scale.iter().any(|u| u.abs() > 1e-6));
        for boundary in ["free", "convex", "penalized"] {
            input.boundary = boundary.to_owned();
            let out = reoptimize(&input).unwrap();
            assert_eq!(out.log_scale.len(), mesh.positions.len(), "{boundary}");
        }
    }

    #[test]
    fn disconnected_union_equals_the_parts() {
        let a = icosphere(1).unwrap();
        let mut union = input_of(&a, &[0, 7]);
        let offset = a.positions.len() as u32;
        let b = icosphere(1).unwrap();
        union
            .positions
            .extend(b.positions.iter().flatten().map(|x| x + 5.0));
        union
            .triangles
            .extend(b.faces.iter().flatten().map(|v| v + offset));
        union
            .cone_vertices
            .extend([offset + 2, offset + 9, offset + 11]);
        let out = reoptimize(&union).unwrap();
        assert_eq!(out.report.n_components, 2);
        let alone_a = reoptimize(&input_of(&a, &[0, 7])).unwrap();
        let alone_b = reoptimize(&input_of(&b, &[2, 9, 11])).unwrap();
        let expected: Vec<f64> = alone_a
            .cone_angles
            .iter()
            .chain(&alone_b.cone_angles)
            .copied()
            .collect();
        for (x, y) in out.cone_angles.iter().zip(&expected) {
            assert!((x - y).abs() < 1e-8, "{x} vs {y}");
        }
    }

    #[test]
    fn invalid_cone_vertices() {
        let mesh = icosphere(1).unwrap();
        let n = mesh.positions.len() as u32;
        let text = input_error_text(&input_of(&mesh, &[1, n]));
        assert!(text.starts_with("cone_vertices[1]"), "{text}");
        assert!(text.contains("out of range"), "{text}");
        let text = input_error_text(&input_of(&mesh, &[4, 2, 4]));
        assert!(text.starts_with("cone_vertices: "), "{text}");
        assert!(text.contains("more than once"), "{text}");
        // Boundary vertex of a disk.
        let disk = grid_disk(4, 4).unwrap();
        let text = input_error_text(&input_of(&disk, &[0]));
        assert!(text.starts_with("cone_vertices: "), "{text}");
        assert!(text.contains("boundary"), "{text}");
        // A prescribed vertex cannot also be an optimized cone.
        let mut input = input_of(&mesh, &[4]);
        input.prescribed_vertices = vec![4];
        input.prescribed_angles = vec![0.1];
        let text = input_error_text(&input);
        assert!(text.starts_with("cone_vertices: "), "{text}");
        assert!(text.contains("prescribed"), "{text}");
    }

    #[test]
    fn closed_mesh_without_cones_is_infeasible() {
        let mesh = icosphere(1).unwrap();
        let text = input_error_text(&input_of(&mesh, &[]));
        assert!(text.starts_with("cone_vertices: "), "{text}");
        assert!(text.contains("closed"), "{text}");
        // Bounds that forbid the total angle.
        let mut input = input_of(&mesh, &[0, 1]);
        input.upper_bounds = Some(vec![0.1; mesh.positions.len()]);
        let text = input_error_text(&input);
        assert!(text.contains("cannot carry"), "{text}");
    }

    #[test]
    fn lambda_matters_only_for_a_penalized_boundary() {
        let mesh = grid_disk(6, 6).unwrap();
        let center = 3 * 7 + 3;
        let run = |boundary: &str, lambda: f64| {
            let mut input = input_of(&mesh, &[center]);
            input.boundary = boundary.to_owned();
            input.lambda = lambda;
            reoptimize(&input).unwrap()
        };
        let curvature = |o: &ConesOutput| o.boundary_curvature.iter().map(|x| x.abs()).sum::<f64>();
        let (cheap, dear) = (run("penalized", 1e-4), run("penalized", 10.0));
        assert!(
            curvature(&cheap) > curvature(&dear) + 1e-6,
            "lambda changes the result"
        );
        assert_ne!(cheap.report.objective, dear.report.objective);
        let (a, b) = (run("free", 1e-4), run("free", 10.0));
        for (x, y) in a.log_scale.iter().zip(&b.log_scale) {
            assert!((x - y).abs() < 1e-9, "{x} vs {y}");
        }
    }

    #[test]
    fn invalid_lambda() {
        let mesh = icosphere(1).unwrap();
        for lambda in [-0.5, f64::NAN, f64::INFINITY] {
            let mut input = input_of(&mesh, &[0]);
            input.lambda = lambda;
            assert!(
                input_error_text(&input).starts_with("lambda_ = "),
                "{lambda}"
            );
        }
    }

    #[test]
    fn invalid_options_and_arrays() {
        let mesh = icosphere(1).unwrap();
        let n = mesh.positions.len();
        let mut input = input_of(&mesh, &[0]);
        input.boundary = "round".to_owned();
        assert!(input_error_text(&input).starts_with("boundary = \"round\""));
        let mut input = input_of(&mesh, &[0]);
        input.method = "magic".to_owned();
        assert!(input_error_text(&input).starts_with("method = \"magic\""));

        let mut input = input_of(&mesh, &[0]);
        input.distortion_weights = Some(vec![1.0; n - 1]);
        assert!(input_error_text(&input).starts_with("distortion_weights: got"));
        let mut weights = vec![1.0; n];
        weights[3] = -1.0;
        input.distortion_weights = Some(weights);
        assert!(input_error_text(&input).starts_with("distortion_weights = [3] = -1"));
        let mut input = input_of(&mesh, &[0]);
        input.lower_bounds = Some(vec![1.0; n]);
        assert!(input_error_text(&input).starts_with("lower_bounds"));
        let mut input = input_of(&mesh, &[0]);
        input.upper_bounds = Some(vec![-1.0; n]);
        assert!(input_error_text(&input).starts_with("upper_bounds"));

        let mut input = input_of(&mesh, &[0]);
        input.prescribed_vertices = vec![1, 2];
        input.prescribed_angles = vec![0.1];
        assert!(input_error_text(&input).starts_with("prescribed_vertices has 2"));
        let mut input = input_of(&mesh, &[0]);
        input.prescribed_vertices = vec![1];
        input.prescribed_angles = vec![f64::NAN];
        assert!(input_error_text(&input).starts_with("prescribed_angles[0]"));
        let mut input = input_of(&mesh, &[0]);
        input.prescribed_vertices = vec![1, 1];
        input.prescribed_angles = vec![0.1, 0.1];
        assert!(input_error_text(&input).starts_with("prescribed_vertices: "));

        let mut input = input_of(&mesh, &[0]);
        input.positions[4] = f64::INFINITY;
        assert!(input_error_text(&input).starts_with("positions[4]"));
        let mut input = input_of(&mesh, &[0]);
        input.triangles[2] = n as u32;
        assert!(input_error_text(&input).starts_with("triangles[2]"));
    }

    #[test]
    fn bad_meshes_are_mesh_errors() {
        // Three triangles on one edge.
        let input = ConesInput {
            positions: vec![
                0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.5, 1.0, 0.0, 0.5, -1.0, 0.0, 0.5, 0.0, 1.0,
            ],
            triangles: vec![0, 1, 2, 1, 0, 3, 0, 1, 4],
            cone_vertices: vec![],
            boundary: "fixed".to_owned(),
            method: "auto".to_owned(),
            ..ConesInput::default()
        };
        assert!(matches!(reoptimize(&input), Err(UnwrapError::Mesh(m)) if m.contains("3 faces")));
        // A degenerate triangle.
        let input = ConesInput {
            positions: vec![0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 2.0, 0.0, 0.0],
            triangles: vec![0, 1, 2],
            boundary: "fixed".to_owned(),
            method: "auto".to_owned(),
            ..ConesInput::default()
        };
        assert!(matches!(reoptimize(&input), Err(UnwrapError::Mesh(_))));
        // An unreferenced vertex.
        let mesh = icosphere(0).unwrap();
        let mut input = input_of(&mesh, &[0]);
        input.positions.extend([9.0, 9.0, 9.0]);
        assert!(matches!(reoptimize(&input), Err(UnwrapError::Mesh(_))));
        // No faces.
        let input = ConesInput {
            boundary: "fixed".to_owned(),
            method: "auto".to_owned(),
            ..ConesInput::default()
        };
        assert!(matches!(reoptimize(&input), Err(UnwrapError::Mesh(_))));
    }
}
