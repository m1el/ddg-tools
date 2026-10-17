//! The scale factor and adjoint state of a *given* cone configuration.

use ddg_core::SparseMatrix;
use ddg_core::io::Cone;
use ddg_core::sparse::{Factorization, Method, Parallelism, RefineOptions, Symbolic};

use crate::error::{ConeError, PrescribedConeIssue};
use crate::problem::{BoundaryTreatment, ConeProblem, Discrete, Instance, split_pieces};

/// A Cholesky factorization of the Laplacian block of the unknowns.
///
/// On a closed mesh the Laplacian annihilates constants; then the equation of the first
/// unknown is dropped and its value set to zero, which is a solution whenever the right-hand
/// side sums to zero.
pub(crate) struct LaplaceSolver {
    matrix: SparseMatrix,
    factorization: Factorization,
    closed: bool,
}

impl LaplaceSolver {
    pub fn new(d: &Discrete, par: Parallelism<'_>) -> Result<LaplaceSolver, ConeError> {
        let matrix = if d.closed {
            let rest: Vec<u32> = (1..d.n as u32).collect();
            d.a.submatrix(&rest, &rest)?
        } else {
            d.a.clone()
        };
        let factorization = Symbolic::analyze(&matrix, Method::Llt)?.factor(&matrix, par)?;
        Ok(LaplaceSolver {
            matrix,
            factorization,
            closed: d.closed,
        })
    }

    /// Solves `A x = b`. On a closed mesh `b` must sum to zero and the solution with `x_0 = 0`
    /// is returned.
    pub fn solve(&self, b: &[f64], par: Parallelism<'_>) -> Result<Vec<f64>, ConeError> {
        let skip = usize::from(self.closed);
        let rhs = &b[skip..];
        let mut x = vec![0.0; b.len()];
        if !rhs.is_empty() {
            let refine = RefineOptions {
                max_steps: 2,
                tolerance: 1e-15,
            };
            self.factorization
                .solve_refined(&self.matrix, rhs, &mut x[skip..], refine, par)?;
        }
        Ok(x)
    }
}

/// The scale factor and adjoint state that a cone configuration determines.
#[derive(Debug, Clone, PartialEq)]
pub struct ConeEvaluation {
    /// The log scale factor `u` per vertex: the solution of `(L u)_i = µ_i − Ω_i` at interior
    /// vertices with `u = 0` on the boundary. On a closed mesh the additive constant is the
    /// one that minimizes the distortion (`Σ wE m u = 0`).
    pub log_scale: Vec<f64>,
    /// The adjoint state `φ` per vertex: the solution of `(L φ)_i = −wE_i m_i u_i`, zero on
    /// the boundary. On a closed mesh it is defined up to an additive constant and returned
    /// with `Σ m φ = 0`.
    ///
    /// If the cones are an optimum of the canonical problem for some `λ` (unit cone costs, no
    /// bounds), then `φ_i = λ·sign(µ_i)` at every cone and `|φ_i| ≤ λ` elsewhere (after adding
    /// the right constant on a closed mesh).
    pub adjoint: Vec<f64>,
    /// The distortion `½ Σ wE m u²`.
    pub distortion: f64,
    /// The total cone angle `Σ |µ_i|`.
    pub total_cone_angle: f64,
    /// Closed meshes: `Σ µ_i − Σ Ω_i`, which Gauss–Bonnet requires to vanish. If it does not,
    /// no scale factor exists for the cones as given; the returned `u` solves
    /// `(L u)_i = µ_i − Ω_i − defect · m_i` instead (the defect spread over the surface in
    /// proportion to area). Zero on meshes with boundary. On a mesh with several components:
    /// the defect of largest magnitude among the closed components (each component has its own,
    /// spread over that component; see `component_defects`).
    pub gauss_bonnet_defect: f64,
    /// The Gauss–Bonnet defect of every connected component, in the order of
    /// `TriMesh::face_components` (zero for a component with boundary). One entry for a
    /// connected mesh.
    pub component_defects: Vec<f64>,
}

pub(crate) fn cones_to_measure(d: &Discrete, cones: &[Cone]) -> Result<Vec<f64>, ConeError> {
    let mut unknown_of = vec![u32::MAX; d.n_vertices];
    for (i, &v) in d.unknowns.iter().enumerate() {
        unknown_of[v as usize] = i as u32;
    }
    let mut mu = vec![0.0; d.n];
    let mut seen = vec![false; d.n];
    for (index, cone) in cones.iter().enumerate() {
        let issue = |reason| ConeError::InvalidCone {
            index,
            vertex: cone.vertex,
            angle: cone.angle,
            reason,
        };
        let v = cone.vertex as usize;
        if v >= d.n_vertices {
            return Err(issue(PrescribedConeIssue::VertexOutOfRange));
        }
        if !cone.angle.is_finite() {
            return Err(issue(PrescribedConeIssue::NonFiniteAngle));
        }
        let i = unknown_of[v];
        if i == u32::MAX {
            return Err(issue(PrescribedConeIssue::BoundaryVertex));
        }
        if seen[i as usize] {
            return Err(issue(PrescribedConeIssue::DuplicateVertex));
        }
        seen[i as usize] = true;
        mu[i as usize] = cone.angle;
    }
    Ok(mu)
}

/// Solves `A u = µ − Ω` for a cone measure on the unknowns and returns `u` (with the
/// distortion-minimizing constant on a closed mesh) and the Gauss–Bonnet defect.
///
/// [formula: discrete-yamabe-linearization]
/// [formula: closed-mesh-optimal-constant]
pub(crate) fn log_scale_of_measure(
    d: &Discrete,
    laplace: &LaplaceSolver,
    mu: &[f64],
    par: Parallelism<'_>,
) -> Result<(Vec<f64>, f64), ConeError> {
    let mut rhs: Vec<f64> = mu.iter().zip(&d.omega).map(|(m, o)| m - o).collect();
    let mut defect = 0.0;
    if d.closed {
        defect = mu.iter().sum::<f64>() - d.total_curvature;
        let total_mass: f64 = d.mass.iter().sum();
        for (r, m) in rhs.iter_mut().zip(&d.mass) {
            *r -= defect * m / total_mass;
        }
    }
    let mut u = laplace.solve(&rhs, par)?;
    if d.closed {
        let weight: f64 = d.w.iter().sum();
        let mean = d.w.iter().zip(&u).map(|(w, u)| w * u).sum::<f64>() / weight;
        for value in &mut u {
            *value -= mean;
        }
    }
    Ok((u, defect))
}

/// `u`, `φ` (per unknown), the distortion and the Gauss–Bonnet defect of a measure on the
/// unknowns of a connected problem.
fn evaluate_discrete(
    d: &Discrete,
    mu: &[f64],
    par: Parallelism<'_>,
) -> Result<(Vec<f64>, Vec<f64>, f64, f64), ConeError> {
    if d.n == 0 {
        return Ok((Vec::new(), Vec::new(), 0.0, 0.0));
    }
    let laplace = LaplaceSolver::new(d, par)?;
    let (u, defect) = log_scale_of_measure(d, &laplace, mu, par)?;
    let rhs: Vec<f64> = d.w.iter().zip(&u).map(|(w, u)| -w * u).collect();
    let mut phi = laplace.solve(&rhs, par)?;
    if d.closed {
        let total_mass: f64 = d.mass.iter().sum();
        let mean = d.mass.iter().zip(&phi).map(|(m, p)| m * p).sum::<f64>() / total_mass;
        for value in &mut phi {
            *value -= mean;
        }
    }
    let distortion = d.distortion(&u);
    Ok((u, phi, distortion, defect))
}

/// Computes the log scale factor `u` and the adjoint state `φ` of a given cone configuration,
/// for the mesh, geometry and distortion weights of `problem` (its `λ`, cone costs, bounds and
/// prescribed cones play no role).
///
/// Use it to evaluate cones from another source, or to check a configuration against the
/// optimality conditions.
///
/// The mesh may have several connected components: each is evaluated on its own (its own
/// additive constant of `u` and of `φ` if it is closed, its own Gauss–Bonnet defect), with the
/// masses of the whole mesh (divided by the total area), so that `distortion` and `adjoint` are
/// those of the whole problem.
///
/// # Errors
/// The errors of the problem's validation; [`ConeError::InvalidCone`] for a cone at a boundary
/// vertex, out of range, listed twice or with a non-finite angle; [`ConeError::Sparse`] if the
/// Laplacian cannot be factorized.
///
/// [formula: component-separation]
pub fn evaluate_cones(
    problem: &ConeProblem<'_>,
    cones: &[Cone],
    par: Parallelism<'_>,
) -> Result<ConeEvaluation, ConeError> {
    // Only the mesh, the geometry and the distortion weights matter here.
    let relaxed = ConeProblem {
        mesh: problem.mesh,
        geometry: problem.geometry,
        lambda: 0.0,
        distortion_weights: problem.distortion_weights.clone(),
        cone_costs: None,
        lower_bounds: None,
        upper_bounds: None,
        prescribed_cones: Vec::new(),
    };
    let whole = Instance::new(&relaxed, BoundaryTreatment::Fixed)?;
    if problem.mesh.n_components() == 1 {
        let d = Discrete::from_instance(&whole)?;
        let mu = cones_to_measure(&d, cones)?;
        let total_cone_angle = mu.iter().map(|m| m.abs()).sum();
        let (u, phi, distortion, defect) = evaluate_discrete(&d, &mu, par)?;
        if d.n == 0 {
            return Ok(ConeEvaluation {
                log_scale: vec![0.0; d.n_vertices],
                adjoint: vec![0.0; d.n_vertices],
                distortion: 0.0,
                total_cone_angle,
                gauss_bonnet_defect: 0.0,
                component_defects: vec![0.0],
            });
        }
        return Ok(ConeEvaluation {
            distortion,
            log_scale: d.to_vertices(&u),
            adjoint: d.to_vertices(&phi),
            total_cone_angle,
            gauss_bonnet_defect: defect,
            component_defects: vec![defect],
        });
    }

    // Several components: the cones are checked against the whole mesh (vertex numbers of the
    // caller), then every component is evaluated with its part of the measure.
    let nv = problem.mesh.n_vertices();
    let mut angle = vec![0.0; nv];
    let mut seen = vec![false; nv];
    for (index, cone) in cones.iter().enumerate() {
        let issue = |reason| ConeError::InvalidCone {
            index,
            vertex: cone.vertex,
            angle: cone.angle,
            reason,
        };
        let v = cone.vertex as usize;
        if v >= nv {
            return Err(issue(PrescribedConeIssue::VertexOutOfRange));
        }
        if !cone.angle.is_finite() {
            return Err(issue(PrescribedConeIssue::NonFiniteAngle));
        }
        if !whole.data.unknown[v] {
            return Err(issue(PrescribedConeIssue::BoundaryVertex));
        }
        if seen[v] {
            return Err(issue(PrescribedConeIssue::DuplicateVertex));
        }
        seen[v] = true;
        angle[v] = cone.angle;
    }
    let total_cone_angle = cones.iter().map(|c| c.angle.abs()).sum();
    let pieces = split_pieces(problem.mesh, problem.geometry)?;
    let mut log_scale = vec![0.0; nv];
    let mut adjoint = vec![0.0; nv];
    let mut distortion = 0.0;
    let mut component_defects = Vec::with_capacity(pieces.len());
    for (c, piece) in pieces.iter().enumerate() {
        let d = Discrete::from_instance(&whole.part(piece, c as u32))?;
        let to_parent = &piece.sub.vertex_to_parent;
        let mu: Vec<f64> = d
            .unknowns
            .iter()
            .map(|&i| angle[to_parent[i as usize] as usize])
            .collect();
        let (u, phi, part_distortion, defect) = evaluate_discrete(&d, &mu, par)?;
        for (k, &i) in d.unknowns.iter().enumerate() {
            let v = to_parent[i as usize] as usize;
            log_scale[v] = u[k];
            adjoint[v] = phi[k];
        }
        distortion += part_distortion;
        component_defects.push(defect);
    }
    let gauss_bonnet_defect = component_defects
        .iter()
        .copied()
        .fold(0.0f64, |a, b| if b.abs() > a.abs() { b } else { a });
    Ok(ConeEvaluation {
        log_scale,
        adjoint,
        distortion,
        total_cone_angle,
        gauss_bonnet_defect,
        component_defects,
    })
}
