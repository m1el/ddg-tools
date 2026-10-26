//! Optimal cone singularities for conformal flattening.
//!
//! An implementation of
//!
//! > Yousuf Soliman, Dejan Slepčev, Keenan Crane.
//! > *Optimal Cone Singularities for Conformal Flattening.*
//! > ACM Transactions on Graphics 37(4), 2018. <https://doi.org/10.1145/3197517.3201367>
//!
//! Given a triangle mesh, it finds where to put cone singularities, and with which angles, so
//! that the mesh can be flattened conformally with the least area distortion for the total
//! cone angle spent. The result is certified: every solve checks its own optimality
//! conditions and duality gap.
//!
//! # Example
//!
//! ```
//! use ddg_core::generators::bumpy_grid_disk;
//! use ddg_core::{Geometry, TriMesh};
//! use ddg_optimal_cones::{ConeProblem, SolveOptions, Status, solve};
//!
//! # fn main() -> Result<(), Box<dyn std::error::Error>> {
//! // A disk with bumps; any connected manifold triangle mesh will do.
//! let shape = bumpy_grid_disk(12, 12, 0.6)?;
//! let mesh = TriMesh::new(shape.positions.len(), shape.faces.clone())?;
//! let geometry = Geometry::from_positions(&mesh, &shape.positions)?;
//!
//! // λ trades distortion against total cone angle: smaller λ, more cones.
//! let problem = ConeProblem::new(&mesh, &geometry, 0.002);
//! let solution = solve(&problem, &SolveOptions::default())?;
//!
//! // A solve that does not reach a certified optimum is not an error: check the status.
//! assert_eq!(solution.report.status, Status::Converged);
//! assert!(solution.report.certificate.holds(1e-9));
//! assert!(!solution.cones.is_empty());
//! for cone in &solution.cones {
//!     // `angle` is the curvature of the flattened metric at the vertex: positive for a tip.
//!     assert!(!mesh.is_boundary_vertex(cone.vertex as usize) && cone.angle != 0.0);
//! }
//! # Ok(())
//! # }
//! ```
//!
//! Steering goes through the fields of [`ConeProblem`]:
//!
//! ```
//! # use ddg_core::generators::bumpy_grid_disk;
//! # use ddg_core::{Geometry, TriMesh};
//! use ddg_optimal_cones::{Cone, ConeProblem, SolveOptions, solve};
//! # fn main() -> Result<(), Box<dyn std::error::Error>> {
//! # let shape = bumpy_grid_disk(12, 12, 0.6)?;
//! # let mesh = TriMesh::new(shape.positions.len(), shape.faces.clone())?;
//! # let geometry = Geometry::from_positions(&mesh, &shape.positions)?;
//! let nv = mesh.n_vertices();
//! let interior = mesh.interior_vertices();
//!
//! let mut cone_costs = vec![1.0; nv];
//! cone_costs[interior[0] as usize] = 0.0; // a free cone site: a cone here costs nothing
//! for &v in &interior[10..40] {
//!     cone_costs[v as usize] = f64::INFINITY; // no cone in this region
//! }
//! let mut distortion_weights = vec![1.0; nv];
//! for &v in &interior[40..60] {
//!     distortion_weights[v as usize] = 4.0; // distortion matters more here
//! }
//! let problem = ConeProblem {
//!     cone_costs: Some(cone_costs),
//!     distortion_weights: Some(distortion_weights),
//!     lower_bounds: Some(vec![-std::f64::consts::FRAC_PI_2; nv]), // bound the angles
//!     upper_bounds: Some(vec![std::f64::consts::FRAC_PI_2; nv]),
//!     // pin a cone: this vertex gets exactly this angle, at no cost
//!     prescribed_cones: vec![Cone { vertex: interior[5], angle: 0.5 }],
//!     ..ConeProblem::new(&mesh, &geometry, 0.002)
//! };
//! let solution = solve(&problem, &SolveOptions::default())?;
//! assert_eq!(solution.measure[interior[5] as usize], 0.5);
//! assert!(interior[10..40].iter().all(|&v| solution.measure[v as usize] == 0.0));
//! assert!(solution.measure.iter().all(|m| m.abs() <= std::f64::consts::FRAC_PI_2));
//! # Ok(())
//! # }
//! ```
//!
//! # The problem
//!
//! With `L` the positive semidefinite cotan Laplacian, `Ω` the angle defects, `m` the lumped
//! vertex masses divided by the total area, and `u` the log conformal scale factor
//! (flattened edge length `≈ e^{(u_i+u_j)/2} ℓ_ij`), [`solve`] finds the cone angles `µ` that
//!
//! ```text
//! minimize     ½ Σ_i wE_i m_i u_i²  +  λ Σ_i wR_i |µ_i|
//! subject to   (L u)_i = µ_i − Ω_i  (i interior),   u_i = 0  (i boundary),   lo_i ≤ µ_i ≤ hi_i
//! ```
//!
//! * `λ ≥ 0` ([`ConeProblem::lambda`]) is dimensionless because the masses sum to one. It is
//!   the paper's `λ` divided by the surface area if the paper's masses are left as they are.
//! * `wE_i ≥ 0` ([`ConeProblem::distortion_weights`]) weights the distortion.
//! * `wR_i ∈ [0, +∞]` ([`ConeProblem::cone_costs`]) weights the cone cost: `0` is a free cone
//!   site, `+∞` forbids a cone.
//! * `lo_i ≤ 0 ≤ hi_i` ([`ConeProblem::lower_bounds`], [`ConeProblem::upper_bounds`]) bound
//!   the angle.
//! * [`ConeProblem::prescribed_cones`] fix the angle of a vertex at no cost.
//! * On a closed mesh `u` is defined up to a constant, which is optimized too, and the angles
//!   sum to `2πχ`; constraints that do not allow that are reported as
//!   [`ConeError::Infeasible`].
//!
//! * By default `u = 0` on the boundary. [`solve_with_boundary`] with a [`BoundaryTreatment`]
//!   optimizes the boundary values too (paper Sec. 8.3): `Free`, `Penalized`, `Convex`; see the
//!   crate's `README.md`, section "The boundary".
//!
//! A cone angle is the curvature of the flattened metric at the vertex,
//! `µ_i = 2π − (angle sum around i after flattening)`: positive for a cone tip.
//!
//! # What comes back
//!
//! [`Solution`] holds the optimum (`measure`, `log_scale_optimal`, `adjoint`, `objective`), the
//! cones extracted from it (`cones`, `log_scale`, `distortion`; they differ from the optimum
//! only if a threshold or the optional merging of cone clusters changed something, which
//! [`ExtractionReport`] lists), and a [`SolveReport`] with the status, the stages, iteration
//! and factorization counts, timings, and the [`Certificate`] of the returned optimum.
//!
//! [`certify`] computes that certificate for any triple; [`evaluate_cones`] computes the scale
//! factor and adjoint state of a given cone configuration; [`reoptimize_cones`] finds the
//! least-distortion angles for given cone positions (by a direct solve of the linear optimality
//! system, or the reduced system for few cones); [`feature_size_weights`] builds the
//! curvature-adapted weights of the paper's Sec. 8.1. [`solve_with_feedback`] (and the state
//! machine [`Feedback`] behind it) lowers the distortion weight where the solver's linearized
//! scale factor is far from the one a flattening realizes, e.g. at needle tips (crate README,
//! "The feedback loop").
//!
//! The mesh may have several connected components. Every function then works on the canonical
//! problem of the whole mesh (masses divided by the total area, one constant of `u` per closed
//! component); the problem separates, each component is solved as its part of it, and
//! [`SolveReport::components`] reports each component (see the crate's `README.md`, section
//! "Disconnected meshes").
//!
//! # Method, deviations, measurements
//!
//! The method is the paper's two-stage semismooth Newton method (regularized problems, then
//! the exact optimality system), with safeguards, a treatment of closed meshes and of angle
//! bounds that the paper leaves open, and a certificate. The crate's `README.md` lists every
//! deviation from the paper with its reason, what was found about the authors' published
//! cones, how `λ` relates to theirs, and measured running times.
//!
//! Functions that implement a formula carry a `[formula: <id>]` line; the ids are those of
//! `lean/proofs`, where the formula is proved or cited.

#![warn(missing_docs)]

mod certificate;
mod error;
mod evaluate;
mod extract;
mod factor_once;
mod feedback;
mod linsys;
mod newton;
mod problem;
mod reoptimize;
mod site;
mod solve;
mod weights;

pub use ddg_core::io::Cone;
pub use ddg_core::sparse::{Parallelism, ThreadPool};

pub use certificate::{Certificate, certify, certify_with_boundary};
pub use error::{ConeError, PrescribedConeIssue};
pub use evaluate::{ConeEvaluation, evaluate_cones};
pub use extract::{ExtractionReport, MergedCluster, UnmergedCluster, UnmergedReason};
pub use factor_once::{
    DenseStats, FactorOnceOptions, FactorOnceReport, FactorOnceTimings, FallbackReason, MethodUsed,
    SolveMethod, Unsupported,
};
pub use feedback::{
    Feedback, FeedbackError, FeedbackOptions, FeedbackReport, FeedbackRound, FeedbackSolution,
    FeedbackStop, ScaleMismatch, every_component_has_boundary, scale_mismatch, solve_with_feedback,
};
pub use linsys::{LinearSolver, LinearSolverReport};
pub use newton::StageStop;
pub use problem::{BoundaryTreatment, ConeProblem};
pub use reoptimize::{
    Reoptimization, ReoptimizeComponentReport, ReoptimizeMethod, ReoptimizeMethodUsed,
    ReoptimizeOptions, ReoptimizeReport, reoptimization_problem, reoptimize_cones,
};
pub use site::soft_threshold;
pub use solve::{
    ComponentReport, MethodReport, RegularizationMetric, Solution, SolveOptions, SolveReport,
    StageKind, StageReport, Status, Timings, factor_once_support, solve, solve_with_boundary,
};
pub use weights::{FeatureSizeConvention, FeatureSizeWeights, feature_size_weights};
