//! Variational surface cutting.
//!
//! A clean-room implementation, from the text of the paper only, of
//!
//! > Nicholas Sharp, Keenan Crane. *Variational Surface Cutting.*
//! > ACM Transactions on Graphics 37(4), 2018. <https://doi.org/10.1145/3197517.3201356>
//!
//! The discrete method, the sign conventions and every difference from the paper as printed are
//! in `SPEC.md` next to this crate's `Cargo.toml` (the differences in its §9); measurements are
//! in the crate's `README.md`.
//!
//! # What is here
//!
//! The *static evaluation* of a cut (SPEC §2, §3.1–§3.7, §6, §7), for two regions described by
//! one level-set function `φ` (region 0 where `φ ≤ 0`, region 1 where `φ > 0`). An
//! [`Evaluator`] is built once per mesh and maps `φ` to an [`Evaluation`]:
//!
//! 1. the cut `γ = {φ = 0}` as points on edges and segments across faces ([`Cut`]);
//! 2. the cut mesh, in which every crossed face is split and the two sides of the cut are
//!    separated ([`CutTriangles`], and as a mesh with geometry [`CutMesh`]);
//! 3. in every region the log conformal scale factor `u` of the flattening with `u = 0` on the
//!    cut and on the boundary of the mesh (the Yamabe equation `Δu = −K`), its normal derivative
//!    on the boundary, and for energies other than the Dirichlet energy the adjoint state;
//! 4. the gradient density `G_r` of each region's energy on the cut and the speed
//!    `s = G_0 − G_1` at which the level-set function must change on the cut to descend
//!    (optionally with the sliver smoothing of the paper);
//! 5. the energy
//!    `E = Σ_r [ w_r J_r + β_r Area_r ] + ∫_γ (α_L + η) ds` with all its terms.
//!
//! The flow that moves the cut (harmonic extension of `s`, the implicit length step,
//! redistancing; SPEC §3.8–§3.12, §4) builds on this and is not part of this module set yet.
//!
//! # Example
//!
//! ```
//! use ddg_core::generators::icosphere;
//! use ddg_core::{Geometry, TriMesh};
//! use ddg_variational_cuts::{CutProblem, Energy, EvalOptions, Evaluator};
//!
//! # fn main() -> Result<(), Box<dyn std::error::Error>> {
//! let sphere = icosphere(3)?;
//! let mesh = TriMesh::new(sphere.positions.len(), sphere.faces.clone())?;
//! let geometry = Geometry::from_positions(&mesh, &sphere.positions)?;
//!
//! // One evaluator per mesh: it owns the symbolic factorization.
//! let mut evaluator = Evaluator::new(&mesh, &geometry)?;
//!
//! // A cut along the circle of latitude z = 0.3; region 0 is the cap above it.
//! let phi: Vec<f64> = sphere.positions.iter().map(|p| 0.3 - p[2]).collect();
//! let problem = CutProblem::new(Energy::Dirichlet, 0.1);
//! let evaluation = evaluator.evaluate(&phi, &problem, &EvalOptions::default())?;
//!
//! // The surface is scaled to unit area; the two regions share it.
//! let area = &evaluation.energy.area;
//! assert!((area[0] + area[1] - 1.0).abs() < 1e-12);
//! // The smaller cap is flatter after cutting, so it has the smaller gradient density, and the
//! // speed G_0 − G_1 is negative: the level-set function decreases, the small cap grows.
//! assert!(evaluation.speed.iter().all(|&s| s < 0.0));
//! // What a flattening consumes: the cut mesh with geometry and maps to the input.
//! let cut_mesh = evaluation.cut_mesh()?;
//! assert_eq!(cut_mesh.mesh.n_components(), 2);
//! # Ok(())
//! # }
//! ```
//!
//! # Conventions
//!
//! * **Unit area.** The method runs on the surface scaled to unit total area
//!   ([`UnitAreaScale`]). The parameters of a [`CutProblem`] and the numbers of an
//!   [`Evaluation`] refer to the scaled surface; [`CutMesh::geometry`] is in the caller's unit.
//! * **Signs.** The cotan Laplacian is positive semidefinite, `u` is the true log scale factor
//!   (negative inside a cap of a sphere, zero on its rim), normals point out of a region. A
//!   region whose own gradient density is large shrinks (SPEC D2).
//! * **Regions.** Region ids are `u32`. The cut mesh gives every cut point one vertex per
//!   side, so each of its vertices lies in exactly one region and per-region quantities are
//!   plain per-vertex arrays.
//! * **Normal derivative on the cut.** The paper divides the integrated normal derivative of
//!   each cut vertex by its dual length. That quotient is wrong by orders of magnitude where the
//!   cut passes close to a mesh vertex, so by default the fluxes are smoothed along the cut
//!   first ([`FluxRecovery`], SPEC §3.5 and D12). [`EvalOptions::paper`] selects the procedure
//!   of the paper.
//! * Bad input never panics: see [`CutError`]. Nothing is clamped or left out silently; what
//!   was moved or excluded is counted in the [`EvalReport`].
//!
//! Functions that implement a formula carry a `[formula: <id>]` line with the ids of SPEC §10.

#![warn(missing_docs)]

mod curvature;
mod cut;
mod cut_mesh;
mod error;
mod evaluate;
mod problem;
mod recovery;
mod system;

pub use ddg_core::sparse::{Parallelism, ThreadPool};
pub use ddg_core::topology::NONE;

pub use curvature::{face_curvature_sum_residual, face_curvatures};
pub use cut::{Cut, CutPoint, CutSegment, edge_crossing, region_of};
pub use cut_mesh::{
    CutMesh, CutTriangles, CutVertex, EdgeKind, aspect_ratio, barycentric_dot, subtriangle_area,
};
pub use error::CutError;
pub use evaluate::{
    EnergyTerms, EvalReport, Evaluation, Evaluator, adjoint_source_pointwise,
    adjoint_source_scaled, cut_speed, dirichlet_energy_term, dual_lengths, gradient_dirichlet,
    gradient_pointwise, gradient_scaled, normal_flux, smooth_slivers,
};
pub use problem::{
    CurvatureSource, CutProblem, Energy, EvalOptions, FluxRecovery, PAPER_FEATURE_SMOOTHING,
    PointwiseEnergy, RegionWeights, SUPPORTED_REGIONS, ScaleWeight, UnitAreaScale,
};
pub use recovery::flux_smoothing_stiffness;
