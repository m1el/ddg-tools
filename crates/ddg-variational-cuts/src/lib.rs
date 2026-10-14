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
//! Two regions described by one level-set function `φ` (region 0 where `φ ≤ 0`, region 1 where
//! `φ > 0`).
//!
//! **The static evaluation of a cut** (SPEC §2, §3.1–§3.7, §6, §7). An [`Evaluator`] is built
//! once per mesh and maps `φ` to an [`Evaluation`]:
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
//! **The flow** (SPEC §3.8–§3.12, §4, §5). A [`Flow`] owns an evaluator and the state `φ`;
//! [`Flow::step`] takes one step and reports it, [`Flow::run`] steps until a [`StopRule`] says
//! so (a run whose energy rises ends with [`StopReason::EnergyRose`] and is set back to its
//! state of lowest energy; every [`RunReport`] says where that was). The problem
//! ([`CutProblem`]: energy, `α_L`, `η`, `w_r`, `β_r`, curvature) and the options
//! ([`FlowOptions`]: step size, symmetry classes, the choice of every stage) are arguments of
//! each step, so all of them can change between steps. [`init`] has helpers for the initial
//! cut and the visibility penalty of the paper, [`Cut::polylines`] returns the cut as ordered
//! polylines.
//!
//! The default stages differ from the paper's in three places, each of which was forced by a law
//! the flow must obey and is selectable back ([`FlowOptions::paper`]): the length term acts on
//! the cut instead of on the mesh ([`LengthStep`], SPEC D14), the rate is extended from the cut
//! outward instead of harmonically ([`Extension`], D15), and it is scaled by the mean slope of
//! the level-set function at the cut ([`SpeedScaling`], D16).
//!
//! The distortion term roughens a cut and the length term smooths it; every step reports the
//! rate `ω` that sets the scale of both ([`StepReport::roughening_rate`], SPEC D18), and
//! [`length_weight_for_ripples`] and [`StepSize::Adaptive`] turn it into a length weight and a
//! time step. `ω` is the 0.9-quantile of the local rates along the cut, weighted by length, and
//! the step rule slows down the cut points beyond it one by one, so that a spike of curvature
//! does not set the weight or the step of the whole cut ([`StepStatistics`], SPEC D22;
//! [`StepStatistics::Extremes`] is the rule with the largest rate and speed).
//!
//! **Several regions** (SPEC §8). [`RegionFlow`] is the same flow for `N` level-set functions,
//! one per region (negative inside; a vertex belongs to the region whose function is smallest).
//! The cut is then a network: [`Cut::from_regions`] and [`Evaluator::evaluate_regions`] add
//! triple points ([`TriplePoint`]) in the faces whose corners lie in three regions, and
//! [`Cut::curves`] returns the curves between them. [`init::regions_from_vertex_labels`],
//! [`init::regions_from_face_labels`] and [`init::normal_clusters`] make a start. With the two
//! functions `φ, −φ` everything reproduces the two-region code. [`FlowOptions::for_regions`]
//! selects the paper's stages for more than two regions; the measurements for and against are
//! in its documentation and in SPEC §8.5. [`RegionFlow::step_with_stages`] returns every stage
//! of a step ([`RegionStages`]), as [`Flow::step_with_stages`] does for two regions.
//!
//! ```
//! use ddg_core::generators::icosphere;
//! use ddg_core::{Geometry, TriMesh};
//! use ddg_variational_cuts::{
//!     CutProblem, Energy, FlowOptions, RegionFlow, StepSize, StopRule, init,
//! };
//!
//! # fn main() -> Result<(), Box<dyn std::error::Error>> {
//! let sphere = icosphere(2)?;
//! let mesh = TriMesh::new(sphere.positions.len(), sphere.faces.clone())?;
//! let geometry = Geometry::from_positions(&mesh, &sphere.positions)?;
//!
//! // Four clusters of the vertex normals as a start: one level-set function per region.
//! let labels = init::normal_clusters(&mesh, &sphere.positions, 4, 50)?;
//! let phis = init::regions_from_vertex_labels(&mesh, &geometry, &labels, 4)?;
//! let mut flow = RegionFlow::new(&mesh, &geometry, &phis)?;
//!
//! let problem = CutProblem::for_regions(Energy::Dirichlet, 3.0, 4);
//! let options = FlowOptions::for_regions(StepSize::default());
//! let stop = StopRule { max_steps: 10, ..StopRule::default() };
//! let report = flow.run(&problem, &options, &stop)?;
//! assert!(report.final_energy.total < report.steps[0].energy.total);
//!
//! // The network: curves between triple points, each between two regions.
//! let network = flow.evaluate(&problem, &options.eval)?;
//! assert!(network.report.n_triple_points >= 4);
//! let curves = network.cut.curves(&mesh, &sphere.positions)?;
//! assert!(curves.iter().all(|curve| curve.regions[0] < curve.regions[1]));
//! // The cut mesh has one piece per patch, and every face knows its region.
//! assert!(network.cut_mesh()?.mesh.n_components() >= 4);
//! # Ok(())
//! # }
//! ```
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
//! The flow, started from a geodesic disk:
//!
//! ```
//! use ddg_core::generators::icosphere;
//! use ddg_core::{Geometry, TriMesh};
//! use ddg_variational_cuts::{
//!     CutProblem, Energy, Flow, FlowOptions, StepSize, StopRule, init,
//!     length_weight_for_ripples,
//! };
//!
//! # fn main() -> Result<(), Box<dyn std::error::Error>> {
//! let sphere = icosphere(2)?;
//! let mesh = TriMesh::new(sphere.positions.len(), sphere.faces.clone())?;
//! let geometry = Geometry::from_positions(&mesh, &sphere.positions)?;
//!
//! // Region 0 is the disk of geodesic radius 1 around vertex 0 (lengths in the caller's unit).
//! let phi0 = init::geodesic_circle(&mesh, &geometry, 0, 1.0)?;
//! let mut flow = Flow::new(&mesh, &geometry, &phi0)?;
//!
//! // A length weight under which ripples of the cut shorter than 0.6 (a third of the equator
//! // of the unit-area sphere) decay.
//! let mut problem = CutProblem::new(Energy::Dirichlet, 0.0);
//! let omega = flow.roughening_rate(&problem, &Default::default())?;
//! problem.length_weight = length_weight_for_ripples(omega, 0.6);
//!
//! let options = FlowOptions::new(StepSize::default());
//! let stop = StopRule { max_steps: 20, ..StopRule::default() };
//! let report = flow.run(&problem, &options, &stop)?;
//! assert!(report.final_energy.total < report.steps[0].energy.total);
//!
//! // The state is the level-set function; the cut is one closed polyline.
//! let cut = flow.evaluate(&problem, &options.eval)?.cut;
//! assert_eq!(cut.polylines(&mesh, &sphere.positions)?.len(), 1);
//! assert_eq!(flow.level_set().len(), mesh.n_vertices());
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
mod flow;
pub mod init;
mod problem;
mod recovery;
mod regions;
mod system;

pub use ddg_core::sparse::{Parallelism, ThreadPool};
pub use ddg_core::topology::NONE;

pub use curvature::{face_curvature_sum_residual, face_curvatures};
pub use cut::{
    CurveEnd, Cut, CutCurve, CutPoint, CutPolyline, CutSegment, TriplePoint, edge_crossing,
    region_of, regions_of,
};
pub use cut_mesh::{
    CutMesh, CutTriangles, CutVertex, EdgeKind, aspect_ratio, barycentric_dot, subtriangle_area,
};
pub use error::CutError;
pub use evaluate::{
    EnergyTerms, EvalReport, Evaluation, Evaluator, adjoint_source_pointwise,
    adjoint_source_scaled, cut_speed, dirichlet_energy_term, dual_lengths, gradient_dirichlet,
    gradient_pointwise, gradient_scaled, normal_flux, smooth_slivers,
};
pub use flow::{
    Extension, Flow, FlowOptions, GradientNorm, LengthStep, PenaltyRule, ROUGHENING_QUANTILE,
    Redistance, RunReport, SpeedScaling, StepReport, StepSize, StepStages, StepStatistics,
    StepTimings, StopReason, StopRule, energy_criterion, explicit_step, length_weight_for_ripples,
    motion_criterion, symmetry_projection, weighted_quantile,
};
pub use problem::{
    CurvatureSource, CutProblem, Energy, EvalOptions, FluxRecovery, PAPER_FEATURE_SMOOTHING,
    PointwiseEnergy, RegionWeights, SUPPORTED_REGIONS, ScaleWeight, UnitAreaScale,
};
pub use recovery::flux_smoothing_stiffness;
pub use regions::{RegionFlow, RegionStages};
