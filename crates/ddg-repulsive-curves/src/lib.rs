//! Repulsive curves: the tangent-point energy of curve networks and its fractional Sobolev
//! gradient flow.
//!
//! An implementation of
//!
//! > Chris Yu, Henrik Schumacher, Keenan Crane. *Repulsive Curves.* ACM Transactions on
//! > Graphics 40(2), Article 10, 2021. <https://doi.org/10.1145/3439429>
//!
//! with the paper's supplemental material. The authors' reference implementation
//! (`repulsive-curves`, MIT, Copyright (c) 2019 Christopher Yu) may be read and ported; its
//! licence notice is kept in `LICENSES/repulsive-curves-MIT.txt`.
//!
//! The contract of this crate is `SPEC.md` next to its `Cargo.toml`: every discrete formula
//! with its `rc-` formula id, every place where the paper or the reference is wrong or silent
//! (decisions `D1`, `D2`, ...), the public API, the defaults, the verification and measurement
//! plans and the work packages.
//!
//! # Implemented so far (work packages WP1–WP5 and the line-search follow-up)
//!
//! * WP1: [`CurveNetwork`]: validated topology, components, edge quantities, geometric
//!   validation (`SPEC.md` §2); [`read_curve_obj`] and [`write_curve_obj`]; [`Exponents`],
//!   [`kernel`](mod@kernel): the tangent-point kernel and gradients (§3); the exact `O(|E|²)`
//!   [`energy`](fn@energy::energy), parallel and deterministic.
//! * WP2: constraints and their projection ([`constraints`], [`solve`], §5–6), the fractional
//!   Sobolev metric and Gram matrices ([`metric`]), subdivision ([`mod@subdivide`]).
//! * WP3: the gradient flow ([`step`]): [`Session`] and [`flow`], the line searches
//!   ([`LineSearch`], with the start rules of [`StartRule`]), the stopping rules, length
//!   schedules, L-BFGS and NCG ([`optimize`]) and the continuous collision detection of
//!   [`ccd`] that guards the isotopy class.
//! * WP4: the Barnes–Hut energy ([`BhTree`]), wired into the session as [`BarnesHutEnergy`]
//!   (the tree is built once per step and refitted in the line search).
//! * WP5: the block cluster tree and the hierarchical metric product ([`BlockClusterTree`],
//!   [`HierarchicalMetric`], [`HierarchicalCgSaddle`]).
//! * WP7: obstacles ([`Obstacle`]), penalties ([`SurfacePotential`], [`TotalLengthPenalty`],
//!   [`LengthDifference`], [`FieldAlignment`]) and the implicit surfaces and vector fields of
//!   [`surfaces`].
//! * Follow-up: the grow-previous start of the guarded search (`SPEC.md` §17f).
//! * WP6 and WP6b: the multigrid saddle solver for curves and curve networks ([`mod@multigrid`]:
//!   [`MgOptions`], [`CurveHierarchy`], [`MultigridSaddle`], [`MultigridSolvers`], the
//!   near-contact correction of [`NearCorrection`]; `SPEC.md` §8.4, §17g).
//! * WP3c: multigrid in the session ([`Solver::Multigrid`]), the dense fallback, and
//!   [`Backend::Auto`] with a measured crossover (`SPEC.md` §17h).
//!
//! The remaining modules are filled by the later work packages of `SPEC.md` §15.

#![warn(missing_docs)]
#![allow(clippy::needless_range_loop)]

pub mod ccd;
pub mod constraints;
pub mod energy;
pub mod error;
pub mod io;
pub mod kernel;
pub mod metric;
pub mod multigrid;
pub mod network;
pub mod obstacle;
pub mod optimize;
pub mod penalties;
pub mod scene;
pub mod solve;
pub mod spatial;
pub mod step;
pub mod subdivide;
pub mod surfaces;

pub use ccd::{
    CcdOptions, CcdReport, CollisionObstacle, Contact, ContactPrimitive, PreparedObstacle,
    first_contact, path_contact,
};
pub use constraints::{
    Constraint, ConstraintSet, ImplicitSurface, RankReport, TangentSite, default_constraints,
};
pub use energy::barnes_hut::{BhStats, BhTree};
pub use energy::{BhOptions, EnergyEval, EnergyReport, energy};
pub use error::RcError;
pub use io::{parse_curve_obj, read_curve_obj, write_curve_obj};
pub use kernel::{Exponents, KernelEval, kernel, kernel_with_gradient};
pub use metric::hierarchical::{CgOptions, HierarchicalCgSaddle, HierarchicalMetric};
pub use metric::{DenseGram, DenseMatrix, DofLayout, GramQuadrature, Metric, gram};
pub use multigrid::{
    Coarsening, CurveHierarchy, CycleOptions, CycleReport, Hierarchy, HierarchyInput, Level,
    MgOptions, MultigridSaddle, MultigridSolvers, NearCorrection, PairCorrection, ProlongationRow,
};
pub use network::{Axis, Component, CurveNetwork, EdgeGeometry, check_planar};
pub use obstacle::{Obstacle, PotentialValue, closest_point_on_triangle};
pub use optimize::{Hencky, Method};
pub use penalties::{
    FieldAlignment, LengthDifference, SurfacePotential, TotalLengthPenalty, calibrate_hencky_weight,
};
pub use solve::{
    DenseSaddle, LinearReport, ProjectedGradient, Projection, ProjectionOptions, ProjectionStatus,
    Projector, SaddleSolution, SaddleSolver, SparseSaddle, Tolerance, project_constraints,
    projected_gradient,
};
pub use spatial::bct::{BctOptions, BctStats, BlockClusterTree, ClusterKernel};
pub use spatial::{Aabb, Bvh, Element, Node};
pub use step::line_search::{LineSearch, StartRule, Trial, TrialOutcome};
pub use step::schedule::{LengthSchedule, ScheduleProgress};
pub use step::stop::{GradNorm, Status, StepStatus, StopRule};
pub use step::{
    Backend, BarnesHutEnergy, Control, DEFAULT_AUTO_DENSE_MAX_VERTICES,
    DEFAULT_DENSE_FALLBACK_MAX_VERTICES, DEFAULT_DENSE_THREADS, DEFAULT_MIN_EDGE_LENGTH_RATIO,
    DirectSolvers, EnergyModel, ExactEnergy, FlowOptions, FlowReport, FlowResult, PenaltyTerm,
    RetargetReport, RetargetStatus, Session, Solver, SolverFactory, SolverFallback, SolverRequest,
    StepReport, Timings, flow, soft_options,
};
pub use subdivide::{SubdivisionMap, subdivide};
pub use surfaces::{
    CircularField, ConstantField, DoubleTorus, Metaballs, Plane, RoundedBox, SmoothUnion, Sphere,
    SphereTangentField, Torus, Union, VectorField,
};
