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
//! # Contents
//!
//! * [`CurveNetwork`]: validated topology, components, edge quantities, geometric validation
//!   (`SPEC.md` §2); [`read_curve_obj`] and [`write_curve_obj`]; [`Exponents`],
//!   [`kernel`](mod@kernel): the tangent-point kernel and gradients (§3); the exact `O(|E|²)`
//!   [`energy`](fn@energy::energy), parallel and deterministic, with the near-pair quadrature
//!   of §3.4 ([`NearPairOptions`], [`SecondNeighbourRule`]; not in the paper, on by default).
//! * Constraints and their projection ([`constraints`], [`solve`], §5–6), the fractional Sobolev
//!   metric and Gram matrices ([`metric`]), subdivision ([`mod@subdivide`]).
//! * The gradient flow ([`step`]): [`Session`] and [`flow`], the line searches ([`LineSearch`],
//!   with the start rules of [`StartRule`]), the stopping rules, length schedules, L-BFGS and
//!   NCG for the soft problems ([`optimize`], [`soft_options`]) and the continuous collision
//!   detection of [`ccd`] that guards the isotopy class.
//! * The Barnes–Hut energy ([`BhTree`], in the session as [`BarnesHutEnergy`]), the block cluster
//!   tree and the hierarchical metric product ([`BlockClusterTree`], [`HierarchicalMetric`],
//!   [`HierarchicalCgSaddle`]).
//! * The multigrid saddle solver for curves and curve networks ([`mod@multigrid`]: [`MgOptions`],
//!   [`CurveHierarchy`], [`MultigridSaddle`], [`MultigridSolvers`], the near-contact correction
//!   of [`NearCorrection`]; §8.4, §17g), in the session as [`Solver::Multigrid`] with a dense
//!   fallback, and [`Backend::Auto`] with a measured crossover (§17h).
//! * Obstacles ([`Obstacle`]), penalties ([`SurfacePotential`], [`TotalLengthPenalty`],
//!   [`LengthDifference`], [`FieldAlignment`], [`Hencky`]) and the implicit surfaces and vector
//!   fields of [`surfaces`].
//! * A reader for the authors' scene format ([`scene`]).
//!
//! Measurements, including the paper's Fig. 20 and Sec. 7 comparisons, are in `SPEC.md` §17
//! (§17s holds the Sec. 7 results of record) and in the crate's README.

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
pub use energy::near_pair::SecondNeighbourRule;
pub use energy::{BhOptions, EnergyEval, EnergyReport, NearPairOptions, energy};
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
    DEFAULT_NEWTON_POLISH, DirectSolvers, EnergyModel, ExactEnergy, FlowOptions, FlowReport,
    FlowResult, GROWTH_BACKOFF_MAX, PenaltyTerm, RetargetReport, RetargetStatus, Session, Solver,
    SolverFactory, SolverFallback, SolverRequest, StepReport, Timings, flow, soft_options,
};
pub use subdivide::{SubdivisionMap, subdivide};
pub use surfaces::{
    CircularField, ConstantField, DoubleTorus, Metaballs, Plane, RoundedBox, SmoothUnion, Sphere,
    SphereTangentField, Torus, Union, VectorField,
};
