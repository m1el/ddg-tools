//! Repulsive surfaces: the tangent-point energy of triangle meshes and its fractional Sobolev
//! gradient flow.
//!
//! An implementation of
//!
//! > Chris Yu, Caleb Brakensiek, Henrik Schumacher, Keenan Crane. *Repulsive Surfaces.* ACM
//! > Transactions on Graphics 40(6), Article 268, 2021.
//! > <https://doi.org/10.1145/3478513.3480521>
//!
//! The authors' reference implementation (`repulsive-surfaces`, MIT, "Copyright (c) 2019 Nick
//! Sharp") may be read and ported; its licence notice is kept in
//! `LICENSES/repulsive-surfaces-MIT.txt`. The generic machinery (BVH, block cluster tree,
//! L-BFGS, implicit surfaces, CCD primitives) comes from `ddg-repulsive-curves` (Yu,
//! Schumacher, Crane, *Repulsive Curves*, ACM TOG 40(2), 2021), a dependency.
//!
//! The contract of this crate is `SPEC.md` next to its `Cargo.toml`: every discrete formula
//! with its `rs-` formula id, every place where the paper or the reference is wrong or silent
//! (decisions `D1`–`D31`), the public API, the verification and measurement plans and the work
//! packages (§15).
//!
//! # Contents
//!
//! * [`Surface`]: a validated manifold triangle mesh, orientation optional (Klein bottles
//!   welcome), with edges, fans, components, boundary loops and the [`Orientation`] class of
//!   every component; [`read_obj`], [`write_obj`].
//! * [`geometry`]: face areas, normals, barycenters and their derivatives, vertex areas, the
//!   cotan stiffness `W` and the lumped mass `M`.
//! * [`Exponents`] and the [`kernel`](mod@kernel) `K(X, N; Y) = |N·(X − Y)|^α |X − Y|^{−β}`.
//! * The exact energy and its differential ([`energy`](fn@energy::energy), [`energy::exact`]),
//!   parallel and bitwise independent of the thread count.
//! * The flow ([`step`]): [`FlowOptions`], [`Session`] (step by step, cancellable, the
//!   surface changes with remeshing), [`flow`], [`run_scene`] for the reference's scene files.
//! * The parts it is made of: [`metric`] (the Sobolev and comparison metrics), [`constraints`],
//!   [`solve`] (dense, hierarchical GMRES, Schur complement, projection), [`spatial`] and
//!   [`energy::barnes_hut`], [`remesh`], [`penalties`], [`obstacle`], [`implicit`], [`ccd`],
//!   [`generate`] (test meshes), [`scene`]; [`optimize`] (the comparison optimizers, WP9).

#![warn(missing_docs)]
#![allow(clippy::needless_range_loop)]

pub mod ccd;
pub mod constraints;
pub mod energy;
pub mod error;
pub mod generate;
pub mod geometry;
pub mod implicit;
pub mod io;
pub mod kernel;
pub mod metric;
pub mod obstacle;
pub mod optimize;
pub mod penalties;
pub mod remesh;
pub mod scene;
pub mod solve;
pub mod spatial;
pub mod step;
pub mod surface;

pub use ccd::CcdOptions;
pub use constraints::{BarycenterMode, Constraint};
pub use energy::{BhOptions, Energy, EnergyEval, EnergyEvaluation, EnergyReport, energy};
pub use error::RsError;
pub use geometry::FaceGeometry;
pub use implicit::{Cylinder, ImplicitSurface};
pub use io::{obj_to_string, parse_obj, read_obj, write_obj};
pub use kernel::{Exponents, KernelEval, kernel, kernel_with_gradient};
pub use metric::dense::{DenseMetric, dense_metric};
pub use metric::{FractionalForm, L2Kind, Metric, SobolevWeights};
pub use obstacle::{ObstacleKind, StaticSet};
pub use penalties::{Penalty, Quadrature};
pub use remesh::{RemeshOptions, RemeshOrder, RemeshReport};
pub use scene::SceneSpec;
pub use solve::iterative::{Backend, HierarchicalOptions};
pub use solve::{GmresOptions, GmresSide, LinearReport, NewtonSteps, ProjectionStatus};
pub use step::{
    ArmijoRule, BcqnVariant, Control, FlowOptions, FlowReport, FlowResult, LbfgsVariant, LineSearch, Method,
    PinSchedule, RunReport, RunSceneOptions, Schedule, ScheduleTiming, Session, StartRule,
    Status, StepReport, StepStatus, StopRule, Timings, Trial, TrialOutcome, flow, run_scene,
    scene_options,
};
pub use surface::{Component, NONE, Orientation, Surface};
