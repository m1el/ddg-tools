//! Optimal cone singularities for conformal flattening.
//!
//! (crate documentation is completed below)

#![warn(missing_docs)]

mod certificate;
mod error;
mod evaluate;
mod extract;
mod linsys;
mod newton;
mod problem;
mod site;
mod solve;
mod weights;

pub use ddg_core::io::Cone;
pub use ddg_core::sparse::{Parallelism, ThreadPool};

pub use certificate::{Certificate, certify};
pub use error::{ConeError, PrescribedConeIssue};
pub use evaluate::{ConeEvaluation, evaluate_cones};
pub use extract::{ExtractionReport, MergedCluster, UnmergedCluster, UnmergedReason};
pub use linsys::{LinearSolver, LinearSolverReport};
pub use newton::StageStop;
pub use problem::ConeProblem;
pub use site::soft_threshold;
pub use solve::{
    RegularizationMetric, Solution, SolveOptions, SolveReport, StageKind, StageReport, Status,
    Timings, solve,
};
pub use weights::{FeatureSizeConvention, FeatureSizeWeights, feature_size_weights};
