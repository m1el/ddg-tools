//! Extrinsic baselines for the comparisons of the paper (Sec. 8.4; SPEC §10.4, S23).
//!
//! * [`qem`]: Garland and Heckbert's quadric error metrics edge collapse, with a
//!   closest-point-projection map of the fine vertices onto the coarse mesh.
//!
//! * [`QemOptions::successive_map`]: the successive self-parameterization of Liu et al. 2021 built
//!   during the collapses (module `successive`).

pub mod qem;
mod successive;

pub use successive::SuccessiveStats;

pub use qem::{
    ProjectedPoint, QemError, QemOptions, QemReport, QemResult, Quadric, optimal_point,
    optimal_point_on_segment, plane_quadric, qem_simplify,
};
