//! Extrinsic baselines for the comparisons of the paper (Sec. 8.4; SPEC §10.4, S23).
//!
//! * [`qem`]: Garland and Heckbert's quadric error metrics edge collapse, with a
//!   closest-point-projection map of the fine vertices onto the coarse mesh.
//!
//! The successive self-parameterization of Liu et al. 2021 is a later work package.

pub mod qem;

pub use qem::{
    ProjectedPoint, QemError, QemOptions, QemReport, QemResult, Quadric, optimal_point,
    plane_quadric, qem_simplify,
};
