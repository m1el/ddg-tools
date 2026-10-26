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
//! # Implemented so far (work package WP1)
//!
//! * [`CurveNetwork`]: validated topology, components, edge quantities, geometric validation
//!   (`SPEC.md` §2); [`read_curve_obj`] and [`write_curve_obj`] for polyline OBJ files.
//! * [`Exponents`], [`kernel`](mod@kernel): the tangent-point kernel and its gradients (§3.1–3.2).
//! * [`energy`](fn@energy::energy): the exact `O(|E|²)` energy and differential, parallel and
//!   deterministic (§3).
//!
//! The other modules are stubs filled by the later work packages of `SPEC.md` §15.

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

pub use energy::{BhOptions, EnergyEval, EnergyReport, energy};
pub use error::RcError;
pub use io::{parse_curve_obj, read_curve_obj, write_curve_obj};
pub use kernel::{Exponents, KernelEval, kernel, kernel_with_gradient};
pub use network::{Axis, Component, CurveNetwork, EdgeGeometry, check_planar};
