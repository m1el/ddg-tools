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
//! The contract of this crate is `SPEC.md` next to its `Cargo.toml`: every discrete formula with
//! its `rc-` formula id, every place where the paper or the reference is wrong or silent
//! (decisions `D1`, `D2`, ...), the public API, the defaults, the verification and measurement
//! plans and the work packages. This file is a placeholder until work package WP1 lands.

#![warn(missing_docs)]
