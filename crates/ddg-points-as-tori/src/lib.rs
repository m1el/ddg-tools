//! Points as Tori: fast pointwise signed distance for oriented point clouds.
//!
//! An implementation of
//!
//! > Nicole Feng, Ioannis Gkioulekas, Keenan Crane. *Points as Tori: Fast Pointwise Signed
//! > Distance for Point Clouds.* ACM Transactions on Graphics 45(4), Article 53, 2026.
//! > <https://doi.org/10.1145/3811385>
//!
//! For every point of an oriented point cloud a small transformer predicts the six coefficients
//! of a local quadratic height function; a torus is fitted to its shift, principal curvatures and
//! directions, and the signed distance at any query point is a distance-weighted blend of the
//! tori's closed-form signed distance functions.
//!
//! The authors' reference implementation (`points-as-tori`, MIT, "Copyright (c) 2026 Nicole
//! Feng") may be read and ported; its licence notice is kept in `LICENSES/points-as-tori-MIT.txt`
//! and covers the converted pretrained weights as well.
//!
//! The contract of this crate is `SPEC.md` next to its `Cargo.toml` (formula ids `pat-…`, the
//! decisions, the public API, the verification and measurement plans and the work packages);
//! `NOTES-REFERENCE.md` records what the reference does where the paper is silent or differs.
//!
//! This is the skeleton created with the contract; the modules are added by the work packages.

#![warn(missing_docs)]
