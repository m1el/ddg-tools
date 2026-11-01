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
//! * Later work packages (`SPEC.md` §15) fill: [`metric`], [`constraints`], [`solve`],
//!   [`spatial`], [`energy::barnes_hut`], [`remesh`], [`step`], [`penalties`], [`obstacle`],
//!   [`implicit`], [`ccd`], [`optimize`], [`generate`], [`scene`].

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

pub use energy::{BhOptions, Energy, EnergyEval, EnergyReport, energy};
pub use error::RsError;
pub use geometry::FaceGeometry;
pub use io::{obj_to_string, parse_obj, read_obj, write_obj};
pub use kernel::{Exponents, KernelEval, kernel, kernel_with_gradient};
pub use penalties::Penalty;
pub use surface::{Component, NONE, Orientation, Surface};
