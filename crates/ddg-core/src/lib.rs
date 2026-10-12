//! Foundation crate of the `ddg-tools` workspace: everything the discrete differential geometry
//! algorithms share.
//!
//! * [`topology`] — a validated, oriented, manifold triangle mesh ([`TriMesh`]) with halfedge
//!   connectivity, boundary loops, components, and explicit repair helpers.
//! * [`triangulate`] — polygon triangulation that matches Blender's, with a map from triangle
//!   corners back to the original polygon corners (cargo feature `blender-tessellation`).
//! * [`geometry`] — intrinsic triangle geometry (corner angles, areas, cotangents) computed from
//!   edge lengths ([`Geometry`]).
//! * [`operators`] — the cotan Laplacian, lumped mass, angle defects and boundary geodesic
//!   curvature, following the conventions in the repository's `CLAUDE.md`.
//! * [`sparse`] — a small sparse matrix type and direct solvers with reusable symbolic analysis.
//! * [`io`] — Wavefront OBJ reading and writing, and the cone files of the optimal-cones dataset.
//! * [`generators`] — small procedural test meshes.
//!
//! # References
//!
//! The discrete operators are those of Sec. 5.1 of
//! Yousuf Soliman, Dejan Slepčev, Keenan Crane, *Optimal Cone Singularities for Conformal
//! Flattening*, ACM Transactions on Graphics 37(4), 2018,
//! <https://doi.org/10.1145/3197517.3201367>. The cotan operator goes back to MacNeal (1949).
//! The numerically careful triangle area and angle formulas are from W. Kahan, *Miscalculating
//! Area and Angles of a Needle-like Triangle* (2014).
//!
//! # Conventions
//!
//! `f64` everywhere, stored indices are `u32`, ids at function boundaries are `usize`. Invalid
//! input is reported through error enums and never panics. Nothing is normalized, clamped or
//! repaired silently. Functions that implement a mathematical formula carry a
//! `[formula: <id>]` line in their documentation.

#![warn(missing_docs)]

pub mod generators;
pub mod geometry;
pub mod io;
pub mod operators;
pub mod sparse;
pub mod topology;
#[cfg(feature = "blender-tessellation")]
pub mod triangulate;

pub use geometry::{Geometry, GeometryError};
pub use sparse::{Factorization, Method, Parallelism, SparseError, SparseMatrix, Symbolic, SymbolicInfo};
pub use topology::{
    CompactedVertices, MeshError, NONE, OrientedFaces, SubMesh, TriMesh, compact_vertices, orient_faces,
};
