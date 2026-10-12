//! Foundation crate of the `ddg-tools` workspace: everything the discrete differential geometry
//! algorithms share.
//!
//! * [`topology`] — a validated, oriented, manifold triangle mesh ([`TriMesh`]) with halfedge
//!   connectivity, boundary loops, components, and explicit repair helpers.
//! * `triangulate` — polygon triangulation that matches Blender's, with a map from triangle
//!   corners back to the original polygon corners (cargo feature `blender-tessellation`).
//! * [`geometry`] — intrinsic triangle geometry (corner angles, areas, cotangents) computed from
//!   edge lengths ([`Geometry`]).
//! * [`operators`] — the cotan Laplacian, lumped mass, angle defects and boundary geodesic
//!   curvature, following the conventions in the repository's `CLAUDE.md`.
//! * [`sparse`] — a small sparse matrix type and direct solvers with reusable symbolic analysis.
//! * [`io`] — Wavefront OBJ reading and writing, and the cone files of the optimal-cones dataset.
//! * [`generators`] — small procedural test meshes.
//!
//! # Example
//!
//! ```
//! use ddg_core::generators::icosphere;
//! use ddg_core::operators::{angle_defects, cotan_laplacian, lumped_mass_normalized};
//! use ddg_core::sparse::{Method, Parallelism, Symbolic};
//! use ddg_core::{Geometry, TriMesh};
//!
//! # fn main() -> Result<(), Box<dyn std::error::Error>> {
//! let sphere = icosphere(2)?;
//! let mesh = TriMesh::new(sphere.positions.len(), sphere.faces.clone())?;
//! assert_eq!((mesh.euler_characteristic(), mesh.genus()), (2, 0));
//!
//! let geometry = Geometry::from_positions(&mesh, &sphere.positions)?;
//! let mut a = cotan_laplacian(&mesh, &geometry)?; // symmetric PSD, constants in the kernel
//! let mass = lumped_mass_normalized(&mesh, &geometry)?; // sums to 1
//! a.add_to_diagonal(&mass)?; // L + M is positive definite
//!
//! let symbolic = Symbolic::analyze(&a, Method::Llt)?; // reusable for this sparsity pattern
//! let factorization = symbolic.factor(&a, Parallelism::Sequential)?;
//! let mut x = angle_defects(&mesh, &geometry)?; // right-hand side, overwritten by the solution
//! factorization.solve_in_place(&mut x, Parallelism::Sequential)?;
//! # Ok(())
//! # }
//! ```
//!
//! # References
//!
//! The discrete operators are those of Sec. 5.1 of
//! Yousuf Soliman, Dejan Slepčev, Keenan Crane, *Optimal Cone Singularities for Conformal
//! Flattening*, ACM Transactions on Graphics 37(4), 2018,
//! <https://doi.org/10.1145/3197517.3201367>. The cotan operator goes back to MacNeal (1949).
//! The numerically careful triangle area and angle formulas are from W. Kahan, *Miscalculating
//! Area and Angles of a Needle-like Triangle* (2014). Polygon triangulation is a port of
//! Blender's face tessellation (<https://www.blender.org>), and the sparse factorizations are
//! those of the `faer` crate.
//!
//! # Conventions
//!
//! `f64` everywhere, stored indices are `u32`, ids at function boundaries are `usize`. Data from
//! outside (face lists, positions, lengths, files, matrices) is validated and rejected through
//! error enums whose messages say what is wrong and where; it never causes a panic. Accessors
//! that take an element id index like slices. Nothing is normalized, clamped or repaired
//! silently. Functions that implement a mathematical formula carry a `[formula: <id>]` line in
//! their documentation.
//!
//! # License
//!
//! MIT, except for the private module `triangulate::blender`, which is ported from Blender and
//! licensed under the GNU GPL, version 2 or later. It is compiled only with the cargo feature
//! `blender-tessellation` (on by default); see the crate README.

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
pub use sparse::{SparseError, SparseMatrix};
pub use topology::{MeshError, TriMesh};
