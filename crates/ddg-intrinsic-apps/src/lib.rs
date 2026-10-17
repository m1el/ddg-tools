//! Applications of intrinsic triangulations, Sec. 5 of
//!
//! > Nicholas Sharp, Yousuf Soliman, Keenan Crane. *Navigating Intrinsic Triangulations.*
//! > ACM Transactions on Graphics 38(4), 2019. <https://doi.org/10.1145/3306346.3322979>
//!
//! built on the crate `ddg-intrinsic` (the signpost data structure, flips, Delaunay refinement,
//! the optimal Delaunay triangulation, tracing and point queries). Each application is a module
//! of its own; the decisions the paper leaves open are in `SPEC.md`, the measurements in
//! `README.md`, both next to this crate's `Cargo.toml`.
//!
//! * [`steiner`] (paper Sec. 5.1): a short cut through all vertices, an approximate Steiner
//!   tree of the vertex set over the intrinsic Delaunay triangulation with Fermat points
//!   inserted as intrinsic vertices (Smith, Lee and Liebman 1981), its paths on the input mesh,
//!   and the check that cutting along it leaves a disk that unfolds into the plane.
//! * [`fem`] (Sec. 5.2): stiffness and mass matrices of piecewise linear finite elements on an
//!   intrinsic triangulation, the regular 4-1 subdivision of the input as the extrinsic
//!   comparison, and the extreme eigenvalues and condition number of the stiffness matrix
//!   (Lanczos, with or without a mass matrix).
//! * [`heat`] (Sec. 5.3): geodesic distance by the heat method (Crane, Weischedel and Wardetzky
//!   2013) on an intrinsic triangulation, factored once and reused for every source, with
//!   values at the vertices and at any point of the input mesh.
//! * [`exact`]: the exact polyhedral distance from vertices (window propagation after Mitchell,
//!   Mount and Papadimitriou 1987 and Surazhsky et al. 2005, with the vertex filter of Xin and
//!   Wang 2009). It is the ground truth against which Sec. 5.3 is measured.
//! * [`amr`] (Sec. 5.4): adaptive mesh refinement driven by local a posteriori error
//!   estimates (Morin, Nochetto, Siebert 2002; Mekchay, Nochetto 2005), for the harmonic
//!   Green's function, the short-time heat kernel and a harmonic map of a disk.
//! * [`vector`] (Sec. 5.5): connection Laplacians on an intrinsic triangulation (finite
//!   differences, Sharp, Soliman, Crane 2019 "The Vector Heat Method", Sec. 5.3; finite
//!   elements, Knöppel, Crane, Pinkall, Schröder 2013, Sec. 6.1.1), harmonic vector fields
//!   with boundary data and the vector maximum principle, smoothest `n`-direction fields.
//!
//! # Conventions
//!
//! Those of the workspace (`CLAUDE.md`): `f64`, stored indices `u32`, ids at function
//! boundaries `usize`; the cotan Laplacian `L` is positive semidefinite. Operations take
//! option structs, return reports, and fail with error enums that name the element; nothing is
//! clamped or normalized silently.

#![warn(missing_docs)]

pub mod amr;
pub mod exact;
pub mod fem;
pub mod heat;
pub mod steiner;
pub mod vector;
