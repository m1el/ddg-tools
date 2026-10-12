//! Boundary First Flattening with cone singularities, and the mesh cutting it needs.
//!
//! This crate is a Rust port of the method and reference implementation of
//!
//! > Rohan Sawhney and Keenan Crane. *Boundary First Flattening.*
//! > ACM Transactions on Graphics 37(1), 2017. <https://doi.org/10.1145/3132705>
//!
//! (reference C++ code: <https://github.com/GeometryCollective/boundary-first-flattening>, MIT,
//! © 2017 Rohan Sawhney and Keenan Crane; the licence notice is in `LICENSES/BFF-MIT.txt`).
//!
//! Given a triangle mesh and a set of cones (vertex, cone angle), [`flatten_with_cones`]
//! produces a conformal map of the cone metric to the plane: it solves for the log scale
//! factor, cuts the surface through the cones into disks, and flattens every disk by building
//! its boundary curve first and extending it harmonically. The flattenings of a single disk —
//! prescribed boundary scale factors, prescribed boundary angles, the automatic map of least
//! area distortion, the unit disk — are available through [`DiskBff`].
//!
//! * [`cut`] — cutting a mesh along edges, with index maps back to the original.
//! * [`plan`] — choosing the cut: every piece becomes a disk with all cones on its boundary;
//!   edges can be forbidden or required.
//! * [`disk`] — Boundary First Flattening of a disk.
//! * [`cones`] — the scale factor of a cone metric and the cone flattening.
//! * [`curve`] — boundary lengths and angles to a closed plane curve.
//! * [`quality`] — angle and area distortion, flipped triangles, overlap, similarity fit.
//! * [`mollify`] — intrinsic mollification of degenerate edge lengths (never applied unasked).
//!
//! The mollification technique is from N. Sharp and K. Crane, *A Laplacian for Nonmanifold
//! Triangle Meshes* (SGP 2020); the cut planner combines K. Mehlhorn's approximation of
//! Steiner trees (1988) with the tree–cotree decomposition of D. Eppstein (2003) and the
//! greedy loops of J. Erickson and K. Whittlesey (2005); see [`plan`].
//!
//! # Example
//!
//! ```
//! use ddg_bff::{Cone, ConeFlattenInput, ConeFlattenOptions, CutSource, flatten_with_cones};
//! use ddg_core::generators::icosphere;
//! use ddg_core::{Geometry, TriMesh};
//!
//! # fn main() -> Result<(), Box<dyn std::error::Error>> {
//! let sphere = icosphere(2)?;
//! let mesh = TriMesh::new(sphere.positions.len(), sphere.faces.clone())?;
//! let geometry = Geometry::from_positions(&mesh, &sphere.positions)?;
//! // The 12 vertices of the icosahedron, each with a cone angle of π/3: they sum to 4π.
//! let cones: Vec<Cone> = (0..12)
//!     .map(|vertex| Cone { vertex, angle: std::f64::consts::PI / 3.0 })
//!     .collect();
//! let flat = flatten_with_cones(
//!     &ConeFlattenInput {
//!         mesh: &mesh,
//!         geometry: &geometry,
//!         cones: &cones,
//!         scale_factors: None,
//!         cut: CutSource::Plan { may_cut: None, must_cut: &[] },
//!         anchors: &[],
//!     },
//!     &ConeFlattenOptions::default(),
//! )?;
//! // One UV per triangle corner, the edges that were cut, the chart of every face.
//! assert_eq!(flat.corner_uvs.len(), 3 * mesh.n_faces());
//! assert!(!flat.cut_edges().is_empty());
//! assert_eq!(flat.face_chart.len(), mesh.n_faces());
//! assert_eq!(flat.n_charts(), 1);
//! // The report: the two sides of every cut edge have the same length, nothing is flipped,
//! // nothing overlaps.
//! assert!(flat.seams.max_length_mismatch < 1e-12);
//! let chart = &flat.charts[0];
//! assert_eq!(chart.distortion.flipped_triangles, 0);
//! assert!(chart.is_injective());
//! assert_eq!(chart.overlap.overlapped_fraction(), Some(0.0));
//! # Ok(())
//! # }
//! ```
//!
//! # Conventions
//!
//! Those of the repository's `CLAUDE.md`: the cotan Laplacian `L` is positive semidefinite,
//! `u` is the true log scale factor (lengths in the plane are about `e^u` times lengths on the
//! surface), a cone angle is `µ = 2π − (angle sum around the vertex in the plane)`, and the
//! scale factor solves `(L u)_i = µ_i − Ω_i` with `u = 0` on the boundary. Indices are those of
//! `ddg_core::topology`: corner `3 f + k` is the corner of face `f` at `faces[f][k]`, and
//! faces and corners keep their ids under cutting.
//!
//! Bad input is reported through [`BffError`], whose messages name the vertex, edge, face or
//! chart concerned; nothing panics and no result contains NaN. Nothing is normalized behind
//! the caller's back: layouts are not recentred, rescaled or packed (every chart starts at the
//! origin, so charts lie on top of each other), edge lengths are not mollified, cone angles
//! that violate Gauss–Bonnet on a closed component are refused unless a policy says what to
//! do with them ([`GaussBonnetPolicy`]), and every place where a free constant had to be
//! chosen says how. All solves run on the calling thread.
//!
//! # What the result says about itself
//!
//! Per chart ([`ChartReport`]): how the boundary curve was closed ([`ClosureReport`]), angle
//! and area distortion with flipped and degenerate triangles ([`DistortionReport`]), and
//! overlap ([`OverlapReport`]: boundary vertices with more than a full turn around them,
//! crossings of the boundary, and an estimate of the overlapped area). Per cone
//! ([`ConeReport`]): the angle asked for, the angle the scale factors imply, the angle of the
//! layout, and its wedges. Across the cut ([`SeamReport`]): agreement of the lengths of the
//! two sides. A cone of negative angle has more than a full turn around it; where the cut
//! leaves all of it in one wedge the layout overlaps itself there, and the report says so
//! (see [`plan`], *Negative cones*).
//!
//! Functions that implement a formula carry a line of the form `[formula: id]` in their
//! documentation; the ids are listed with their statements in the crate README.
//!
//! # Differences from the reference implementation
//!
//! See the crate README for the complete list. The main ones: Poisson problems are solved
//! exactly instead of with a diagonal shift of `1e-8` (the shift is available as an option,
//! and with it the two programs agree to rounding error); edge lengths are not mollified
//! unless the caller does it; surfaces with several boundary loops are cut between the loops
//! instead of having their holes filled; cuts for handles are part of the same cut graph as
//! the cuts through cones and are seamless too; the cut planner accepts forbidden and
//! required edges; cone angles that violate Gauss–Bonnet are refused by default; the layout is
//! not recentred, rotated or packed.

#![warn(missing_docs)]

mod chart;
pub mod cones;
pub mod curve;
pub mod cut;
pub mod disk;
pub mod error;
mod laplace;
pub mod mollify;
pub mod plan;
pub mod quality;

pub use chart::{ChartReport, OverlapReport};
pub use cones::{
    ClosedComponentReport, Cone, ConeFlattenInput, ConeFlattenOptions, ConeFlattening, ConeReport,
    CutSource, GaussBonnetPolicy, ScaleFactorOptions, ScaleFactors, SeamReport, Timings,
    flatten_with_cones, solve_scale_factors,
};
pub use curve::ClosureReport;
pub use cut::{
    ComponentTopology, CutMesh, component_topology, cut_euler_characteristic, cut_mesh,
    discontinuous_edges,
};
pub use disk::{
    DEFAULT_OVERLAP_SAMPLES, DiskBff, DiskFlattening, DiskOptions, Extension,
    REFERENCE_DIAGONAL_SHIFT,
};
pub use error::BffError;
pub use mollify::{Mollified, REFERENCE_MOLLIFY_DELTA, mollify_edge_lengths};
pub use plan::{CutPlan, CutPlanInput, PlannedChart, plan_cut};
pub use quality::{Coverage, DistortionReport, Similarity, Stats, coverage, fit_similarity};
