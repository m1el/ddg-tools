//! Intrinsic triangulations of a triangle mesh, with the signpost data structure of
//!
//! > Nicholas Sharp, Yousuf Soliman, Keenan Crane. *Navigating Intrinsic Triangulations.*
//! > ACM Transactions on Graphics 38(4), 2019. <https://doi.org/10.1145/3306346.3322979>
//!
//! An *intrinsic triangulation* is another triangulation of the same polyhedral surface: its
//! vertices are points of the input surface, its edges are straight paths on it, its triangles
//! are flat. Flipping edges until the triangulation is Delaunay leaves the surface exactly as
//! it is and gives non-negative cotan weights; that is what the rest of this workspace uses it
//! for. What is in scope, the conventions and every decision the paper leaves open are in
//! `SPEC.md` next to this crate's `Cargo.toml`.
//!
//! # What is here
//!
//! * [`IntrinsicTriangulation`]: the halfedge mesh that flips can produce (double edges,
//!   self-edges, triangles that use a vertex twice, vertices of degree one), edge lengths,
//!   signposts and the place of every vertex on the input mesh ([`SurfacePoint`]).
//! * [`IntrinsicTriangulation::flip_edge`] with every refusal ([`FlipRefusal`]) and
//!   [`IntrinsicTriangulation::flip_to_delaunay`] with a report ([`DelaunayReport`]).
//! * Tracing queries over the input mesh ([`InputSurface::trace`]) and over the intrinsic
//!   triangulation ([`IntrinsicTriangulation::trace`]), from a vertex, a point of an edge or a
//!   point of a face ([`Trace`]).
//! * The path of every intrinsic edge over the input mesh
//!   ([`IntrinsicTriangulation::trace_edges`]), checked: a trace that does not arrive at the
//!   edge's far vertex is an error.
//! * The cotan Laplacian and the lumped vertex areas of the triangulation
//!   ([`IntrinsicTriangulation::cotan_laplacian`], [`IntrinsicTriangulation::lumped_mass`]),
//!   with the conventions of `ddg_core::operators`.
//!
//! * Vertex insertion and edge split ([`IntrinsicTriangulation::insert_vertex`],
//!   [`IntrinsicTriangulation::split_edge`]), point queries both ways
//!   ([`IntrinsicTriangulation::point_on_input`],
//!   [`IntrinsicTriangulation::point_on_triangulation`]).
//! * The overlay, the common subdivision with the input mesh
//!   ([`IntrinsicTriangulation::overlay`]), and a straight segment inside a face as its pieces
//!   on input faces ([`IntrinsicTriangulation::segments_on_input`]).
//! * Delaunay refinement to an angle bound ([`IntrinsicTriangulation::refine`]); the report
//!   lists every corner left below the bound with the reason ([`RefineReport`]).
//! * The triangulation as a simplicial `TriMesh` with a `Geometry` for code that takes one
//!   ([`IntrinsicTriangulation::export_simplicial`]), and what identifies that state
//!   ([`SimplicialExport::check_state`]). Splitting does not keep an angle bound; the export
//!   reports the angles it has.
//!
//! * Vertex repositioning inside the star of an inserted vertex
//!   ([`IntrinsicTriangulation::move_vertex`], paper Sec. 3.3.3) and the intrinsic optimal
//!   Delaunay triangulation ([`IntrinsicTriangulation::optimal_delaunay`], Sec. 4.3, after Chen
//!   and Holst 2011), with angle and area statistics
//!   ([`IntrinsicTriangulation::shape_statistics`]).
//!
//! * The connectivity on its own, [`Complex`]: the halfedge Δ-complex with its topological
//!   operations (flip, face split, edge split, removal of a vertex of degree three), each with
//!   its index bookkeeping, and a mask of constrained edges that flips leave alone. The
//!   [`Triangulation`] trait is the interface this representation shares with integer
//!   coordinates (Gillespie, Sharp, Crane 2021); Delaunay flipping is written once over it.
//!   [`AngleFrames`] traces over a triangulation that stores no signposts, computing them
//!   from the corner angles (SPEC N50).
//!
//! Not implemented: the paper's general route for a move beyond the star (SPEC N20); only
//! inserted vertices of degree three can be removed
//! ([`IntrinsicTriangulation::remove_vertex`]); refinement with constrained edges.
//!
//! # Example
//!
//! ```
//! use ddg_core::generators::grid_disk;
//! use ddg_core::Geometry;
//! use ddg_intrinsic::{DelaunayOptions, DelaunayStatus, EdgeTraceOptions, IntrinsicTriangulation};
//!
//! # fn main() -> Result<(), Box<dyn std::error::Error>> {
//! // A flat grid, sheared so that every triangle has an angle of 124 degrees.
//! let mut shape = grid_disk(8, 8)?;
//! for p in &mut shape.positions {
//!     p[0] += 1.5 * p[1];
//! }
//! let mesh = shape.topology()?;
//! let geometry = Geometry::from_positions(&mesh, &shape.positions)?;
//!
//! let mut triangulation = IntrinsicTriangulation::new(&mesh, &geometry)?;
//! let report = triangulation.flip_to_delaunay(&DelaunayOptions::default())?;
//! assert_eq!(report.status, DelaunayStatus::Delaunay);
//! assert!(report.n_flips > 0);
//! // No interior edge has a negative cotan weight any more.
//! assert!(report.min_cotan_weight > -1e-9);
//!
//! // It is the same surface: same area, same cone angles.
//! let consistency = triangulation.check()?;
//! assert!(consistency.relative_area_error.abs() < 1e-12);
//! assert!(consistency.max_angle_sum_error < 1e-12);
//! let laplacian = triangulation.cotan_laplacian(); // positive semidefinite, as in ddg_core
//! assert_eq!(laplacian.nrows(), mesh.n_vertices());
//!
//! // Where the flipped edges run on the input mesh: the input edges each of them crosses.
//! let traces = triangulation.trace_edges(&EdgeTraceOptions::default())?;
//! assert!(traces.report.n_crossings >= traces.report.n_traced);
//! assert!(traces.report.max_relative_mismatch < 1e-12);
//! # Ok(())
//! # }
//! ```
//!
//! # The reference implementation
//!
//! The authors' implementation (the demo `navigating-intrinsic-triangulations-demo`, © 2019
//! Nick Sharp, on geometry-central, © 2017–2019 Nicholas Sharp and the geometry-central
//! contributors; both MIT, notices in this crate's `LICENSES/`) was read for the details the
//! paper leaves out. No code was translated from it. Taken from it are conventions and
//! constants: the normalization of angles at boundary vertices (SPEC N2), the cotan-weight
//! form of the Delaunay test with its tolerance `1e-6`
//! ([`DelaunayCriterion::REFERENCE`]), the order of the flip queue, and the flag for edges
//! that still lie along an input edge. The tracing is a different algorithm (see
//! [`trace`]), as are the length of a flipped edge and the convexity test of a flip
//! ([`formulas`]).
//!
//! # Conventions
//!
//! Those of the workspace (`CLAUDE.md`): `f64`, stored indices `u32`, ids at function
//! boundaries `usize`. Accessors that take an element id index like slices. Operations return
//! error enums that name the element; an operation that can do a variable amount of work
//! returns a report. Nothing is perturbed, repaired or normalized silently.

#![warn(missing_docs)]

pub mod complex;
pub mod delaunay;
pub mod edge_trace;
mod error;
pub mod export;
mod flip;
pub mod formulas;
pub mod frames;
pub mod generic_ops;
mod input;
pub mod insert;
pub mod interface;
pub mod merge;
pub mod odt;
mod operators;
pub mod overlay;
mod point;
pub mod query;
pub mod refine;
pub mod remove;
pub mod reposition;
mod surface;
pub mod trace;
pub mod triangulation;

pub use complex::Complex;
pub use delaunay::{DelaunayCriterion, DelaunayOptions, DelaunayReport, DelaunayStatus, FlipLimit};
pub use edge_trace::{EdgeTrace, EdgeTraceOptions, EdgeTraceReport, EdgeTraces};
pub use error::{
    EdgeTraceFailure, EdgeTraceFailureReason, FlipRefusal, InsertRefusal, IntrinsicError,
    MoveRefusal, RemovalRefusal,
};
pub use export::{ExportOptions, SimplicialExport};
pub use flip::{Flip, FlipPlan, FlipRequest, FlipShapeTest, plan_flip, plan_flip_with};
pub use frames::AngleFrames;
pub use input::InputSurface;
pub use insert::{EdgeSplitPlan, FaceSplitPlan, Insertion, plan_edge_split, plan_face_split};
pub use interface::Triangulation;
pub use odt::{
    BoundaryVertices, OdtCenters, OdtIteration, OptimalDelaunayOptions, OptimalDelaunayReport,
    OptimalDelaunayStatus, TargetLength,
};
pub use operators::{Quality, ShapeStatistics, Summary};
pub use overlay::{Overlay, OverlayVertex};
pub use point::SurfacePoint;
pub use query::{PathOnInput, PathsOnInput, Segment};
pub use refine::{
    BelowBound, BelowReason, InsertionLimit, RefineOptions, RefineReport, RefineRules, RefineStatus,
};
pub use merge::{EdgeMergePlan, MERGE_STRAIGHT_TOLERANCE, plan_edge_merge, straight_pair};
pub use remove::{InsertedVertexRemoval, RemovalStep, removal_flip_choice};
pub use reposition::VertexMove;
pub use trace::{Crossing, Trace, TraceOptions, TraceStop};
pub use triangulation::{ConsistencyReport, IntrinsicTriangulation, Irregularities};
