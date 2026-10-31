//! Intrinsic triangulations stored by integer coordinates.
//!
//! An implementation of
//!
//! > Mark Gillespie, Nicholas Sharp, Keenan Crane. *Integer Coordinates for Intrinsic
//! > Geometry Processing.* ACM Transactions on Graphics 40(6), Article 252, 2021.
//! > <https://doi.org/10.1145/3478513.3480522>
//!
//! The contract of this crate is `SPEC.md` next to its `Cargo.toml`.
//!
//! An intrinsic triangulation `T¹` of an input mesh `T⁰` is stored as its connectivity (the
//! [`Complex`](ddg_intrinsic::Complex) shared with the signpost representation of
//! `ddg_intrinsic`), its edge lengths, and integers that say how the input edges run over it:
//! one **normal coordinate** `n` per edge (the number of input edges crossing it, `−1` if it runs
//! along one) and one **roundabout** `r` per halfedge leaving an input vertex (which input
//! halfedge comes next counter-clockwise), plus the place `q⁰` of every inserted vertex. The
//! connectivity of the common subdivision is recovered from the integers exactly.
//!
//! # What is here (work package WP1)
//!
//! * [`IntegerTriangulation`]: the data, its initialization from an input mesh (with optional
//!   intrinsic mollification, reported), read-only accessors, and the derived quantities of
//!   every face: corner coordinates (Eq. 1) and emanating counts (Eq. 2)
//!   ([`IntegerTriangulation::face_coordinates`]).
//! * [`formulas`]: the integer formulas on plain integers, compared with the Lean oracle by
//!   exact equality.
//! * [`trace`]: the combinatorial part of Supp. Alg. 1–4: [`trace_from`] (`TraceFrom`),
//!   [`extract_curve`] (`ExtractCurve`), [`IntegerTriangulation::identify_traced_edge`]
//!   (`IdentifyTracedEdge`, with the walk back past inserted vertices),
//!   [`IntegerTriangulation::extract_edge`] (`ExtractEdge`).
//! * [`IntegerTriangulation::certify`]: the certificate of SPEC §2.4.
//! * [`IntegerTriangulation::from_intrinsic`]: integer coordinates read off a signpost
//!   triangulation by tracing its edges (floating point decides it; the result is certified).
//!
//!
//! # Operations (work package WP3)
//!
//! * [`flip`]: [`IntegerTriangulation::flip_edge`] (Eq. 3 with the corrected first terms,
//!   Eq. 4), [`IntegerTriangulation::can_flip`], and Delaunay flipping (the shared
//!   [`Triangulation::flip_to_delaunay`](ddg_intrinsic::Triangulation::flip_to_delaunay);
//!   [`IntegerTriangulation::flip_to_delaunay_certified`] certifies the result).
//! * [`split`]: [`IntegerTriangulation::split_edge`] for edges along input edges (Eq. 6);
//!   [`IntegerTriangulation::split_edge_given`] for crossed edges with the geometric decisions
//!   supplied (SPEC I20).
//! * [`remove`]: [`IntegerTriangulation::remove_vertex`] (Alg. 8, I8, I24, I27).
//! * [`interface`]: `impl Triangulation for IntegerTriangulation`.
//!
//! # Geometry of the crossings (work package WP2)
//!
//! * [`geometry`]: the positions of the crossings, from straight segments through unfolded
//!   strips (clamps counted, SPEC I5): [`IntegerTriangulation::extract_geometric_curve`]
//!   (Alg. 2), [`IntegerTriangulation::extract_geometric_crossing`],
//!   [`IntegerTriangulation::extract_edge_geometry`] (Alg. 4 with positions), and
//!   [`IntegerTriangulation::extract_intrinsic_edge`] (`ExtractIntrinsicEdge`, the path of an
//!   edge over the input mesh; the shared interface's `trace_edge_on_input`).
//!
//! Not yet: face splits (WP4), splits of crossed edges through the shared interface (they need
//! the geometric crossings of WP2 wired into the split), the common subdivision (WP5).
//!
//! # The reference implementation
//!
//! geometry-central (MIT, © 2017–2019 Nicholas Sharp and the geometry-central contributors) at
//! commit `ab116f1` was read for the details the paper leaves out; nothing was translated from it.
//! Its notice is in this crate's `LICENSES/`.
//!
//! # Conventions
//!
//! Those of the workspace (`CLAUDE.md`) and of `ddg_intrinsic` (halfedge `h = 3f + k` from corner
//! `k` to corner `k + 1`). Normal coordinates are `i64` with checked arithmetic (SPEC I2); Eq. 1
//! never rounds an odd numerator (I14). Bad data is an error naming the face, edge, halfedge,
//! vertex or crossing, never a panic.

#![warn(missing_docs)]

pub mod certify;
pub mod convert;
pub mod coords;
mod error;
pub mod flip;
pub mod formulas;
pub mod geometry;
pub mod interface;
pub mod remove;
pub mod split;
pub mod trace;

pub use certify::CertifyReport;
pub use convert::ConversionReport;
pub use coords::{
    InitOptions, InitReport, IntegerTriangulation, LocalIndices, NO_ROUNDABOUT,
    PAPER_MOLLIFY_EPSILON, Parts,
};
pub use error::{IntegerError, InvalidCoordinates, Place, ValidityItem};
pub use formulas::FaceCoordinates;
pub use geometry::{GeometricCrossing, GeometricCurve, InputEdgeGeometry, StripReport};
pub use interface::{REPRESENTATION, to_intrinsic_error};
pub use remove::Removal;
pub use trace::{
    Arrival, Crossing, CurveEnd, CurveSegment, HalfTrace, InputEdgePath, NormalGluing, PathPiece,
    TracedInputEdge, extract_curve, trace_from, trace_with,
};
