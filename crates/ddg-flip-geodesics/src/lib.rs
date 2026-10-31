//! Geodesic paths, loops and networks on triangle meshes by edge flips.
//!
//! An implementation of
//!
//! > Nicholas Sharp, Keenan Crane. *You Can Find Geodesic Paths in Triangle Meshes by Just
//! > Flipping Edges.* ACM Transactions on Graphics 39(6), Article 249, 2020.
//! > <https://doi.org/10.1145/3414685.3417839>
//!
//! The contract of this crate is `SPEC.md` next to its `Cargo.toml`.
//!
//! Everything is generic over [`ddg_intrinsic::Triangulation`] (F1), so it runs on the
//! signposts of `ddg-intrinsic` and on the integer coordinates of `ddg-integer-coords`
//! alike: it reads the connectivity ([`ddg_intrinsic::Complex`]) and the edge lengths, and
//! changes the triangulation only by its flips.
//!
//! # What is here
//!
//! * Edge networks ([`Network`], [`GeodesicNetwork`]): paths and loops as sequences of
//!   directed edges ([`Segment`], [`EdgePath`]), the order of the segments that share an edge
//!   (the layering, paper Fig. 11), pinned vertices, and the validity check including
//!   non-crossing ([`Network::check`]).
//! * Initial paths: through a vertex sequence ([`path_from_vertices`]), the shortest path in the
//!   edge graph ([`dijkstra_path`]), through waypoints ([`waypoint_path`]) and from a set of
//!   edges with junctions ([`network_from_edges`]; seams, cut graphs).
//! * The geometry of a joint ([`JointGeometry`]): both wedge angles summed corner by corner
//!   (Eq. 1), boundary wedges, the Δ-complex cases (self-edges incident on a wedge twice,
//!   out-and-back joints of degree zero), the classification of SPEC §3, flexible and
//!   blocked wedges.
//! * FlipOut ([`GeodesicNetwork::shorten`], modules [`shorten`] and [`flipout`], SPEC §4):
//!   Alg. 1 on one wedge (flips until the outer arc is an edge path, the out-and-back
//!   contraction, the App. B move of a single-segment loop;
//!   [`GeodesicNetwork::shorten_joint`]) and Alg. 2 on paths, loops and networks (the queue
//!   of F4, stale and blocked entries, stalls, the final sweep, the stop rules), with a
//!   report ([`ShortenReport`]) and an observer of every iteration
//!   ([`GeodesicNetwork::shorten_with`]).
//! * The geodesic certificate ([`certify()`], SPEC §7 layer 1a), computed on the result.
//!
//! * Extraction on the input mesh ([`GeodesicNetwork::extract`], module [`extract`], SPEC §5):
//!   one polyline per path, each segment the trace of its edge (parallel per edge, assembled in
//!   order), the certificate on `T⁰` computed from the input mesh alone ([`InputCertificate`]),
//!   3D polylines ([`GeodesicNetwork::polylines_3d`]).
//!
//! * Geodesic Bézier curves by midpoint subdivision ([`geodesic_bezier`],
//!   [`GeodesicNetwork::bezier_subdivide`], module [`bezier`], SPEC §6.1), with edge splits
//!   that carry the network along ([`GeodesicNetwork::split_edge`]).
//! * The intrinsic constrained Delaunay triangulation of a network
//!   ([`GeodesicNetwork::make_delaunay`], module [`cdt`], SPEC §6.3).
//!
//! # Example
//!
//! ```
//! use ddg_core::Geometry;
//! use ddg_core::generators::grid_disk;
//! use ddg_flip_geodesics::{
//!     CertifyOptions, GeodesicNetwork, ShortenOptions, ShortenStatus, dijkstra_path,
//!     path_from_vertices,
//! };
//! use ddg_intrinsic::IntrinsicTriangulation;
//!
//! # fn main() -> Result<(), Box<dyn std::error::Error>> {
//! // A flat 4 x 4 grid; vertex (i, j) is 5 i + j.
//! let shape = grid_disk(4, 4)?;
//! let mesh = shape.topology()?;
//! let geometry = Geometry::from_positions(&mesh, &shape.positions)?;
//! let tri = IntrinsicTriangulation::new(&mesh, &geometry)?;
//!
//! // The shortest edge path between two corners, and a path that turns at (2, 2).
//! let diagonal = dijkstra_path(&tri, 0, 24)?;
//! let turn = path_from_vertices(&tri, &[7, 12, 13], false)?;
//! let mut network = GeodesicNetwork::new(tri, &[diagonal, turn], &[], None)?;
//!
//! // Not geodesic yet: the turn at (2, 2) bends by 90 degrees.
//! let certificate = network.certify(&CertifyOptions::default())?;
//! assert!(!certificate.geodesic);
//! assert!(certificate.failures.iter().any(|f| f.vertex == 12));
//!
//! // FlipOut: flip edges until every joint is straight; the paths stay edge paths.
//! let before = network.length();
//! let report = network.shorten(&ShortenOptions::default())?;
//! assert_eq!(report.status, ShortenStatus::Geodesic);
//! assert!(network.length() < before);
//! assert!(network.certify(&CertifyOptions::default())?.geodesic);
//! # Ok(())
//! # }
//! ```
//!
//! # The reference implementation
//!
//! The authors' implementation in geometry-central (© 2017–2019 Nicholas Sharp and the
//! geometry-central contributors, MIT; notice in `LICENSES/geometry-central-MIT.txt`) was read
//! for the details the paper leaves out; places that follow it are cited `gc: <file>:<line>`
//! (SPEC header). Taken from it: the default layering, the outermost-segment convention, the
//! test for a wedge blocked by a path end, the tolerances and Dijkstra's tie order, the flip
//! scan of Alg. 1 with its length guard, the move of a single-segment loop, stale entries.
//! Different: wedge angles are summed from corner angles rather than read from signposts (F2),
//! the out-and-back joint is resolved by the layering and contracted (F12, F26), crossing
//! inputs are refused, `closed` and pinning are explicit (F9, F10), the queue has a total
//! order (F4), a final sweep follows the empty queue (F5), edges carrying segments are never
//! flipped by construction (F7), and the stop rules include an iteration limit for paths
//! too (F13).
//!
//! # Conventions
//!
//! Those of the workspace (`CLAUDE.md`): `f64`, stored indices `u32`, ids at function
//! boundaries `usize`; errors name the element; nothing is repaired or normalized silently.

#![warn(missing_docs)]

pub mod bezier;
pub mod cdt;
pub mod certify;
mod error;
pub mod extract;
pub mod flipout;
pub mod init;
pub mod network;
pub mod report;
mod scratch;
pub mod shorten;
pub mod wedge;

pub use bezier::{
    BezierMidpoint, BezierOptions, BezierReport, MAX_BEZIER_ROUNDS, MidpointKind, geodesic_bezier,
};
pub use cdt::CdtReport;
pub use certify::{CertifyOptions, GeodesicCertificate, JointFailure, JointStatus, certify};
pub use error::FlipError;
pub use extract::{
    EdgeExtraction, ExtractOptions, ExtractReport, Extraction, InputCertificate, crossing_sines,
    polyline_3d,
};
pub use init::{
    InitialPaths, dijkstra_path, network_from_edges, path_from_vertices, waypoint_path,
};
pub use network::{
    EdgePath, GeodesicNetwork, Joint, Layering, Network, NetworkSplit, Segment, SegmentId,
    SegmentInfo, SegmentRef,
};
pub use report::{ShortenOptions, ShortenReport, ShortenStatus, ShortenStep, Shortening};
pub use shorten::ShortenSession;
pub use wedge::{EdgeEnd, JointGeometry, Side, Turn, VertexStar, Wedge, wedge_angle};
