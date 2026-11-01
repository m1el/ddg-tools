//! Surface simplification using intrinsic error metrics.
//!
//! An implementation of
//!
//! > Hsueh-Ti Derek Liu, Mark Gillespie, Benjamin Chislett, Nicholas Sharp, Alec Jacobson,
//! > Keenan Crane. *Surface Simplification using Intrinsic Error Metrics.* ACM Transactions on
//! > Graphics 42(4), Article 118, 2023. <https://doi.org/10.1145/3592403>
//!
//! with the authors' reference code (<https://github.com/HTDerekLiu/intrinsic-simplification>,
//! MIT, Copyright (c) 2023 Hsueh-Ti Derek Liu; notice in `LICENSES/`) as a reference. The
//! contract is `SPEC.md` in this crate.
//!
//! Work in progress: see `SPEC.md` §15 for the work packages.
//!
//! # What is here (WP0–WP2)
//!
//! * [`CoarseMesh`] (SPEC §2): the connectivity ([`ddg_intrinsic::Complex`], a Δ-complex),
//!   one length per edge, the angular coordinate of every halfedge at its tail with a
//!   reference halfedge per vertex (S4), the tracked points (the fine vertices, then any extra
//!   [`InputPoint`]s) with a per-face index, the fixed mask, and stable fine vertex ids (S2):
//!   every public vertex id is the input vertex id, the complex's renumbering is hidden.
//!   Constructed from an input mesh ([`CoarseMesh::from_positions`],
//!   [`CoarseMesh::from_lengths`], [`CoarseMesh::from_mesh`]); checked by
//!   [`CoarseMesh::certify`] (SPEC §2.4).
//! * Edge flips with point tracking ([`CoarseMesh::flip_edge`], SPEC §3.2, S16) and Delaunay
//!   flipping ([`CoarseMesh::flip_to_delaunay`], [`CoarseMesh::flip_to_delaunay_from`], §3.3).
//! * Vertex removal (paper Sec. 4; SPEC §3.4–3.6): flattening by a bracketed Newton
//!   iteration with Springborn's flips as the fallback, flips to degree three (two on the
//!   boundary; ears flipped first), the removal from the complex with the points of the star
//!   re-expressed, local Delaunay flipping; in two phases
//!   ([`CoarseMesh::begin_removal`], then [`CoarseMesh::finish_removal`] or
//!   [`CoarseMesh::revert_removal`]) on an undo log, so that the tentative removal of the
//!   error metric ([`CoarseMesh::tentative_removal`]) runs the same code as the committed one
//!   and leaves the mesh bitwise unchanged (S5). [`PreparedRemoval`] carries what the ICE
//!   metric (WP3) needs.
//! * [`baseline`]: the extrinsic baseline of Sec. 8.4, quadric error metrics (Garland, Heckbert 1997)
//!   with a closest-point projection map (SPEC §10.4).

pub mod baseline;
mod error;
mod flatten;
mod flip;
pub mod formulas;
mod mesh;
mod removal;
mod undo;

pub use error::{IemError, RemovalFailure, RemovalRefusal};
pub use flatten::{FLATTEN_MAX_FACTOR, FlattenReport};
pub use flip::{DelaunayFlipReport, FlippedEdge};
pub use mesh::{
    ANGLE_CHECK_TOLERANCE, Certificate, CoarseMesh, CoarseMeshOptions, InputPoint, PointLocation,
};
pub use removal::{NeighbourChange, PreparedRemoval, RemovalKind, RemovalReport};
