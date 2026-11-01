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
//! # What is here (WP0–WP3)
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
//! * The ICE error metric ([`metric`], SPEC §4: channels, Eq. 5–9, memoryless or memory-based)
//!   and intrinsic coarsening ([`simplify()`], paper Alg. 1, SPEC §5: the queue, targets, fixed
//!   vertices, snapshots, the report), with [`SimplifyOptions::reference_compat`] for
//!   differential tests against the reference (SPEC §12.2).
//! * [`baseline`]: the extrinsic baseline of Sec. 8.4, quadric error metrics (Garland, Heckbert 1997)
//!   with a closest-point projection map (SPEC §10.4).
//! * [`curvature`]: principal curvatures (Panozzo, Puppo, Rocca 2010) for anisotropic
//!   coarsening (SPEC §8), and 3D directions as angles in the intrinsic frames.
//! * [`CoarseMesh::compact`] and the map (WP4, SPEC §7): [`CoarseOutput`] with every tracked
//!   point in the compacted coarse ids, the scalar prolongation `P`
//!   ([`CoarseOutput::scalar_prolongation`]), the vector prolongation `P^vec`
//!   ([`CoarseOutput::vector_prolongation`], module [`prolongation`]) and the greedy colouring of
//!   the coarse faces ([`CoarseOutput::colour_faces`]).
//! * [`anisotropy`]: anisotropic edge-length scaling from principal directions (paper Sec. 8.6.2,
//!   SPEC §8) and the coarse mesh of the scaled metric.
//! * [`applications`] (WP6, SPEC §10.1): Poisson, single-source geodesic distance (exact and heat),
//!   all pairs, Voronoi and the exponential map on the coarse mesh, prolongated with `P`.
//! * [`refine`] (WP5, SPEC §6): pre-refinement ([`CoarseMesh::from_refined`],
//!   [`CoarseMesh::from_intrinsic`]) and post-refinement ([`CoarseMesh::post_refine`]) with the
//!   workspace's Delaunay refinement, tracked points composed through the `T⁰ → T¹` query.
//! * [`multigrid`]: the Galerkin multigrid solver of Liu et al. 2021 (SPEC §9.2).

pub mod anisotropy;
pub mod applications;
pub mod baseline;
pub mod curvature;
mod error;
mod flatten;
mod flip;
pub mod formulas;
mod mesh;
pub mod metric;
pub mod multigrid;
mod output;
pub mod prolongation;
pub mod refine;
mod removal;
mod simplify;
mod undo;

pub use error::{IemError, RemovalFailure, RemovalRefusal};
pub use flatten::{FLATTEN_MAX_FACTOR, FlattenReport};
pub use flip::{DelaunayFlipReport, FlippedEdge};
pub use mesh::{
    ANGLE_CHECK_TOLERANCE, Certificate, CoarseMesh, CoarseMeshOptions, InputPoint, PointLocation,
};
pub use metric::{
    BoundaryFrames, Channel, ChannelSource, ChannelState, ChannelWeights, FlatRule, Frames,
    MetricOptions,
};
pub use output::{CoarseOutput, CoarsePoint, FaceColouring};
pub use prolongation::{
    C64, ComplexSparse, FineFrames, FineNeighbour, ProlongationError, VectorProlongationReport,
    normalize, unnormalize,
};
pub use removal::{NeighbourChange, PreparedRemoval, RemovalKind, RemovalReport};
pub use simplify::{
    ComponentCount, KeyCheck, ReferenceCompat, RemovalRecord, Requeue, SimplifyOptions,
    SimplifyReport, SimplifyTimings, Snapshot, StopReason, Target, refusal_name, simplify,
};
