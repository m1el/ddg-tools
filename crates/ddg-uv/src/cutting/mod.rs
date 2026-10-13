//! The cut flow on polygon meshes, end to end.
//!
//! "Variational Surface Cutting" (Sharp, Crane 2018, crate `ddg-variational-cuts`, implemented
//! clean-room from the text of the paper) moves a cut on a triangle mesh toward low flattening
//! distortion. This module makes it usable from polygon meshes as an editor provides them:
//! polygons, a triangulation of them, marked seams.
//!
//! ```text
//! Session::new(mesh, triangulation, seams, start)      the flow on the triangulation cut open at the seams
//!   .step / .run(problem, options)                     the flow; η and symmetry classes per polygon-mesh vertex
//!   .state()  /  .set_state(..)                        region functions per polygon corner
//!   .preview(δ)                                        the cut as polylines in space
//!   .describe_cut(δ) -> CutOnPolygons                  edge points and per-polygon paths
//! insert(mesh, &CutOnPolygons) -> Inserted             the crossed polygons split along the paths
//! unwrap_with_cut(..)                                  insert, then `unwrap` with seams = old seams + cut
//! ```
//!
//! The pieces:
//!
//! * [`Session`] owns the polygon mesh, the triangle mesh the flow runs on (the triangulation
//!   cut along the seams, so the seams are boundary for the flow, as `prepare` cuts them), the
//!   maps between polygon corners and the vertices of that mesh, and a
//!   [`Flow`](ddg_variational_cuts::Flow). The state crosses the API as **region functions per
//!   polygon corner**, in the caller's length unit; corners of one vertex that no seam
//!   separates must carry the same value (an error names the vertex). With two regions the
//!   functions are `φ` and `−φ`.
//! * [`Start`] says how to begin: a geodesic circle, polygon labels, or corner values.
//! * [`CutOnPolygons`] is the cut on the caller's polygons ([`EdgePoint`], [`PolygonPath`]).
//! * [`insert`] builds a new polygon mesh in which the crossed polygons are split along the
//!   paths, with the new seam edges and the maps to the old mesh ([`Inserted`]).
//! * [`unwrap_with_cut`] inserts and unwraps; [`ConePlacement`] says whether `unwrap` may place
//!   cones as well.
//!
//! More than two regions are not available yet: the core crate does not have them. The session
//! keeps regions as an array of functions and converts to and from the flow's level-set function
//! in two places only (`Maps::phi_from_regions`, `Maps::regions_from_phi`), and the start from
//! labels in one; see the crate README for what enabling them takes.

mod describe;
mod endtoend;
mod error;
mod insert;
mod session;

pub use describe::{
    CutOnPolygons, EdgePoint, InteriorPoint, JunctionPoint, PathEnd, PolygonPath,
};
pub use endtoend::{ConePlacement, CutUnwrap, CutUnwrapOptions, unwrap_with_cut};
pub use error::CuttingError;
pub use insert::{Inserted, VertexSource, insert};
pub use session::{DEFAULT_CROSSING_OFFSET, PreviewPolyline, Session, Start};
