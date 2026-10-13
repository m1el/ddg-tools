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
//! # Example
//!
//! ```
//! use ddg_uv::cutting::{ConePlacement, CutUnwrapOptions, DEFAULT_CROSSING_OFFSET, Session, Start};
//! use ddg_uv::primitives::uv_sphere;
//! use ddg_variational_cuts::{CutProblem, Energy, FlowOptions, StepSize};
//!
//! # fn main() -> Result<(), Box<dyn std::error::Error>> {
//! let mesh = uv_sphere(24, 12, 1.0).expect("valid parameters");
//! // Two regions, started as the geodesic disk of radius 0.9 around vertex 5.
//! let start = Start::GeodesicCircle { vertex: 5, radius: 0.9 };
//! let mut session = Session::new(&mesh.input(), 2, &start)?;
//!
//! let problem = CutProblem::new(Energy::Hencky, 0.05);
//! let options = FlowOptions::new(StepSize::default());
//! for _ in 0..3 {
//!     session.step(&problem, &options)?;
//! }
//!
//! // The state: one function per region and per polygon corner.
//! let state = session.state();
//! assert_eq!((state.len(), state[0].len()), (2, mesh.corner_vertices.len()));
//! // The cut on the caller's polygons, and the polygons split along it.
//! let cut = session.describe_cut(DEFAULT_CROSSING_OFFSET)?;
//! let inserted = ddg_uv::cutting::insert(&mesh, &cut)?;
//! assert!(inserted.mesh.n_polygons() > mesh.n_polygons());
//! // The UVs: insert, then unwrap with the cut as seams.
//! let (_, result) = session.unwrap(
//!     DEFAULT_CROSSING_OFFSET,
//!     &CutUnwrapOptions { cones: ConePlacement::Off, ..Default::default() },
//! )?;
//! assert_eq!(result.map.uvs.len(), result.inserted.mesh.corner_vertices.len());
//! # Ok(())
//! # }
//! ```
//!
//! # The pieces
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
//! **Several regions.** With more than two regions the session runs the core crate's
//! `RegionFlow` (so does `FlowKind::Regions` for two). The state is one function per region and
//! corner; the cut is a network: where three regions meet in a triangle of the triangulation the
//! description has a `JunctionPoint`, and three paths end there. `insert` builds a vertex
//! inside that triangle and splits it into three quadrilaterals (the exact construction is in
//! the documentation of `insert`). `Session::recommended_options` returns the options the
//! core crate recommends for the number of regions (`FlowOptions::for_regions` for more than
//! two). Every region component then becomes a chart in `unwrap_with_cut`.

mod describe;
mod endtoend;
mod error;
mod insert;
mod session;

pub use describe::{CutOnPolygons, EdgePoint, InteriorPoint, JunctionPoint, PathEnd, PolygonPath};
pub use endtoend::{ConePlacement, CutUnwrap, CutUnwrapOptions, unwrap_with_cut};
pub use error::CuttingError;
pub use insert::{Inserted, VertexSource, insert};
pub use session::{DEFAULT_CROSSING_OFFSET, FlowKind, PreviewPolyline, Session, Start};
