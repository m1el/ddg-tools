//! End-to-end UV unwrapping: a polygon mesh goes in, one UV per polygon corner comes out.
//!
//! The pipeline ties the workspace together. It places cone singularities with
//!
//! > Yousuf Soliman, Dejan Slepčev, Keenan Crane.
//! > *Optimal Cone Singularities for Conformal Flattening.*
//! > ACM Transactions on Graphics 37(4), 2018. <https://doi.org/10.1145/3197517.3201367>
//!
//! (crate `ddg-optimal-cones`) and flattens the surface conformally through those cones with
//!
//! > Rohan Sawhney, Keenan Crane. *Boundary First Flattening.*
//! > ACM Transactions on Graphics 37(1), 2017. <https://doi.org/10.1145/3132705>
//!
//! (crate `ddg-bff`, which also plans and makes the cuts). What this crate adds is everything
//! around the two: validation with messages in the caller's numbering, polygons with a
//! supplied triangulation, seams, steering, the policies for the cases the papers leave open
//! (negative cones, closed surfaces, cone clusters), chart packing, and a report.
//!
//! # Example
//!
//! ```
//! use ddg_uv::primitives::grid_cube;
//! use ddg_uv::{MeshInput, UnwrapOptions, unwrap};
//!
//! # fn main() -> Result<(), Box<dyn std::error::Error>> {
//! // A cube of 4 x 4 quads per side. Any polygon mesh will do: positions, the vertex of every
//! // polygon corner, the corner offsets of the polygons, and a triangulation of the polygons
//! // given as triples of corner indices.
//! let mesh = grid_cube(4, 2.0);
//! // Seams are vertex pairs; here the first edge of the first quad.
//! let quad = mesh.polygon(0);
//! let seams = [[quad[0], quad[1]]];
//! let input = MeshInput { seam_edges: &seams, ..mesh.input() };
//!
//! let map = unwrap(&input, &UnwrapOptions::default())?;
//!
//! assert_eq!(map.uvs.len(), mesh.corner_vertices.len()); // one UV per polygon corner
//! assert!(map.uvs.iter().all(|uv| uv.iter().all(|x| (0.0..=1.0).contains(x)))); // packed
//! assert_eq!(map.polygon_charts.len(), mesh.n_polygons());
//! assert!(map.cut_edges.contains(&[quad[0].min(quad[1]), quad[0].max(quad[1])]));
//! for cone in &map.cones {
//!     println!("vertex {} carries a cone of {} rad", cone.vertex, cone.angle);
//! }
//! for warning in &map.report.warnings {
//!     println!("warning: {warning}");
//! }
//! # Ok(())
//! # }
//! ```
//!
//! The module [`cutting`] puts the cut flow of "Variational Surface Cutting" (Sharp, Crane 2018,
//! crate `ddg-variational-cuts`) in front of this pipeline for polygon meshes: a session that
//! moves a cut on the seam-cut triangulation, the cut expressed on the caller's polygons, the
//! insertion of the cut into a polygon mesh, and `unwrap` with the cut as seams.
//!
//! # Stages
//!
//! [`unwrap`] runs these in order; [`prepare`], [`solve_cones`], [`flatten`] and [`assemble`]
//! are the same stages as separate functions (the benchmarks time them apart).
//!
//! 1. **Validation** ([`prepare`]). The arrays are checked ([`MeshInput::validate`]). Vertices
//!    that no polygon uses are ignored. The triangulated surface must be a manifold: every
//!    refusal names the vertex, edge and polygons in the caller's numbering ([`UvError`]).
//!    Inconsistent winding is repaired per connected component with a warning
//!    ([`Orientation`]); a surface that cannot be oriented is an error. Triangles without
//!    area are handled as [`Degenerate`] says, by default by mollifying the edge lengths,
//!    and reported.
//! 2. **Seams** ([`prepare`]). The surface is cut along the marked seams before anything else,
//!    so that the cone solver and the flattener both work on the cut surface: the seams are
//!    boundary, the scale factor is zero along them, and the two sides of a seam keep their
//!    length. The connected pieces that result are the charts. A seam that is not an edge of a
//!    polygon is ignored with a warning; a seam of a single edge between two interior vertices
//!    is extended by one edge, because a slit of one edge has no two sides.
//! 3. **Cones, then flattening**, per piece ([`solve_cones`], [`flatten`]). First the spike
//!    check ([`spikes`], [`SpikeCheck`]): if a handful of vertices carry almost all of the
//!    distortion the piece would have without cones (thin spikes, which wreck the result), the
//!    report names them and warns. The cone problem is solved with the steering translated into its terms (pinned vertices are free cone
//!    sites, prescribed cones fixed angles, forbidden vertices excluded) and with the policies
//!    of [`UnwrapOptions`]: [`ConeSigns`], [`ClosedSurfaces`], [`MergeClusters`]. The piece is
//!    then cut through its cones **along polygon edges only** and flattened. Because cuts never
//!    run along a diagonal of the triangulation, all triangle corners that belong to one
//!    polygon corner receive the same UV; this is checked, and a violation is
//!    [`UvError::Internal`].
//! 4. **Packing** ([`assemble`], [`pack`]). All charts go into the unit square with one common
//!    scale and a margin ([`PackOptions`]).
//!
//! Steering applies at vertices that are interior to their piece. On a boundary, a seam or a
//! slit there are no cones: pinned and forbidden vertices have no effect there, and a
//! prescribed cone that cannot be imposed is listed in the report with a warning.
//!
//! # What comes back
//!
//! A [`UvMap`]: the UV of every polygon corner, the chart of every polygon, the cones, the
//! edges that were cut, and a [`Report`] with timings per stage, cone statistics, the
//! solver's certificate summary, distortion and overlap per chart, and warnings. Nothing is
//! changed, clamped or repaired without a line in the report. The result is deterministic,
//! nothing panics, no UV is NaN, and the input is only read.
//!
//! # Conventions
//!
//! Those of the repository's `CLAUDE.md`: a cone angle is `µ = 2π − (angle sum around the
//! vertex in the plane)`, positive for a cone tip; `λ` is taken with the area of every piece
//! normalized to 1. `ddg-core` has no triangulator; that is why the triangulation is an input.
//!
//! The crate README holds the measurements behind every default.

#![warn(missing_docs)]

pub mod cones;
pub mod cutting;
pub mod error;
pub mod flatten;
pub mod input;
pub mod measure;
pub mod options;
pub mod pack;
pub mod pipeline;
pub mod prepare;
pub mod primitives;
pub mod report;
mod slit;
pub mod spikes;

pub use ddg_core::io::Cone;

pub use cones::{ConeStage, PieceCones, solve_cones};
pub use error::{ErrorKind, UvError};
pub use flatten::{FlatStage, PieceFlat, flatten};
pub use input::{MeshInput, PolygonMesh, fan_triangulation};
pub use measure::{AreaDistortion, area_distortion};
pub use options::{
    BoundaryValues, ClosedSurfaces, ConeSigns, Degenerate, MergeClusters, Orientation, PackOptions,
    SlitPlacement, SpikeAction, SpikeCheck, UnwrapOptions,
};
pub use pack::{PackReport, pack_charts};
pub use pipeline::{UvMap, assemble, unwrap};
pub use prepare::{Piece, Prepared, prepare};
pub use report::{
    ChartReport, ClosedHandling, ConeSolveReport, ConeStats, DegenerateReport, FlattenReport,
    MergeOutcome, Report, SeamLengths, SeamReport, SteeringReport, Timings, ValidationReport,
    Value,
};
pub use spikes::SpikeReport;
