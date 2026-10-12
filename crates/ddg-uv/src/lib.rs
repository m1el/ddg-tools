//! End-to-end UV unwrapping (work in progress): polygon mesh in, one UV per polygon corner out.
//!
//! Built on `ddg-optimal-cones` (Soliman, Slepčev, Crane 2018) and `ddg-bff` (Sawhney, Crane
//! 2017). The stages that exist so far are validation and seam cutting ([`prepare`]), the
//! chart packer ([`pack`]) and the distortion measure ([`measure`]).

pub mod cones;
pub mod error;
pub mod flatten;
pub mod input;
pub mod measure;
pub mod options;
pub mod pack;
pub mod pipeline;
pub mod prepare;
pub mod report;
mod slit;

pub use ddg_core::io::Cone;
pub use error::{ErrorKind, UvError};
pub use input::{MeshInput, PolygonMesh, fan_triangulation};
pub use options::{
    ClosedSurfaces, ConeSigns, Degenerate, MergeClusters, Orientation, PackOptions,
    SlitPlacement, UnwrapOptions,
};
pub use cones::{ConeStage, PieceCones, solve_cones};
pub use flatten::{FlatStage, PieceFlat, flatten};
pub use measure::{AreaDistortion, area_distortion};
pub use pack::{PackReport, pack_charts};
pub use pipeline::{UvMap, assemble, unwrap};
pub use prepare::{Piece, Prepared, prepare};
pub use report::{
    ChartReport, ClosedHandling, ConeSolveReport, ConeStats, DegenerateReport, FlattenReport,
    MergeOutcome, Report, SeamReport, SteeringReport, Timings, ValidationReport, Value,
};
