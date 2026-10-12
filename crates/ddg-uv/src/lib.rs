//! End-to-end UV unwrapping (work in progress): polygon mesh in, one UV per polygon corner out.
//!
//! Built on `ddg-optimal-cones` (Soliman, Slepčev, Crane 2018) and `ddg-bff` (Sawhney, Crane
//! 2017). The stages that exist so far are validation and seam cutting ([`prepare`]), the
//! chart packer ([`pack`]) and the distortion measure ([`measure`]).

pub mod error;
pub mod input;
pub mod measure;
pub mod options;
pub mod pack;
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
pub use pack::{PackReport, pack_charts};
pub use prepare::{Piece, Prepared, prepare};
