//! Errors of the cut session and of the insertion of a cut into a polygon mesh.

use std::fmt;

use ddg_variational_cuts::CutError;

use crate::error::UvError;

/// Why a cut session, a cut description or an insertion was refused. Every message names the
/// vertex, corner, edge or polygon concerned, in the caller's numbering.
#[derive(Debug, Clone, PartialEq)]
pub enum CuttingError {
    /// The polygon mesh, its triangulation or its seams were refused (see [`UvError`]).
    Mesh(UvError),
    /// The cut flow refused something (see [`CutError`]); its messages name vertices and faces
    /// of the triangle mesh the flow runs on, which [`crate::cutting::Session`] numbers (use the
    /// maps of the session to translate).
    Flow(CutError),
    /// Fewer than two regions.
    TooFewRegions {
        /// The number asked for.
        got: usize,
    },
    /// More than two regions: the cut flow handles two regions only for now (SPEC §8 of
    /// `ddg-variational-cuts` is not implemented).
    RegionsNotAvailable {
        /// The number asked for.
        n_regions: usize,
    },
    /// An array has the wrong number of entries.
    WrongLength {
        /// Which array.
        what: &'static str,
        /// What it has one entry for.
        per: &'static str,
        /// The required length.
        expected: usize,
        /// The length that was passed.
        got: usize,
    },
    /// A value is NaN or infinite.
    NonFinite {
        /// Which array.
        what: &'static str,
        /// The region (or 0).
        region: usize,
        /// The polygon corner (or vertex) where it is.
        position: usize,
        /// The value.
        value: f64,
    },
    /// Two corners of one vertex that no seam separates carry different values of a region
    /// function.
    CornersDisagree {
        /// The vertex.
        vertex: u32,
        /// The region function.
        region: usize,
        /// The two corners.
        corners: [u32; 2],
        /// Their values.
        values: [f64; 2],
    },
    /// A polygon label is not below the number of regions.
    LabelOutOfRange {
        /// The polygon.
        polygon: u32,
        /// Its label.
        label: u32,
        /// The number of regions.
        n_regions: usize,
    },
    /// A vertex index is out of range.
    VertexOutOfRange {
        /// What the index was for.
        what: &'static str,
        /// The index.
        vertex: u32,
        /// The number of vertices.
        n_vertices: usize,
    },
    /// An option or argument is outside its range.
    InvalidOption {
        /// The name.
        name: &'static str,
        /// What is wrong.
        reason: String,
    },
    /// A cut description does not fit the polygon mesh it is inserted into.
    InvalidDescription {
        /// The polygon, if the problem is local to one.
        polygon: Option<u32>,
        /// What is wrong.
        problem: String,
    },
    /// The description uses something that is not implemented (a path that ends at a junction
    /// of three or more regions).
    Unsupported {
        /// What.
        what: String,
    },
    /// A check on a result failed. This is a bug.
    Internal {
        /// What failed.
        message: String,
    },
}

impl From<UvError> for CuttingError {
    fn from(e: UvError) -> Self {
        CuttingError::Mesh(e)
    }
}

impl From<CutError> for CuttingError {
    fn from(e: CutError) -> Self {
        CuttingError::Flow(e)
    }
}

pub(crate) fn internal(message: impl Into<String>) -> CuttingError {
    CuttingError::Internal {
        message: message.into(),
    }
}

pub(crate) fn invalid(polygon: Option<u32>, problem: impl Into<String>) -> CuttingError {
    CuttingError::InvalidDescription {
        polygon,
        problem: problem.into(),
    }
}

impl fmt::Display for CuttingError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        use CuttingError::*;
        match self {
            Mesh(e) => write!(f, "{e}"),
            Flow(e) => write!(f, "cut flow: {e}"),
            TooFewRegions { got } => write!(f, "{got} regions asked for; at least 2 are needed"),
            RegionsNotAvailable { n_regions } => write!(
                f,
                "{n_regions} regions asked for; the cut flow handles exactly 2 regions so far \
                 (N regions are not available yet)"
            ),
            WrongLength {
                what,
                per,
                expected,
                got,
            } => write!(
                f,
                "{what} has {got} entries but needs one per {per} ({expected})"
            ),
            NonFinite {
                what,
                region,
                position,
                value,
            } => write!(f, "{what} {region} is {value} at entry {position}"),
            CornersDisagree {
                vertex,
                region,
                corners: [c0, c1],
                values: [v0, v1],
            } => write!(
                f,
                "region function {region} differs at vertex {vertex}: corner {c0} has {v0} and \
                 corner {c1} has {v1}, but no seam separates them"
            ),
            LabelOutOfRange {
                polygon,
                label,
                n_regions,
            } => write!(
                f,
                "polygon {polygon} has the label {label}, but there are only {n_regions} regions"
            ),
            VertexOutOfRange {
                what,
                vertex,
                n_vertices,
            } => write!(f, "{what} {vertex} is out of range ({n_vertices} vertices)"),
            InvalidOption { name, reason } => write!(f, "{name}: {reason}"),
            InvalidDescription { polygon, problem } => match polygon {
                Some(p) => write!(f, "cut description, polygon {p}: {problem}"),
                None => write!(f, "cut description: {problem}"),
            },
            Unsupported { what } => write!(f, "not available yet: {what}"),
            Internal { message } => write!(f, "internal error: {message}"),
        }
    }
}

impl std::error::Error for CuttingError {}
