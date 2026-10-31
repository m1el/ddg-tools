//! The error type of this crate.

use std::fmt;

use ddg_intrinsic::IntrinsicError;

/// Why an operation on a geodesic network failed.
///
/// Every message names the path, segment, edge or vertex concerned. A refused input changes
/// nothing. Things that are legal but worth knowing (a joint that is blocked, an iteration
/// limit that was reached) are not errors; they are in the report or certificate of the
/// operation.
#[derive(Debug, Clone, PartialEq)]
pub enum FlipError {
    /// An operation of the triangulation failed.
    Intrinsic(IntrinsicError),
    /// An element id is out of range.
    IndexOutOfRange {
        /// What kind of element ("vertex", "edge", "path", ...).
        what: &'static str,
        /// The id.
        index: usize,
        /// The exclusive upper bound it violates.
        bound: usize,
    },
    /// A path has no segment.
    EmptyPath {
        /// The path.
        path: usize,
    },
    /// Two consecutive segments of a path do not share a vertex.
    PathNotConnected {
        /// The path.
        path: usize,
        /// The position of the first of the two segments in the path.
        index: usize,
        /// The vertex where segment `index` ends.
        end: usize,
        /// The vertex where the next segment starts.
        start: usize,
    },
    /// A path marked closed does not end where it starts.
    LoopNotClosed {
        /// The path.
        path: usize,
        /// The vertex where its first segment starts.
        first: usize,
        /// The vertex where its last segment ends.
        last: usize,
    },
    /// Two consecutive vertices of a vertex sequence are not joined by an edge.
    NotAdjacent {
        /// The position of the first of the two vertices in the sequence.
        index: usize,
        /// The first vertex.
        from: usize,
        /// The second vertex.
        to: usize,
    },
    /// A vertex sequence is too short to describe a path (an open path needs two vertices, a
    /// closed one one).
    TooFewVertices {
        /// The number of vertices given.
        n: usize,
        /// Whether the path was to be closed.
        closed: bool,
    },
    /// A shortest path was asked from a vertex to itself.
    SameVertex {
        /// The vertex.
        vertex: usize,
    },
    /// No edge path joins the two vertices.
    NotConnected {
        /// The start.
        from: usize,
        /// The end.
        to: usize,
    },
    /// An explicit layering has not one list per edge.
    LayeringSize {
        /// The number of lists given.
        len: usize,
        /// The number of edges.
        n_edges: usize,
    },
    /// An explicit layering lists a segment on an edge it does not lie on.
    LayeringWrongEdge {
        /// The edge whose list names the segment.
        edge: usize,
        /// The path of the segment.
        path: usize,
        /// Its position in the path.
        index: usize,
        /// The edge the segment lies on.
        segment_edge: usize,
    },
    /// An explicit layering lists a segment twice.
    LayeringDuplicate {
        /// The edge.
        edge: usize,
        /// The path of the segment.
        path: usize,
        /// Its position in the path.
        index: usize,
    },
    /// An explicit layering does not list a segment.
    LayeringMissing {
        /// The edge the segment lies on.
        edge: usize,
        /// The path of the segment.
        path: usize,
        /// Its position in the path.
        index: usize,
    },
    /// Two joints at a vertex cross (SPEC §2.1): in the cyclic order of the segment ends at the
    /// vertex, the two pairs interleave. Inputs that cross are refused, never repaired.
    Crossing {
        /// The vertex.
        vertex: usize,
        /// The paths of the two joints (the same path twice if a path crosses itself).
        paths: [usize; 2],
    },
    /// The end of an open path is not pinned.
    UnpinnedEnd {
        /// The path.
        path: usize,
        /// The vertex.
        vertex: usize,
    },
    /// An option is out of its range.
    InvalidOptions {
        /// What is wrong.
        what: String,
    },
    /// An internal invariant does not hold: a bug, reported instead of a panic.
    Internal {
        /// What does not hold, naming the elements concerned.
        what: String,
    },
}

impl fmt::Display for FlipError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            FlipError::Intrinsic(e) => write!(f, "{e}"),
            FlipError::IndexOutOfRange { what, index, bound } => {
                write!(f, "{what} {index} is out of range (there are {bound})")
            }
            FlipError::EmptyPath { path } => write!(f, "path {path} has no segment"),
            FlipError::PathNotConnected {
                path,
                index,
                end,
                start,
            } => write!(
                f,
                "path {path} is not connected: segment {index} ends at vertex {end}, but \
                 segment {} starts at vertex {start}",
                index + 1
            ),
            FlipError::LoopNotClosed { path, first, last } => write!(
                f,
                "path {path} is marked closed, but it starts at vertex {first} and ends at \
                 vertex {last}"
            ),
            FlipError::NotAdjacent { index, from, to } => write!(
                f,
                "vertices {index} and {} of the sequence ({from} and {to}) are not joined by \
                 an edge",
                index + 1
            ),
            FlipError::TooFewVertices { n, closed } => {
                if *closed {
                    write!(f, "a closed path needs at least one vertex, got {n}")
                } else {
                    write!(f, "an open path needs at least two vertices, got {n}")
                }
            }
            FlipError::SameVertex { vertex } => write!(
                f,
                "a shortest path from vertex {vertex} to itself was asked for; it is empty"
            ),
            FlipError::NotConnected { from, to } => write!(
                f,
                "vertices {from} and {to} are not connected by edges (different components)"
            ),
            FlipError::LayeringSize { len, n_edges } => write!(
                f,
                "the layering has {len} lists, but the triangulation has {n_edges} edges"
            ),
            FlipError::LayeringWrongEdge {
                edge,
                path,
                index,
                segment_edge,
            } => write!(
                f,
                "the layering of edge {edge} lists segment {index} of path {path}, which lies \
                 on edge {segment_edge}"
            ),
            FlipError::LayeringDuplicate { edge, path, index } => write!(
                f,
                "the layering of edge {edge} lists segment {index} of path {path} more than once"
            ),
            FlipError::LayeringMissing { edge, path, index } => write!(
                f,
                "the layering of edge {edge} does not list segment {index} of path {path}, \
                 which lies on it"
            ),
            FlipError::Crossing { vertex, paths } => {
                if paths[0] == paths[1] {
                    write!(f, "path {} crosses itself at vertex {vertex}", paths[0])
                } else {
                    write!(
                        f,
                        "paths {} and {} cross at vertex {vertex}",
                        paths[0], paths[1]
                    )
                }
            }
            FlipError::UnpinnedEnd { path, vertex } => write!(
                f,
                "path {path} is open and ends at vertex {vertex}, which is not pinned"
            ),
            FlipError::InvalidOptions { what } => write!(f, "invalid options: {what}"),
            FlipError::Internal { what } => write!(f, "internal error (a bug): {what}"),
        }
    }
}

impl std::error::Error for FlipError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            FlipError::Intrinsic(e) => Some(e),
            _ => None,
        }
    }
}

impl From<IntrinsicError> for FlipError {
    fn from(e: IntrinsicError) -> Self {
        FlipError::Intrinsic(e)
    }
}
