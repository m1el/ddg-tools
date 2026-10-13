//! The error type of this crate, and the reasons an operation can be refused.

use std::fmt;

use ddg_core::{GeometryError, MeshError};

use crate::point::SurfacePoint;

/// Why an operation on an intrinsic triangulation failed.
///
/// Every message names the vertex, edge or face concerned. Things that are legal but worth
/// knowing (a flip limit that was reached, an edge that could not be flipped during Delaunay
/// flipping) are not errors; they are in the report of the operation.
#[derive(Debug, Clone, PartialEq)]
pub enum IntrinsicError {
    /// The geometry was not built for the mesh.
    Geometry(GeometryError),
    /// The triangulation cannot be handed out as a `ddg_core::TriMesh` (see
    /// [`IntrinsicTriangulation::to_tri_mesh`](crate::IntrinsicTriangulation::to_tri_mesh)).
    NotSimplicial {
        /// What stands in the way, naming the first offending element.
        reason: String,
        /// The error of `TriMesh::new`, if that is what refused the faces.
        mesh_error: Option<MeshError>,
    },
    /// An element id is out of range.
    IndexOutOfRange {
        /// What kind of element ("edge", "vertex", ...).
        what: &'static str,
        /// The id.
        index: usize,
        /// The exclusive upper bound it violates.
        bound: usize,
    },
    /// An edge flip was refused.
    FlipRefused {
        /// The edge.
        edge: u32,
        /// Its end points (tail and head of its halfedge).
        vertices: [u32; 2],
        /// Why.
        reason: FlipRefusal,
    },
    /// A point passed to a query is not a point of the surface it refers to.
    InvalidPoint {
        /// The point.
        point: SurfacePoint,
        /// What is wrong with it.
        reason: &'static str,
    },
    /// The direction of a tracing query is not a finite angle in `[0, 2π)`.
    InvalidDirection {
        /// The offending angle.
        direction: f64,
    },
    /// The length of a tracing query is not positive and finite.
    InvalidLength {
        /// The offending length.
        length: f64,
    },
    /// A tracing query from a point on the boundary of a face was given a direction that does
    /// not lead into that face.
    DirectionLeavesFace {
        /// The start point.
        point: SurfacePoint,
        /// The direction, in the frame of the face.
        direction: f64,
    },
    /// A tracing query crossed more edges than its limit allows.
    CrossingLimit {
        /// Where the query started.
        start: SurfacePoint,
        /// The limit that was reached.
        limit: usize,
        /// The length travelled when it was reached.
        travelled: f64,
    },
    /// The trace of one intrinsic edge over the input mesh failed its check.
    EdgeTrace(EdgeTraceFailure),
    /// The traces of one or more intrinsic edges over the input mesh failed their check. All of
    /// them are listed, in ascending order of the edge.
    EdgeTraces {
        /// Number of edges of the triangulation.
        n_edges: usize,
        /// The failures.
        failures: Vec<EdgeTraceFailure>,
    },
    /// An option is out of range.
    InvalidOptions {
        /// Which option and why.
        reason: String,
    },
    /// A face has edge lengths that do not form a triangle of positive area. Every operation
    /// keeps the triangle inequalities, so this cannot happen for a triangulation built by this
    /// crate; it is reported instead of panicking.
    DegenerateFace {
        /// Which surface: `"input mesh"` or `"intrinsic triangulation"`.
        surface: &'static str,
        /// The face.
        face: u32,
        /// The lengths of its three halfedges.
        lengths: [f64; 3],
    },
    /// The data of the triangulation is inconsistent (found by
    /// [`IntrinsicTriangulation::check`](crate::IntrinsicTriangulation::check)). This is a bug
    /// of this crate, reported instead of panicking.
    Inconsistent {
        /// What is wrong and where.
        what: String,
    },
}

/// Why an edge cannot be flipped (paper Sec. 3.3.1).
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum FlipRefusal {
    /// The edge lies on the boundary: it has one triangle only.
    BoundaryEdge,
    /// Both sides of the edge are the same triangle, which happens exactly when one of its end
    /// points has degree one (paper Fig. 6, right).
    DegreeOneEndpoint {
        /// The end point of degree one.
        vertex: u32,
    },
    /// The two triangles of the edge, unfolded into the plane, do not form a strictly convex
    /// quadrilateral: the other diagonal does not cross the edge in its interior.
    NotConvex {
        /// Twice the signed area of the triangle the flip would create at the edge's tail.
        area_at_tail: f64,
        /// Twice the signed area of the triangle the flip would create at the edge's head.
        area_at_head: f64,
    },
    /// The quadrilateral is strictly convex, but so barely that a triangle the flip would
    /// create is too thin to be described by its edge lengths: its shape number is below
    /// [`FLIP_SHAPE_LIMIT`](crate::formulas::FLIP_SHAPE_LIMIT), or its lengths fail the
    /// triangle inequality by rounding.
    TooThin {
        /// The length the new edge would have.
        new_length: f64,
        /// The smaller shape number of the two triangles the flip would create: twice the
        /// area over the product of the new edge and the triangle's longest side.
        shape: f64,
    },
}

impl fmt::Display for FlipRefusal {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            FlipRefusal::BoundaryEdge => write!(f, "it is a boundary edge"),
            FlipRefusal::DegreeOneEndpoint { vertex } => write!(
                f,
                "its end point {vertex} has degree one (both sides of the edge are the same \
                 triangle)"
            ),
            FlipRefusal::NotConvex {
                area_at_tail,
                area_at_head,
            } => write!(
                f,
                "its two triangles do not form a strictly convex quadrilateral (twice the signed \
                 areas of the triangles a flip would create: {area_at_tail} and {area_at_head})"
            ),
            FlipRefusal::TooThin { new_length, shape } => write!(
                f,
                "the flip would create a triangle too thin to be described by its edge lengths \
                 (new edge of length {new_length}; twice the area is {shape:e} times the product \
                 of the new edge and the longest side, the limit is {:e})",
                crate::formulas::FLIP_SHAPE_LIMIT
            ),
        }
    }
}

/// The trace of an intrinsic edge over the input mesh that did not arrive at the edge's far
/// vertex.
#[derive(Debug, Clone, PartialEq)]
pub struct EdgeTraceFailure {
    /// The intrinsic edge.
    pub edge: u32,
    /// Its end points: the vertex the trace started from and the one it should have reached.
    pub vertices: [u32; 2],
    /// Its length.
    pub length: f64,
    /// The largest mismatch that would have been accepted (an absolute length).
    pub tolerance: f64,
    /// What happened.
    pub reason: EdgeTraceFailureReason,
}

/// How the trace of an intrinsic edge failed.
#[derive(Debug, Clone, PartialEq)]
pub enum EdgeTraceFailureReason {
    /// The path was followed for the edge's whole length without meeting the far vertex.
    Missed {
        /// Where the path ended on the input mesh.
        end: SurfacePoint,
        /// The smallest distance between the end of the path and the far vertex, over all
        /// places where the far vertex's place on the input mesh lay in a face the path went
        /// through (measured in the unfolding of those faces). `None` if it never did.
        mismatch: Option<f64>,
        /// Number of input edges crossed.
        n_crossings: usize,
    },
    /// The path left the surface through a boundary edge of the input mesh.
    HitBoundary {
        /// The point of the boundary edge.
        at: SurfacePoint,
        /// The length travelled until then.
        travelled: f64,
    },
    /// The path ran along an input edge into an input vertex that is not the far vertex.
    HitVertex {
        /// The input vertex.
        vertex: u32,
        /// The length travelled until then.
        travelled: f64,
    },
    /// The path crossed more input edges than the limit allows.
    CrossingLimit {
        /// The limit.
        limit: usize,
        /// The length travelled until then.
        travelled: f64,
    },
    /// The tracing query itself was refused.
    Query(Box<IntrinsicError>),
}

impl fmt::Display for EdgeTraceFailure {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let [a, b] = self.vertices;
        write!(
            f,
            "the trace of intrinsic edge {} = ({a}, {b}) of length {} over the input mesh ",
            self.edge, self.length
        )?;
        match &self.reason {
            EdgeTraceFailureReason::Missed {
                end,
                mismatch,
                n_crossings,
            } => {
                write!(
                    f,
                    "ended at {end} after {n_crossings} crossings without reaching vertex {b}: "
                )?;
                match mismatch {
                    Some(d) => write!(
                        f,
                        "the closest it came is {d}, the tolerance is {}",
                        self.tolerance
                    ),
                    None => write!(
                        f,
                        "the place of vertex {b} is in none of the input faces it went through"
                    ),
                }
            }
            EdgeTraceFailureReason::HitBoundary { at, travelled } => {
                write!(f, "left the surface at {at} after a length of {travelled}")
            }
            EdgeTraceFailureReason::HitVertex { vertex, travelled } => write!(
                f,
                "ran into input vertex {vertex} after a length of {travelled}"
            ),
            EdgeTraceFailureReason::CrossingLimit { limit, travelled } => write!(
                f,
                "crossed {limit} input edges (the limit) within a length of {travelled}"
            ),
            EdgeTraceFailureReason::Query(e) => write!(f, "was refused: {e}"),
        }
    }
}

impl fmt::Display for IntrinsicError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            IntrinsicError::Geometry(e) => write!(f, "{e}"),
            IntrinsicError::NotSimplicial { reason, mesh_error } => {
                write!(
                    f,
                    "the intrinsic triangulation is not a simplicial complex: {reason}"
                )?;
                if let Some(e) = mesh_error {
                    write!(f, " ({e})")?;
                }
                Ok(())
            }
            IntrinsicError::IndexOutOfRange { what, index, bound } => {
                write!(f, "{what} {index} is out of range (must be < {bound})")
            }
            IntrinsicError::FlipRefused {
                edge,
                vertices: [a, b],
                reason,
            } => write!(f, "edge {edge} = ({a}, {b}) cannot be flipped: {reason}"),
            IntrinsicError::InvalidPoint { point, reason } => {
                write!(f, "{point} is not a valid point: {reason}")
            }
            IntrinsicError::InvalidDirection { direction } => write!(
                f,
                "the direction {direction} of a tracing query is not a finite angle in [0, 2π)"
            ),
            IntrinsicError::InvalidLength { length } => write!(
                f,
                "the length {length} of a tracing query is not positive and finite"
            ),
            IntrinsicError::DirectionLeavesFace { point, direction } => write!(
                f,
                "the direction {direction} does not lead from {point} into its face"
            ),
            IntrinsicError::CrossingLimit {
                start,
                limit,
                travelled,
            } => write!(
                f,
                "the tracing query from {start} crossed {limit} edges (the limit) within a \
                 length of {travelled}"
            ),
            IntrinsicError::EdgeTrace(failure) => write!(f, "{failure}"),
            IntrinsicError::EdgeTraces { n_edges, failures } => {
                const MAX: usize = 4;
                write!(
                    f,
                    "{} of {n_edges} intrinsic edges could not be traced over the input mesh: ",
                    failures.len()
                )?;
                for (i, failure) in failures.iter().take(MAX).enumerate() {
                    if i > 0 {
                        write!(f, "; ")?;
                    }
                    write!(f, "{failure}")?;
                }
                if failures.len() > MAX {
                    write!(f, "; ... ({} more)", failures.len() - MAX)?;
                }
                Ok(())
            }
            IntrinsicError::InvalidOptions { reason } => write!(f, "invalid options: {reason}"),
            IntrinsicError::DegenerateFace {
                surface,
                face,
                lengths: [a, b, c],
            } => write!(
                f,
                "face {face} of the {surface} has edge lengths {a}, {b}, {c}, which do not form \
                 a triangle of positive area"
            ),
            IntrinsicError::Inconsistent { what } => {
                write!(f, "inconsistent intrinsic triangulation: {what}")
            }
        }
    }
}

impl std::error::Error for IntrinsicError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            IntrinsicError::Geometry(e) => Some(e),
            IntrinsicError::NotSimplicial {
                mesh_error: Some(e),
                ..
            } => Some(e),
            _ => None,
        }
    }
}

impl From<GeometryError> for IntrinsicError {
    fn from(e: GeometryError) -> Self {
        IntrinsicError::Geometry(e)
    }
}
