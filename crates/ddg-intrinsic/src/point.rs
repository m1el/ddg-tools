//! Points of a triangulated surface.

use std::fmt;

/// A point of a triangulated surface: a vertex, a point of an edge, or a point of a face.
///
/// The ids refer to whichever triangulation the point is used with: the input mesh
/// ([`InputSurface`](crate::InputSurface)) or the intrinsic triangulation
/// ([`IntrinsicTriangulation`](crate::IntrinsicTriangulation)). Each kind of point has a
/// *frame*, the direction that angle zero refers to in tracing queries and in signposts
/// (paper Sec. 2.3); angles grow counter-clockwise.
///
/// * [`SurfacePoint::Vertex`]: the frame of the vertex, in normalized angles (see
///   [`IntrinsicTriangulation::signpost`](crate::IntrinsicTriangulation::signpost)).
/// * [`SurfacePoint::Edge`]: angle zero is the direction of the edge. On the input mesh an
///   edge runs from its lower-numbered to its higher-numbered vertex; on the intrinsic
///   triangulation from the tail to the head of
///   [`edge_halfedge`](crate::IntrinsicTriangulation::edge_halfedge). `t` is measured the same
///   way.
/// * [`SurfacePoint::Face`]: angle zero is the direction of halfedge `3 face`, from
///   `faces[face][0]` to `faces[face][1]`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SurfacePoint {
    /// A vertex.
    Vertex(u32),
    /// The point at parameter `t ∈ [0, 1]` of an edge: `(1 − t) · tail + t · head`.
    Edge {
        /// The edge.
        edge: u32,
        /// The parameter along the edge, from its tail.
        t: f64,
    },
    /// A point of a face in barycentric coordinates: `bary[k]` is the weight of the vertex
    /// `faces[face][k]`. The coordinates are non-negative and sum to one.
    Face {
        /// The face.
        face: u32,
        /// The barycentric coordinates.
        bary: [f64; 3],
    },
}

impl SurfacePoint {
    /// The vertex, if the point is one.
    pub fn as_vertex(&self) -> Option<u32> {
        match *self {
            SurfacePoint::Vertex(v) => Some(v),
            _ => None,
        }
    }
}

impl fmt::Display for SurfacePoint {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            SurfacePoint::Vertex(v) => write!(f, "vertex {v}"),
            SurfacePoint::Edge { edge, t } => write!(f, "the point t = {t:?} of edge {edge}"),
            SurfacePoint::Face {
                face,
                bary: [a, b, c],
            } => write!(f, "the point ({a:?}, {b:?}, {c:?}) of face {face}"),
        }
    }
}

/// How far the sum of barycentric coordinates may be from one for a point to be accepted.
pub(crate) const BARY_SUM_TOLERANCE: f64 = 1e-9;

/// What is wrong with `point` on a surface with the given element counts, if anything.
pub(crate) fn point_defect(
    point: &SurfacePoint,
    n_vertices: usize,
    n_edges: usize,
    n_faces: usize,
) -> Option<&'static str> {
    match *point {
        SurfacePoint::Vertex(v) => {
            (v as usize >= n_vertices).then_some("the vertex does not exist")
        }
        SurfacePoint::Edge { edge, t } => {
            if edge as usize >= n_edges {
                Some("the edge does not exist")
            } else if !(t.is_finite() && (0.0..=1.0).contains(&t)) {
                Some("the edge parameter is not in [0, 1]")
            } else {
                None
            }
        }
        SurfacePoint::Face { face, bary } => {
            if face as usize >= n_faces {
                Some("the face does not exist")
            } else if !bary.iter().all(|b| b.is_finite() && *b >= 0.0) {
                Some("a barycentric coordinate is negative or not finite")
            } else if ((bary[0] + bary[1] + bary[2]) - 1.0).abs() > BARY_SUM_TOLERANCE {
                Some("the barycentric coordinates do not sum to one")
            } else {
                None
            }
        }
    }
}
