//! The shared interface of intrinsic triangulations ([`Triangulation`], `ddg_intrinsic` SPEC
//! N50) for the integer representation, so that code written against it (Delaunay flipping,
//! FlipOut, the applications) runs on integer coordinates.
//!
//! Flips, edge splits and the constrained mask are this crate's operations
//! ([`IntegerTriangulation::flip_edge`], [`IntegerTriangulation::split_edge`]); Delaunay
//! flipping is the trait's provided [`Triangulation::flip_to_delaunay`], the same function the
//! signpost representation runs, so both make the same flips in the same order on the same input
//! (SPEC §4.3). Errors of this crate are passed through as
//! [`IntrinsicError::Representation`] with `representation: "integer coordinates"`, refusals
//! of the shared kind (flips, splits, removals) as themselves, and operations this crate does
//! not provide yet as [`IntrinsicError::Unsupported`]: insertion in a face (work package WP4),
//! splits of edges that cross input edges and the path of an edge over the input mesh (both
//! need the geometric crossings of WP2).

use ddg_intrinsic::{
    Complex, EdgeTrace, EdgeTraceOptions, Flip, FlipRefusal, InputSurface, Insertion,
    IntrinsicError, SurfacePoint, Triangulation,
};

use crate::coords::IntegerTriangulation;
use crate::error::IntegerError;

/// The name of this representation in [`IntrinsicError::Representation`] and
/// [`IntrinsicError::Unsupported`].
pub const REPRESENTATION: &str = "integer coordinates";

/// An error of this crate as an error of the shared interface: refusals and index errors of the
/// shared code as themselves, [`IntegerError::Unsupported`] as
/// [`IntrinsicError::Unsupported`], everything else as [`IntrinsicError::Representation`]
/// with its message.
pub fn to_intrinsic_error(error: IntegerError, operation: &'static str) -> IntrinsicError {
    match error {
        IntegerError::Intrinsic(e) => e,
        IntegerError::IndexOutOfRange { what, index, bound } => {
            IntrinsicError::IndexOutOfRange { what, index, bound }
        }
        IntegerError::Unsupported { operation, reason } => IntrinsicError::Unsupported {
            representation: REPRESENTATION,
            operation,
            reason,
        },
        other => IntrinsicError::Representation {
            representation: REPRESENTATION,
            operation,
            message: other.to_string(),
        },
    }
}

impl Triangulation for IntegerTriangulation {
    #[inline]
    fn complex(&self) -> &Complex {
        &self.complex
    }

    #[inline]
    fn input(&self) -> &InputSurface {
        &self.input
    }

    #[inline]
    fn edge_lengths(&self) -> &[f64] {
        &self.lengths
    }

    #[inline]
    fn vertex_angle_sum(&self, v: usize) -> f64 {
        IntegerTriangulation::vertex_angle_sum(self, v)
    }

    #[inline]
    fn vertex_location(&self, v: usize) -> SurfacePoint {
        self.locations[v]
    }

    fn can_flip(&self, e: usize) -> Result<Result<(), FlipRefusal>, IntrinsicError> {
        IntegerTriangulation::can_flip(self, e).map_err(|x| to_intrinsic_error(x, "can_flip"))
    }

    fn flip_edge(&mut self, e: usize) -> Result<Flip, IntrinsicError> {
        IntegerTriangulation::flip_edge(self, e).map_err(|x| to_intrinsic_error(x, "flip_edge"))
    }

    /// Not provided yet: insertion in a face is work package WP4 (Supp. Alg. 5, 6, App. A).
    fn insert_vertex(&mut self, f: usize, _bary: [f64; 3]) -> Result<Insertion, IntrinsicError> {
        Err(IntrinsicError::Unsupported {
            representation: REPRESENTATION,
            operation: "insert_vertex",
            reason: format!(
                "inserting a vertex in face {f} needs the face split of Supp. Alg. 5-6 (work \
                 package WP4)"
            ),
        })
    }

    fn split_edge(&mut self, e: usize, t: f64) -> Result<Insertion, IntrinsicError> {
        IntegerTriangulation::split_edge(self, e, t)
            .map_err(|x| to_intrinsic_error(x, "split_edge"))
    }

    fn set_constrained(&mut self, e: usize, constrained: bool) -> Result<(), IntrinsicError> {
        self.complex.set_constrained(e, constrained)
    }

    /// Not provided yet: the positions of the crossings of an intrinsic edge with the input
    /// edges come from `ExtractIntrinsicEdge` (SPEC §3.1, work package WP2).
    fn trace_edge_on_input(
        &self,
        e: usize,
        _options: &EdgeTraceOptions,
    ) -> Result<EdgeTrace, IntrinsicError> {
        Err(IntrinsicError::Unsupported {
            representation: REPRESENTATION,
            operation: "trace_edge_on_input",
            reason: format!(
                "the path of edge {e} over the input mesh needs the geometric crossings of \
                 ExtractIntrinsicEdge (work package WP2)"
            ),
        })
    }
}
