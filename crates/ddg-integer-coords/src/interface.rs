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
//! of the shared kind (flips, splits, insertions, removals) as themselves, and operations this
//! crate does not provide as [`IntrinsicError::Unsupported`]. Insertion in a face is the face
//! split of WP4 ([`IntegerTriangulation::insert_vertex`]).
//! The path of an edge over the input mesh is `ExtractIntrinsicEdge`
//! ([`IntegerTriangulation::extract_intrinsic_edge`], WP2).

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

    /// The face split of Supp. Alg. 5, 6 with App. A (WP4):
    /// [`IntegerTriangulation::insert_vertex`], with the default options (SPEC I6 never
    /// refuses). Lengths and refusals as the signposts' (`plan_face_split`, I25).
    fn insert_vertex(&mut self, f: usize, bary: [f64; 3]) -> Result<Insertion, IntrinsicError> {
        IntegerTriangulation::insert_vertex(self, f, bary)
            .map_err(|x| to_intrinsic_error(x, "insert_vertex"))
    }

    fn split_edge(&mut self, e: usize, t: f64) -> Result<Insertion, IntrinsicError> {
        IntegerTriangulation::split_edge(self, e, t)
            .map_err(|x| to_intrinsic_error(x, "split_edge"))
    }

    fn set_constrained(&mut self, e: usize, constrained: bool) -> Result<(), IntrinsicError> {
        self.complex.set_constrained(e, constrained)
    }

    /// `ExtractIntrinsicEdge` (SPEC §3.1, FlipOut contract R2):
    /// [`IntegerTriangulation::extract_intrinsic_edge`]. The crossed input edges are exact; the
    /// mismatch and clamps are reported in the [`EdgeTrace`], never a reason to fail.
    fn trace_edge_on_input(
        &self,
        e: usize,
        options: &EdgeTraceOptions,
    ) -> Result<EdgeTrace, IntrinsicError> {
        self.extract_intrinsic_edge(e, options)
            .map_err(|err| to_intrinsic_error(err, "trace_edge_on_input"))
    }
}
