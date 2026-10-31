//! The interface shared by the representations of an intrinsic triangulation (SPEC N50).
//!
//! [`Triangulation`] is what code that only flips, splits, measures and maps back to the input
//! mesh programs against, so that it runs on every representation: the signposts of this crate
//! ([`IntrinsicTriangulation`]) and the integer coordinates of Gillespie, Sharp and Crane
//! (2021). Its connectivity is the shared [`Complex`]; lengths, angles, the place of every
//! vertex on the input mesh and the path of every edge over it come from the representation.
//!
//! Delaunay flipping ([`Triangulation::flip_to_delaunay`]) is provided here once, on top of the
//! trait's flips, so that every representation flips the same edges in the same order given
//! the same lengths (SPEC N6; the flip of [`IntrinsicTriangulation`] uses exactly this code).
//! What a flip, a face split and an edge split do to the lengths, and when they are refused, is
//! decided from the connectivity and the lengths alone by [`plan_flip`](crate::plan_flip),
//! [`plan_face_split`](crate::plan_face_split) and [`plan_edge_split`](crate::plan_edge_split);
//! [`IntrinsicTriangulation`] uses them, and another representation that uses them too refuses
//! the same operations and computes bitwise the same lengths.

use ddg_core::geometry::{corner_angle, corner_cotan, triangle_area};

use crate::complex::Complex;
use crate::delaunay::{DelaunayCriterion, DelaunayOptions, DelaunayReport};
use crate::edge_trace::{EdgeTrace, EdgeTraceOptions};
use crate::error::{FlipRefusal, IntrinsicError};
use crate::flip::{Flip, FlipShapeTest};
use crate::input::InputSurface;
use crate::insert::Insertion;
use crate::point::SurfacePoint;
use crate::refine::{RefineOptions, RefineReport, RefineRules};
use crate::remove::InsertedVertexRemoval;
use crate::trace::{Trace, TraceOptions};
use crate::triangulation::IntrinsicTriangulation;

/// An intrinsic triangulation of an input mesh, whatever its representation.
///
/// **Connectivity** is [`complex`](Self::complex) (ids and conventions: [`Complex`], SPEC N9).
/// **Guarantees** every implementation gives:
///
/// * A flip ([`flip_edge`](Self::flip_edge)) changes the connectivity exactly as
///   [`Complex::flip_edge`] describes, the length of the flipped edge by `it-flip-length`, and
///   nothing else: every other edge keeps its id, its length and its direction
///   ([`Complex::edge_halfedge`] may change id, never tail and head); vertex and face ids are
///   kept. Refusals are those of [`FlipRefusal`] (SPEC §4, N5), decided from the lengths alone,
///   so they are the same in every representation, plus constrained edges.
/// * A split ([`split_edge`](Self::split_edge), [`insert_vertex`](Self::insert_vertex))
///   changes the connectivity as [`Complex::split_edge`] and [`Complex::split_face`] describe,
///   and gives the new edges the lengths of SPEC N10 (`it-point-distances`; the halves of a
///   split edge `t ℓ` and `(1 − t) ℓ`), so representations stay bitwise identical in lengths.
/// * Constrained edges ([`Complex::is_constrained`]) are never flipped and count as Delaunay,
///   like boundary edges; splitting one gives two constrained halves.
/// * Nothing changes when an operation is refused.
pub trait Triangulation {
    /// The connectivity.
    fn complex(&self) -> &Complex;

    /// The input mesh `T⁰`. Its vertices are the vertices `0..n_input_vertices()` of the
    /// triangulation, with the same ids.
    fn input(&self) -> &InputSurface;

    /// The edge lengths, indexed by edge id.
    fn edge_lengths(&self) -> &[f64];

    /// The angle sum `Θ` of vertex `v`: that of the input mesh at an input vertex, `2π` at a
    /// vertex inserted in the interior, `π` at one inserted on the boundary.
    fn vertex_angle_sum(&self, v: usize) -> f64;

    /// The place of vertex `v` on the input mesh.
    fn vertex_location(&self, v: usize) -> SurfacePoint;

    /// Whether [`flip_edge`](Self::flip_edge) would flip edge `e`, without changing anything.
    ///
    /// # Errors
    /// [`IntrinsicError::IndexOutOfRange`]; a refusal is the `Err` of the inner result.
    fn can_flip(&self, e: usize) -> Result<Result<(), FlipRefusal>, IntrinsicError>;

    /// Flips edge `e` (see the guarantees above), with the shape test of SPEC N5
    /// ([`FlipShapeTest::N5`]). [`Flip::new_corner_angles`] are the corner angles at the new
    /// edge's end points that the flip creates (SPEC N52), bitwise those of
    /// [`corner_angle`](Self::corner_angle).
    ///
    /// # Errors
    /// [`IntrinsicError::IndexOutOfRange`], [`IntrinsicError::FlipRefused`] with the reason.
    fn flip_edge(&mut self, e: usize) -> Result<Flip, IntrinsicError>;

    /// Whether [`flip_edge_with`](Self::flip_edge_with) would flip edge `e` under `test`,
    /// without changing anything (SPEC N52). Provided: [`can_flip`](Self::can_flip) for the
    /// default test, [`IntrinsicError::Unsupported`] for any other; both representations of
    /// this workspace implement it.
    ///
    /// # Errors
    /// As [`can_flip`](Self::can_flip); [`IntrinsicError::Unsupported`].
    fn can_flip_with(
        &self,
        e: usize,
        test: FlipShapeTest,
    ) -> Result<Result<(), FlipRefusal>, IntrinsicError> {
        if test == FlipShapeTest::N5 {
            self.can_flip(e)
        } else {
            Err(unsupported_shape_test(e, test))
        }
    }

    /// Flips edge `e` as [`flip_edge`](Self::flip_edge) does, but refuses by the shape test
    /// `test` instead of that of SPEC N5 (SPEC N52; every other refusal and the arithmetic are
    /// the same, [`plan_flip_with`](crate::plan_flip_with)). Provided:
    /// [`flip_edge`](Self::flip_edge) for the default test, [`IntrinsicError::Unsupported`]
    /// for any other; both representations of this workspace implement it.
    ///
    /// # Errors
    /// As [`flip_edge`](Self::flip_edge); [`IntrinsicError::Unsupported`].
    fn flip_edge_with(&mut self, e: usize, test: FlipShapeTest) -> Result<Flip, IntrinsicError> {
        if test == FlipShapeTest::N5 {
            self.flip_edge(e)
        } else {
            Err(unsupported_shape_test(e, test))
        }
    }

    /// Inserts a vertex at the point `bary` of face `f` and joins it to the three corners
    /// (paper Alg. 7; ids as [`Complex::split_face`]).
    ///
    /// # Errors
    /// As [`IntrinsicTriangulation::insert_vertex`].
    fn insert_vertex(&mut self, f: usize, bary: [f64; 3]) -> Result<Insertion, IntrinsicError>;

    /// Splits edge `e` at the parameter `t` from the tail of its halfedge (paper Alg. 8; ids
    /// as [`Complex::split_edge`]). [`Insertion::halves`] gives both halves in order from the
    /// old edge's tail, directed like it; their lengths are `t ℓ` and `(1 − t) ℓ`.
    ///
    /// # Errors
    /// As [`IntrinsicTriangulation::split_edge`].
    fn split_edge(&mut self, e: usize, t: f64) -> Result<Insertion, IntrinsicError>;

    /// Marks edge `e` as constrained or not ([`Complex::set_constrained`]).
    ///
    /// # Errors
    /// [`IntrinsicError::IndexOutOfRange`].
    fn set_constrained(&mut self, e: usize, constrained: bool) -> Result<(), IntrinsicError>;

    /// The path of edge `e` over the input mesh, from the tail of its halfedge to its head: the
    /// input edges it crosses, in order ([`EdgeTrace`]).
    ///
    /// # Errors
    /// As [`IntrinsicTriangulation::trace_edge`].
    fn trace_edge_on_input(
        &self,
        e: usize,
        options: &EdgeTraceOptions,
    ) -> Result<EdgeTrace, IntrinsicError>;

    // ----- provided -----

    /// Number of vertices of the input mesh.
    fn n_input_vertices(&self) -> usize {
        self.input().mesh().n_vertices()
    }

    /// The length of edge `e`.
    fn edge_length(&self, e: usize) -> f64 {
        self.edge_lengths()[e]
    }

    /// The lengths of the three halfedges of face `f`, by slot.
    fn face_lengths(&self, f: usize) -> [f64; 3] {
        face_lengths(self.complex(), self.edge_lengths(), f)
    }

    /// The interior angle of corner `h` (at the tail of halfedge `h`), from the lengths.
    ///
    /// # Panics
    /// Never for a triangulation whose faces keep the strict triangle inequalities, which
    /// every operation does.
    fn corner_angle(&self, h: usize) -> f64 {
        let l = self.face_lengths(h / 3);
        let k = h % 3;
        corner_angle(l[(k + 1) % 3], l[k], l[(k + 2) % 3])
            .expect("every face keeps the triangle inequalities")
    }

    /// The area of face `f`, from the lengths.
    ///
    /// # Panics
    /// As [`corner_angle`](Self::corner_angle).
    fn face_area(&self, f: usize) -> f64 {
        let [a, b, c] = self.face_lengths(f);
        triangle_area(a, b, c).expect("every face keeps the triangle inequalities")
    }

    /// The sum of all face areas: the area of the input mesh up to rounding.
    fn total_area(&self) -> f64 {
        (0..self.complex().n_faces())
            .map(|f| self.face_area(f))
            .sum()
    }

    /// Whether edge `e` is constrained.
    fn is_constrained(&self, e: usize) -> bool {
        self.complex().is_constrained(e)
    }

    /// The sum of the two angles opposite edge `e` minus `π`; `None` for a boundary edge (see
    /// [`IntrinsicTriangulation::edge_angle_excess`]).
    fn edge_angle_excess(&self, e: usize) -> Option<f64> {
        crate::delaunay::edge_angle_excess(self.complex(), self.edge_lengths(), e)
    }

    /// The cotan weight of edge `e` (see [`IntrinsicTriangulation::edge_cotan_weight`]).
    fn edge_cotan_weight(&self, e: usize) -> f64 {
        crate::delaunay::edge_cotan_weight(self.complex(), self.edge_lengths(), e)
    }

    /// Whether every interior edge that is not constrained satisfies the criterion.
    fn is_delaunay(&self, criterion: &DelaunayCriterion) -> bool {
        let (complex, lengths) = (self.complex(), self.edge_lengths());
        (0..complex.n_edges())
            .all(|e| !crate::delaunay::fails_delaunay(complex, lengths, e, criterion))
    }

    /// Delaunay flipping (paper Sec. 4.1, SPEC N6) with this trait's flips; constrained edges
    /// are left alone and count as Delaunay (see
    /// [`IntrinsicTriangulation::flip_to_delaunay`], which is this function).
    ///
    /// # Errors
    /// [`IntrinsicError::InvalidOptions`] for a tolerance that is negative or not finite.
    fn flip_to_delaunay(
        &mut self,
        options: &DelaunayOptions,
    ) -> Result<DelaunayReport, IntrinsicError>
    where
        Self: Sized,
    {
        crate::delaunay::flip_to_delaunay(self, options)
    }

    /// Follows the straight path that leaves `start` in `direction` for `length` over the
    /// triangulation (the exponential map; SPEC N7, [`SurfacePoint`] for the frames), as
    /// [`IntrinsicTriangulation::trace`] does. Provided: [`AngleFrames`](crate::AngleFrames),
    /// whose signposts are computed from the corner angles (IC-SPEC I9);
    /// [`IntrinsicTriangulation`] uses its stored signposts instead. From a point of a face or
    /// an edge the path, its crossings and its end depend on the lengths alone, so both give
    /// the same result bit for bit (SPEC N51).
    ///
    /// # Errors
    /// As [`IntrinsicTriangulation::trace`].
    fn trace_geodesic(
        &self,
        start: SurfacePoint,
        direction: f64,
        length: f64,
        options: &TraceOptions,
    ) -> Result<Trace, IntrinsicError>
    where
        Self: Sized,
    {
        crate::frames::AngleFrames::of(self).trace(start, direction, length, options)
    }

    /// Removes the inserted interior vertex `v` of any degree: flips edges at it by
    /// [`removal_flip_choice`](crate::remove::removal_flip_choice) until its degree is three,
    /// then merges its three faces; ids are renumbered as [`Complex::remove_vertex`] says
    /// (SPEC N51; Gillespie, Sharp, Crane 2021, Alg. 8). Refinement by
    /// [`RefineRules::Chew`] needs it. A star in which no flip lowers
    /// the degree is [`RemovalRefusal::Stuck`](crate::RemovalRefusal::Stuck); every refusal
    /// leaves the triangulation unchanged. Provided: [`IntrinsicError::Unsupported`].
    ///
    /// # Errors
    /// [`IntrinsicError::IndexOutOfRange`], [`IntrinsicError::RemovalRefused`] (input vertex,
    /// boundary vertex, constrained edge, stuck, the refusals of [`Complex::remove_vertex`]),
    /// [`IntrinsicError::Unsupported`], and the representation's errors of its flips.
    fn remove_inserted_vertex(
        &mut self,
        v: usize,
    ) -> Result<InsertedVertexRemoval, IntrinsicError> {
        Err(IntrinsicError::Unsupported {
            representation: "this representation",
            operation: "remove_inserted_vertex",
            reason: format!("cannot remove vertex {v}: removal is not implemented"),
        })
    }

    /// Delaunay refinement by the given rules (SPEC N1, N4, N14 for [`RefineRules::Ruppert`],
    /// what [`IntrinsicTriangulation::refine`] does; IC-SPEC §5 for [`RefineRules::Chew`]),
    /// written once over this trait ([`crate::refine::refine`]), so every representation
    /// inserts the same points on the same input where its floating-point decisions agree.
    ///
    /// # Errors
    /// As [`IntrinsicTriangulation::refine`]; for [`RefineRules::Chew`] also the errors of
    /// [`remove_inserted_vertex`](Self::remove_inserted_vertex) other than refusals.
    fn refine_with_rules(
        &mut self,
        options: &RefineOptions,
        rules: RefineRules,
    ) -> Result<RefineReport, IntrinsicError>
    where
        Self: Sized,
    {
        crate::refine::refine(self, options, rules, &|_: &Self, _| false)
    }
}

/// The error of the provided [`Triangulation::flip_edge_with`] and
/// [`Triangulation::can_flip_with`] for a test other than the default.
fn unsupported_shape_test(e: usize, test: FlipShapeTest) -> IntrinsicError {
    IntrinsicError::Unsupported {
        representation: "this representation",
        operation: "flip_edge_with",
        reason: format!(
            "cannot flip edge {e} under the shape test {test:?}: only the default (SPEC N5) is \
             implemented"
        ),
    }
}

/// The lengths of the three halfedges of face `f`, by slot.
#[inline]
pub(crate) fn face_lengths(complex: &Complex, lengths: &[f64], f: usize) -> [f64; 3] {
    [
        lengths[complex.halfedge_edge[3 * f] as usize],
        lengths[complex.halfedge_edge[3 * f + 1] as usize],
        lengths[complex.halfedge_edge[3 * f + 2] as usize],
    ]
}

/// The cotangent of the corner opposite halfedge `h`.
#[inline]
pub(crate) fn opposite_cotan(complex: &Complex, lengths: &[f64], h: usize) -> f64 {
    let l = face_lengths(complex, lengths, h / 3);
    let k = h % 3;
    corner_cotan(l[k], l[(k + 1) % 3], l[(k + 2) % 3])
        .expect("every face keeps the triangle inequalities")
}

impl Triangulation for IntrinsicTriangulation {
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
        &self.edge_lengths
    }
    #[inline]
    fn vertex_angle_sum(&self, v: usize) -> f64 {
        self.angle_sums[v]
    }
    #[inline]
    fn vertex_location(&self, v: usize) -> SurfacePoint {
        self.locations[v]
    }
    fn can_flip(&self, e: usize) -> Result<Result<(), FlipRefusal>, IntrinsicError> {
        IntrinsicTriangulation::can_flip(self, e)
    }
    #[inline]
    fn flip_edge(&mut self, e: usize) -> Result<Flip, IntrinsicError> {
        IntrinsicTriangulation::flip_edge(self, e)
    }
    fn can_flip_with(
        &self,
        e: usize,
        test: FlipShapeTest,
    ) -> Result<Result<(), FlipRefusal>, IntrinsicError> {
        IntrinsicTriangulation::can_flip_with(self, e, test)
    }
    #[inline]
    fn flip_edge_with(&mut self, e: usize, test: FlipShapeTest) -> Result<Flip, IntrinsicError> {
        IntrinsicTriangulation::flip_edge_with(self, e, test)
    }
    fn insert_vertex(&mut self, f: usize, bary: [f64; 3]) -> Result<Insertion, IntrinsicError> {
        IntrinsicTriangulation::insert_vertex(self, f, bary)
    }
    fn split_edge(&mut self, e: usize, t: f64) -> Result<Insertion, IntrinsicError> {
        IntrinsicTriangulation::split_edge(self, e, t)
    }
    fn set_constrained(&mut self, e: usize, constrained: bool) -> Result<(), IntrinsicError> {
        IntrinsicTriangulation::set_constrained(self, e, constrained)
    }
    fn trace_edge_on_input(
        &self,
        e: usize,
        options: &EdgeTraceOptions,
    ) -> Result<EdgeTrace, IntrinsicError> {
        IntrinsicTriangulation::trace_edge(self, e, options)
    }
    #[inline]
    fn face_lengths(&self, f: usize) -> [f64; 3] {
        IntrinsicTriangulation::face_lengths(self, f)
    }
    fn corner_angle(&self, h: usize) -> f64 {
        IntrinsicTriangulation::corner_angle(self, h)
    }
    fn face_area(&self, f: usize) -> f64 {
        IntrinsicTriangulation::face_area(self, f)
    }
    fn trace_geodesic(
        &self,
        start: SurfacePoint,
        direction: f64,
        length: f64,
        options: &TraceOptions,
    ) -> Result<Trace, IntrinsicError> {
        IntrinsicTriangulation::trace(self, start, direction, length, options)
    }
    fn remove_inserted_vertex(
        &mut self,
        v: usize,
    ) -> Result<InsertedVertexRemoval, IntrinsicError> {
        IntrinsicTriangulation::remove_inserted_vertex(self, v)
    }
}
