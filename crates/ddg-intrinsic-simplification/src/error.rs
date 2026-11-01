//! Errors and refusals.
//!
//! [`IemError`] is what goes wrong with the caller's input or the crate's own invariants;
//! [`RemovalFailure`] is why one vertex cannot be removed (SPEC §3.5: such a vertex is skipped,
//! it gets cost `∞`, nothing changes).

use std::fmt;

use ddg_core::{GeometryError, MeshError};
use ddg_intrinsic::FlipRefusal;

/// What went wrong.
#[derive(Debug, Clone, PartialEq)]
pub enum IemError {
    /// The input mesh is not a manifold oriented triangle mesh.
    Mesh(MeshError),
    /// The input lengths or positions are not valid (a face without positive area by its
    /// lengths, a non-finite value, a count mismatch).
    Geometry(GeometryError),
    /// An extra tracked point is not a point of a face of the input mesh.
    InvalidPoint {
        /// Its index among the extra points (its point id is this plus the number of input
        /// vertices).
        index: usize,
        /// What is wrong.
        what: String,
    },
    /// An option of [`simplify`](crate::simplify()) is not valid (a negative or non-finite
    /// weight, a per-vertex array of the wrong length, snapshot counts not descending, ...).
    InvalidOptions {
        /// What is wrong.
        what: String,
    },
    /// An index is out of range.
    IndexOutOfRange {
        /// What kind of element.
        what: &'static str,
        /// The index.
        index: usize,
        /// The number of such elements.
        bound: usize,
    },
    /// An edge flip was refused (topologically or geometrically, IT-SPEC N5).
    FlipRefused {
        /// The edge (internal id).
        edge: u32,
        /// Its end points (fine ids).
        vertices: [u32; 2],
        /// Why.
        reason: FlipRefusal,
    },
    /// A vertex cannot be removed (SPEC §3.5); the mesh is unchanged.
    Removal(RemovalFailure),
    /// A removal is pending ([`CoarseMesh::begin_removal`](crate::CoarseMesh::begin_removal)):
    /// finish or revert it first.
    RemovalPending,
    /// No removal is pending.
    NoRemovalPending,
    /// A check of [`CoarseMesh::certify`](crate::CoarseMesh::certify) failed: the numbered item
    /// of SPEC §2.4 and what is wrong.
    Certify {
        /// The item of SPEC §2.4 (1 to 7).
        item: u8,
        /// What is wrong, naming the element.
        what: String,
    },
    /// An operation met a state that its invariants exclude (for example a tracked point that
    /// lies outside both faces of a flip by more than the slack, SPEC S16). Reported, never
    /// repaired; the mesh is left as it was before the operation where that is possible (see
    /// the operation).
    Internal {
        /// What happened, naming the elements.
        what: String,
    },
}

impl fmt::Display for IemError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            IemError::Mesh(e) => write!(f, "invalid input mesh: {e}"),
            IemError::Geometry(e) => write!(f, "invalid input geometry: {e}"),
            IemError::InvalidPoint { index, what } => {
                write!(f, "extra tracked point {index}: {what}")
            }
            IemError::InvalidOptions { what } => {
                write!(f, "invalid simplification options: {what}")
            }
            IemError::IndexOutOfRange { what, index, bound } => {
                write!(f, "{what} {index} is out of range (there are {bound})")
            }
            IemError::FlipRefused {
                edge,
                vertices,
                reason,
            } => write!(
                f,
                "edge {edge} ({}, {}) cannot be flipped: {reason}",
                vertices[0], vertices[1]
            ),
            IemError::Removal(failure) => write!(f, "{failure}"),
            IemError::RemovalPending => write!(
                f,
                "a vertex removal is pending; finish or revert it before changing the mesh"
            ),
            IemError::NoRemovalPending => write!(f, "no vertex removal is pending"),
            IemError::Certify { item, what } => {
                write!(f, "certificate item {item} (SPEC 2.4) fails: {what}")
            }
            IemError::Internal { what } => write!(f, "internal error: {what}"),
        }
    }
}

impl std::error::Error for IemError {}

impl From<MeshError> for IemError {
    fn from(e: MeshError) -> Self {
        IemError::Mesh(e)
    }
}

impl From<GeometryError> for IemError {
    fn from(e: GeometryError) -> Self {
        IemError::Geometry(e)
    }
}

/// Why a vertex could not be removed (SPEC §3.5). The mesh is unchanged (bitwise) when a
/// removal fails.
#[derive(Debug, Clone, PartialEq)]
pub struct RemovalFailure {
    /// The vertex (fine id).
    pub vertex: u32,
    /// Why.
    pub reason: RemovalRefusal,
}

/// The reasons of [`RemovalFailure`].
#[derive(Debug, Clone, PartialEq)]
pub enum RemovalRefusal {
    /// The vertex is not a vertex of the coarse mesh (already removed, or never one).
    NotAVertex,
    /// The vertex is fixed (SPEC §2.2, paper Sec. 8.6.3).
    Fixed,
    /// A face at the vertex has one distinct vertex (paper Sec. 4, "self-faces ... are
    /// skipped").
    SelfFace {
        /// The face (internal id).
        face: u32,
    },
    /// The vertex is an end point of a boundary self-edge: it is the last vertex of its
    /// boundary loop (paper Sec. 4).
    BoundarySelfEdge,
    /// An interior vertex on fewer than three faces cannot be flipped to degree three.
    InteriorDegree {
        /// The number of faces at it.
        faces: usize,
    },
    /// An ear (a boundary vertex on one face) whose opposite edge cannot be flipped (paper
    /// Sec. 4, "Special cases").
    EarNotFlippable {
        /// Why the opposite edge cannot be flipped.
        reason: FlipRefusal,
    },
    /// The flattening did not converge (SPEC §3.4): no root was bracketed, or the iteration
    /// cap was reached.
    FlattenFailed {
        /// The best `|Θ − Θ̃|` reached, in radians.
        residual: f64,
        /// Newton and bisection steps taken in the last attempt.
        iterations: usize,
    },
    /// The flattened star violates the triangle inequality and Springborn's flips (SPEC §3.4)
    /// could not repair it.
    NotFlattenable {
        /// The number of fallback flips made before giving up.
        fallback_flips: usize,
    },
    /// No spoke can be flipped and the vertex still has more faces than the removal needs
    /// (three inside, two on the boundary; SPEC §3.5, S7).
    Stuck {
        /// The number of faces at the vertex when no flip was left.
        faces: usize,
    },
    /// At the degree the removal needs, an edge at the vertex joins it to itself.
    SelfEdge,
    /// The triangle that would replace the star has no positive area by its lengths.
    Degenerate {
        /// Its three lengths.
        lengths: [f64; 3],
    },
    /// The removal met a state its invariants exclude (a tracked point outside both faces of a
    /// flip, SPEC S16), or another removal is pending. Everything was reverted. Never expected:
    /// tests assert that it does not occur.
    Internal {
        /// What happened.
        what: String,
    },
}

impl fmt::Display for RemovalRefusal {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            RemovalRefusal::NotAVertex => write!(f, "it is not a vertex of the coarse mesh"),
            RemovalRefusal::Fixed => write!(f, "it is fixed"),
            RemovalRefusal::SelfFace { face } => {
                write!(f, "face {face} at it has one distinct vertex (a self-face)")
            }
            RemovalRefusal::BoundarySelfEdge => write!(
                f,
                "it is on a boundary self-edge (the last vertex of its boundary loop)"
            ),
            RemovalRefusal::InteriorDegree { faces } => write!(
                f,
                "it is an interior vertex on {faces} faces, which cannot be flipped to three"
            ),
            RemovalRefusal::EarNotFlippable { reason } => write!(
                f,
                "it is an ear whose opposite edge cannot be flipped: {reason}"
            ),
            RemovalRefusal::FlattenFailed {
                residual,
                iterations,
            } => write!(
                f,
                "its flattening did not converge (residual {residual:e} rad after {iterations} \
                 steps)"
            ),
            RemovalRefusal::NotFlattenable { fallback_flips } => write!(
                f,
                "its flattened star violates the triangle inequality and {fallback_flips} \
                 fallback flips did not repair it"
            ),
            RemovalRefusal::Stuck { faces } => write!(
                f,
                "it is on {faces} faces and no edge at it can be flipped to reduce that"
            ),
            RemovalRefusal::SelfEdge => write!(
                f,
                "at the degree the removal needs, an edge joins it to itself"
            ),
            RemovalRefusal::Degenerate { lengths } => write!(
                f,
                "the triangle replacing its star would have the lengths {lengths:?}, which form \
                 no triangle of positive area"
            ),
            RemovalRefusal::Internal { what } => write!(f, "internal error: {what}"),
        }
    }
}

impl fmt::Display for RemovalFailure {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "vertex {} cannot be removed: {}",
            self.vertex, self.reason
        )
    }
}

impl std::error::Error for RemovalFailure {}
