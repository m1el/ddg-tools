//! The error type of this crate (SPEC §7.1).

use std::fmt;

use ddg_bff::BffError;
use ddg_core::geometry::GeometryError;
use ddg_core::topology::MeshError;

use crate::eig::EigError;

/// Why an evaluation, optimization, subdivision, remeshing or flattening could not be done.
///
/// Every variant names the offending element (vertex, face, option) so that the message can be
/// shown to a user as is. Nothing in this crate panics on bad input; it returns one of these.
#[derive(Debug, Clone, PartialEq)]
pub enum DevError {
    /// A mesh failed validation (input, or a mesh built by this crate).
    Mesh(MeshError),
    /// A geometry could not be built.
    Geometry(GeometryError),
    /// Boundary First Flattening or the cut planner failed.
    Bff(BffError),
    /// An eigensolver failed (non-finite matrix entry).
    Eig(EigError),
    /// An array has the wrong length.
    LengthMismatch {
        /// What the array is (e.g. `"positions"`).
        what: &'static str,
        /// The length it must have.
        expected: usize,
        /// The length it has.
        got: usize,
    },
    /// A vertex position has a NaN or infinite coordinate.
    NonFinitePosition {
        /// The vertex.
        vertex: u32,
    },
    /// A face has a zero or non-finite normal `(f1 − f0) × (f2 − f0)`; its unit normal is undefined.
    DegenerateFace {
        /// The face.
        face: u32,
    },
    /// A face has a corner angle below the admissibility floor (SPEC §4.9, D31): it is
    /// degenerate for all practical purposes (a needle or a cap of zero area), although its normal
    /// is not exactly zero.
    NearDegenerateFace {
        /// The face.
        face: u32,
        /// Its smallest corner angle (radians).
        min_angle: f64,
        /// The floor (radians).
        limit: f64,
    },
    /// The energy is undefined at a vertex that counts toward the total (SPEC §3.6, D12).
    UndefinedEnergy {
        /// The vertex.
        vertex: u32,
        /// Why.
        reason: &'static str,
    },
    /// An option has an invalid value.
    InvalidOption {
        /// The option.
        name: &'static str,
        /// Its value (converted to `f64`).
        value: f64,
        /// What it must satisfy.
        requirement: &'static str,
    },
    /// A vertex id is out of range.
    VertexOutOfRange {
        /// Where the id came from (e.g. `"pinned"`).
        what: &'static str,
        /// The offending id.
        vertex: u32,
        /// Number of vertices of the mesh.
        n_vertices: usize,
    },
    /// An internal invariant failed. This is a bug; it is reported instead of panicking.
    Internal(String),
}

impl fmt::Display for DevError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            DevError::Mesh(e) => write!(f, "invalid mesh: {e}"),
            DevError::Geometry(e) => write!(f, "invalid geometry: {e}"),
            DevError::Bff(e) => write!(f, "boundary first flattening failed: {e}"),
            DevError::Eig(e) => write!(f, "eigensolver failed: {e}"),
            DevError::LengthMismatch {
                what,
                expected,
                got,
            } => {
                write!(
                    f,
                    "{what}: expected {expected} entries (one per element), got {got}"
                )
            }
            DevError::NonFinitePosition { vertex } => {
                write!(
                    f,
                    "vertex {vertex} has a non-finite (NaN or infinite) position"
                )
            }
            DevError::DegenerateFace { face } => write!(
                f,
                "face {face} is degenerate: its normal (f1 - f0) x (f2 - f0) is zero or not finite, \
                 so its unit normal and corner angles are undefined"
            ),
            DevError::NearDegenerateFace { face, min_angle, limit } => write!(
                f,
                "face {face} is nearly degenerate: its smallest corner angle is {min_angle:e} rad, below the floor {limit:e} rad"
            ),
            DevError::UndefinedEnergy { vertex, reason } => {
                write!(f, "the energy is undefined at vertex {vertex}: {reason}")
            }
            DevError::InvalidOption {
                name,
                value,
                requirement,
            } => {
                write!(
                    f,
                    "option `{name}` = {value} is invalid: it must satisfy {requirement}"
                )
            }
            DevError::VertexOutOfRange {
                what,
                vertex,
                n_vertices,
            } => write!(
                f,
                "{what}: vertex id {vertex} is out of range (the mesh has {n_vertices} vertices)"
            ),
            DevError::Internal(msg) => write!(f, "internal error (a bug): {msg}"),
        }
    }
}

impl std::error::Error for DevError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            DevError::Mesh(e) => Some(e),
            DevError::Geometry(e) => Some(e),
            DevError::Bff(e) => Some(e),
            DevError::Eig(e) => Some(e),
            _ => None,
        }
    }
}

impl From<MeshError> for DevError {
    fn from(e: MeshError) -> Self {
        DevError::Mesh(e)
    }
}

impl From<GeometryError> for DevError {
    fn from(e: GeometryError) -> Self {
        DevError::Geometry(e)
    }
}

impl From<BffError> for DevError {
    fn from(e: BffError) -> Self {
        DevError::Bff(e)
    }
}

impl From<EigError> for DevError {
    fn from(e: EigError) -> Self {
        DevError::Eig(e)
    }
}
