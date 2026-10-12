//! The error type of this crate.

use std::fmt;

use ddg_core::{GeometryError, SparseError};

/// Why a problem could not be set up or evaluated.
///
/// These are errors of the *input*. A solve that runs but does not reach a certified optimum is
/// not an error: it returns a [`Solution`](crate::Solution) whose
/// [`report.status`](crate::SolveReport::status) says so.
#[derive(Debug, Clone, PartialEq)]
pub enum ConeError {
    /// The geometry was not built for the mesh.
    Geometry(GeometryError),
    /// The mesh has more than one connected component. Split it first
    /// (`TriMesh::split_components`) and solve each component on its own.
    NotConnected {
        /// Number of connected components.
        components: usize,
    },
    /// `λ` is negative, NaN or infinite.
    InvalidLambda {
        /// The offending value.
        lambda: f64,
    },
    /// A per-vertex array does not have one entry per mesh vertex.
    WrongLength {
        /// Which array.
        what: &'static str,
        /// Number of mesh vertices.
        expected: usize,
        /// Length of the array.
        got: usize,
    },
    /// A distortion weight `wE` is negative, NaN or infinite.
    InvalidDistortionWeight {
        /// The vertex.
        vertex: u32,
        /// Its weight.
        value: f64,
    },
    /// A cone cost `wR` is negative or NaN (`+∞` is allowed and forbids a cone).
    InvalidConeCost {
        /// The vertex.
        vertex: u32,
        /// Its cost.
        value: f64,
    },
    /// Angle bounds must satisfy `lower ≤ 0 ≤ upper` (infinite values allowed, NaN not).
    InvalidBounds {
        /// The vertex.
        vertex: u32,
        /// Its lower bound.
        lower: f64,
        /// Its upper bound.
        upper: f64,
    },
    /// A prescribed cone is not usable.
    InvalidPrescribedCone {
        /// Position of the cone in the list.
        index: usize,
        /// Its vertex.
        vertex: u32,
        /// Its angle.
        angle: f64,
        /// What is wrong with it.
        reason: PrescribedConeIssue,
    },
    /// On a closed mesh every distortion weight is zero: the objective no longer depends on
    /// `u`, so its additive constant (and the solution) is undetermined.
    NoDistortionWeight,
    /// On a closed mesh the cone angles must sum to `2πχ` (Gauss–Bonnet), and the bounds,
    /// forbidden vertices and prescribed cones do not allow that.
    Infeasible {
        /// The required sum `Σ Ω_i = 2πχ`.
        required: f64,
        /// Smallest sum the constraints allow.
        min_sum: f64,
        /// Largest sum the constraints allow.
        max_sum: f64,
        /// Euler characteristic of the mesh.
        euler_characteristic: i64,
        /// Number of vertices at which a cone is forbidden.
        n_forbidden: usize,
        /// Number of prescribed cones.
        n_prescribed: usize,
        /// Sum of the prescribed cone angles.
        prescribed_sum: f64,
    },
    /// A cone configuration handed to [`evaluate_cones`](crate::evaluate_cones) is not usable.
    InvalidCone {
        /// Position of the cone in the list.
        index: usize,
        /// Its vertex.
        vertex: u32,
        /// Its angle.
        angle: f64,
        /// What is wrong with it.
        reason: PrescribedConeIssue,
    },
    /// On a closed mesh the cone angles handed to [`evaluate_cones`](crate::evaluate_cones)
    /// do not sum to `2πχ`, so no scale factor exists for them.
    ConeSumMismatch {
        /// The required sum `2πχ`.
        required: f64,
        /// The sum of the given angles.
        got: f64,
    },
    /// A value handed to a function is out of its documented range.
    InvalidArgument {
        /// Which argument and why.
        reason: String,
    },
    /// An option is out of range.
    InvalidOptions {
        /// Which option and why.
        reason: String,
    },
    /// A sparse factorization needed for the set-up failed (out of memory, or a Laplacian
    /// that is numerically not positive definite).
    Sparse(SparseError),
}

/// What is wrong with a prescribed (or evaluated) cone.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PrescribedConeIssue {
    /// The vertex index is not below the number of vertices.
    VertexOutOfRange,
    /// The vertex lies on the boundary, where `u = 0` is imposed and no cone can be placed.
    BoundaryVertex,
    /// The vertex already appeared earlier in the list.
    DuplicateVertex,
    /// The angle is NaN or infinite.
    NonFiniteAngle,
}

impl fmt::Display for ConeError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ConeError::Geometry(e) => write!(f, "{e}"),
            ConeError::NotConnected { components } => write!(
                f,
                "the mesh has {components} connected components; the solver needs a connected \
                 mesh (solve each component separately)"
            ),
            ConeError::InvalidLambda { lambda } => {
                write!(f, "lambda must be finite and non-negative, got {lambda}")
            }
            ConeError::WrongLength {
                what,
                expected,
                got,
            } => write!(
                f,
                "{what} has {got} entries, but the mesh has {expected} vertices"
            ),
            ConeError::InvalidDistortionWeight { vertex, value } => write!(
                f,
                "distortion weight of vertex {vertex} is {value}; it must be finite and \
                 non-negative"
            ),
            ConeError::InvalidConeCost { vertex, value } => write!(
                f,
                "cone cost of vertex {vertex} is {value}; it must be non-negative (+inf forbids \
                 a cone)"
            ),
            ConeError::InvalidBounds {
                vertex,
                lower,
                upper,
            } => write!(
                f,
                "angle bounds of vertex {vertex} are [{lower}, {upper}]; they must satisfy \
                 lower <= 0 <= upper"
            ),
            ConeError::InvalidPrescribedCone {
                index,
                vertex,
                angle,
                reason,
            } => write!(
                f,
                "prescribed cone {index} (vertex {vertex}, angle {angle}): {reason}"
            ),
            ConeError::NoDistortionWeight => write!(
                f,
                "all distortion weights are zero on a closed mesh: the additive constant of the \
                 scale factor is undetermined"
            ),
            ConeError::Infeasible {
                required,
                min_sum,
                max_sum,
                euler_characteristic,
                n_forbidden,
                n_prescribed,
                prescribed_sum,
            } => write!(
                f,
                "infeasible: on this closed mesh (Euler characteristic {euler_characteristic}) \
                 the cone angles must sum to 2 pi chi = {required}, but the angle bounds, the \
                 {n_forbidden} forbidden vertices and the {n_prescribed} prescribed cones \
                 (summing to {prescribed_sum}) only allow sums in [{min_sum}, {max_sum}]"
            ),
            ConeError::InvalidCone {
                index,
                vertex,
                angle,
                reason,
            } => write!(f, "cone {index} (vertex {vertex}, angle {angle}): {reason}"),
            ConeError::ConeSumMismatch { required, got } => write!(
                f,
                "on a closed mesh the cone angles must sum to 2 pi chi = {required}, but the \
                 given angles sum to {got}"
            ),
            ConeError::InvalidArgument { reason } => write!(f, "invalid argument: {reason}"),
            ConeError::InvalidOptions { reason } => write!(f, "invalid options: {reason}"),
            ConeError::Sparse(e) => write!(f, "sparse linear algebra failed during set-up: {e}"),
        }
    }
}

impl fmt::Display for PrescribedConeIssue {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            PrescribedConeIssue::VertexOutOfRange => "the vertex index is out of range",
            PrescribedConeIssue::BoundaryVertex => {
                "the vertex is on the boundary, where no cone can be placed"
            }
            PrescribedConeIssue::DuplicateVertex => "the vertex is listed more than once",
            PrescribedConeIssue::NonFiniteAngle => "the angle is not finite",
        })
    }
}

impl std::error::Error for ConeError {}

impl From<GeometryError> for ConeError {
    fn from(e: GeometryError) -> Self {
        ConeError::Geometry(e)
    }
}

impl From<SparseError> for ConeError {
    fn from(e: SparseError) -> Self {
        ConeError::Sparse(e)
    }
}
