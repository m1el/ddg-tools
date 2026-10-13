//! The error type of this crate.

use std::fmt;

use ddg_core::distance::DistanceError;
use ddg_core::{GeometryError, MeshError, SparseError};

/// Why a cut could not be evaluated.
///
/// These are errors of the *input*: the mesh, the level-set function, the problem description or
/// the options. Every message names the vertex, edge or face concerned. Things that are legal but
/// worth knowing (crossings that were moved away from a vertex, closed components without a cut)
/// are not errors; they are counted in the [`EvalReport`](crate::EvalReport).
#[derive(Debug, Clone, PartialEq)]
pub enum CutError {
    /// The geometry was not built for the mesh.
    Geometry(GeometryError),
    /// The cut mesh is not a valid manifold mesh. This cannot happen for a valid input mesh; it
    /// is reported instead of panicking.
    Mesh(MeshError),
    /// A sparse factorization or solve failed (out of memory, or a cut-mesh Laplacian that is
    /// numerically not positive definite because a cut triangle is nearly degenerate).
    Sparse(SparseError),
    /// The redistancing was refused. This cannot happen for a finite level-set function; it is
    /// reported instead of panicking.
    Distance(DistanceError),
    /// The total area of the mesh is not a positive finite number (edge lengths near the limits
    /// of `f64`), so the surface cannot be scaled to unit area.
    InvalidTotalArea {
        /// The computed total area.
        total_area: f64,
    },
    /// An array does not have one entry per vertex, face or region.
    WrongLength {
        /// Which array.
        what: &'static str,
        /// What it must have one entry for ("vertices", "faces", ...).
        per: &'static str,
        /// The required length.
        expected: usize,
        /// The length that was passed.
        got: usize,
    },
    /// A level-set value is NaN or infinite.
    NonFiniteLevelSet {
        /// The vertex.
        vertex: u32,
        /// Its value.
        value: f64,
    },
    /// The number of regions is not supported (this version handles exactly two).
    UnsupportedRegionCount {
        /// Number of entries of [`CutProblem::regions`](crate::CutProblem::regions).
        got: usize,
    },
    /// A face has vertices in three different regions. With two regions this cannot happen; the
    /// triple points of the `N`-region method are not implemented yet.
    UnsupportedTriplePoint {
        /// The face.
        face: u32,
        /// Its vertices.
        vertices: [u32; 3],
        /// Their regions.
        regions: [u32; 3],
    },
    /// A distortion weight `w_r` is negative or not finite, or an area weight `β_r` is not finite.
    InvalidRegionWeight {
        /// The region.
        region: u32,
        /// Which weight: `"distortion"` (`w_r`) or `"area"` (`β_r`).
        what: &'static str,
        /// Its value.
        value: f64,
    },
    /// The length weight `α_L` is negative or not finite.
    InvalidLengthWeight {
        /// The offending value.
        value: f64,
    },
    /// A penalty `η_i` is negative or not finite.
    InvalidPenalty {
        /// The vertex.
        vertex: u32,
        /// Its penalty.
        value: f64,
    },
    /// A weight `a_i` of the scaled Dirichlet energy is not positive and finite.
    InvalidScaleWeight {
        /// The vertex.
        vertex: u32,
        /// Its weight.
        value: f64,
    },
    /// A curvature override is not finite.
    InvalidCurvature {
        /// The face (`None` for a uniform override).
        face: Option<u32>,
        /// The offending density.
        value: f64,
    },
    /// A parameter of the curvature-based weight is out of range.
    InvalidFeatureSize {
        /// Which parameter and what it must satisfy.
        what: &'static str,
        /// Its value.
        value: f64,
    },
    /// The smoothed feature size of a vertex is not positive and finite, so `a = r³` is not a
    /// valid weight. The smoothing keeps `r` within the range of the unsmoothed feature sizes
    /// (SPEC D17), so this is reported only if rounding or overflow breaks that.
    NonPositiveFeatureSize {
        /// The vertex.
        vertex: u32,
        /// Its smoothed feature size.
        value: f64,
    },
    /// An option is out of range.
    InvalidOptions {
        /// Which option and why.
        reason: String,
    },
    /// A triangle of the cut mesh has no positive area: a crossing sits (numerically) on a
    /// vertex. Raise [`EvalOptions::min_crossing_offset`](crate::EvalOptions::min_crossing_offset).
    DegenerateCutTriangle {
        /// The face of the input mesh that contains the triangle.
        face: u32,
        /// The vertices of that face.
        vertices: [u32; 3],
        /// What is wrong.
        detail: String,
    },
    /// The pointwise energy density or its derivative returned a non-finite value.
    NonFiniteEnergyDensity {
        /// The vertex of the cut mesh at which it was evaluated.
        vertex: u32,
        /// The log scale factor `u` passed to the function.
        log_scale: f64,
        /// Which function: `"b"` or `"b'"`.
        what: &'static str,
        /// What it returned.
        value: f64,
    },
    /// An index passed to a helper is out of range.
    IndexOutOfRange {
        /// What kind of index.
        what: &'static str,
        /// The index.
        index: usize,
        /// The exclusive upper bound it violates.
        bound: usize,
    },
    /// A number of a flow step is NaN or infinite: a term of the energy, a speed, the step
    /// size, or a value of a level-set function after the step (weights, curvatures or step
    /// sizes near the limits of `f64`). The step is not taken, the flow stays in the state
    /// before it, and a run ends with this error instead of a stop reason.
    NonFiniteStep {
        /// The number of steps that the flow had taken before this one.
        step: usize,
        /// Which number.
        what: String,
        /// Its value.
        value: f64,
    },
}

impl fmt::Display for CutError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            CutError::Geometry(e) => write!(f, "{e}"),
            CutError::Mesh(e) => write!(f, "the cut mesh is not a valid mesh: {e}"),
            CutError::Sparse(e) => write!(f, "sparse linear algebra failed: {e}"),
            CutError::Distance(e) => write!(f, "redistancing failed: {e}"),
            CutError::InvalidTotalArea { total_area } => write!(
                f,
                "the total area of the mesh is {total_area}; it must be positive and finite to \
                 scale the surface to unit area"
            ),
            CutError::WrongLength {
                what,
                per,
                expected,
                got,
            } => write!(
                f,
                "{what} has {got} entries, but there are {expected} {per}"
            ),
            CutError::NonFiniteLevelSet { vertex, value } => write!(
                f,
                "level-set value of vertex {vertex} is {value}; it must be finite"
            ),
            CutError::UnsupportedRegionCount { got } => write!(
                f,
                "the problem describes {got} regions; this version supports exactly 2"
            ),
            CutError::UnsupportedTriplePoint {
                face,
                vertices: [a, b, c],
                regions: [ra, rb, rc],
            } => write!(
                f,
                "face {face} = ({a}, {b}, {c}) has vertices in three regions ({ra}, {rb}, {rc}); \
                 triple points are not supported yet"
            ),
            CutError::InvalidRegionWeight {
                region,
                what,
                value,
            } => {
                if *what == "distortion" {
                    write!(
                        f,
                        "distortion weight of region {region} is {value}; it must be finite and \
                         non-negative"
                    )
                } else {
                    write!(
                        f,
                        "{what} weight of region {region} is {value}; it must be finite"
                    )
                }
            }
            CutError::InvalidLengthWeight { value } => write!(
                f,
                "length weight is {value}; it must be finite and non-negative"
            ),
            CutError::InvalidPenalty { vertex, value } => write!(
                f,
                "penalty of vertex {vertex} is {value}; it must be finite and non-negative"
            ),
            CutError::InvalidScaleWeight { vertex, value } => write!(
                f,
                "scaled Dirichlet weight of vertex {vertex} is {value}; it must be finite and \
                 positive"
            ),
            CutError::InvalidCurvature { face, value } => match face {
                Some(face) => write!(
                    f,
                    "curvature override of face {face} is {value}; it must be finite"
                ),
                None => write!(
                    f,
                    "uniform curvature override is {value}; it must be finite"
                ),
            },
            CutError::InvalidFeatureSize { what, value } => {
                write!(f, "feature size: {what}, got {value}")
            }
            CutError::NonPositiveFeatureSize { vertex, value } => write!(
                f,
                "smoothed feature size of vertex {vertex} is {value}; it must be positive (reduce \
                 the smoothing or supply the weight per vertex)"
            ),
            CutError::InvalidOptions { reason } => write!(f, "invalid options: {reason}"),
            CutError::DegenerateCutTriangle {
                face,
                vertices: [a, b, c],
                detail,
            } => write!(
                f,
                "face {face} = ({a}, {b}, {c}) is cut into a triangle without positive area: \
                 {detail}"
            ),
            CutError::NonFiniteEnergyDensity {
                vertex,
                log_scale,
                what,
                value,
            } => write!(
                f,
                "pointwise energy: {what}({log_scale}) = {value} at cut-mesh vertex {vertex}; it \
                 must be finite"
            ),
            CutError::IndexOutOfRange { what, index, bound } => {
                write!(f, "{what} {index} is out of range (must be < {bound})")
            }
            CutError::NonFiniteStep { step, what, value } => write!(
                f,
                "after {step} steps of the flow {what} is {value}; it must be finite, and the \
                 next step was not taken"
            ),
        }
    }
}

impl std::error::Error for CutError {}

impl From<GeometryError> for CutError {
    fn from(e: GeometryError) -> Self {
        CutError::Geometry(e)
    }
}

impl From<MeshError> for CutError {
    fn from(e: MeshError) -> Self {
        CutError::Mesh(e)
    }
}

impl From<DistanceError> for CutError {
    fn from(e: DistanceError) -> Self {
        CutError::Distance(e)
    }
}

impl From<SparseError> for CutError {
    fn from(e: SparseError) -> Self {
        CutError::Sparse(e)
    }
}
