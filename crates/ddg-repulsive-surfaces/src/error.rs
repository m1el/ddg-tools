//! The error type of the crate (`SPEC.md` §9.1).
//!
//! Every variant names the element that is wrong (a face, an edge, a vertex, a constraint) and
//! says what is wrong with it. Bad input never panics. The enum is `#[non_exhaustive]`; the
//! variants for constraints, solves and obstacles are declared here (WP1) for the later work
//! packages of `SPEC.md` §15.

use std::fmt;

use ddg_repulsive_curves::RcError;

/// All errors of the crate (`SPEC.md` §2.1, §2.3, §9.1).
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum RsError {
    // ---- topology (Surface::new, read_obj) ----
    /// The face list is empty.
    NoFaces,
    /// More vertices or faces than `u32` ids can address.
    TooLarge {
        /// Number of vertices requested.
        n_vertices: usize,
        /// Number of faces requested.
        n_faces: usize,
    },
    /// A face references a vertex index `>= n_vertices`.
    VertexOutOfRange {
        /// The face.
        face: u32,
        /// Corner slot (0, 1 or 2).
        corner: u8,
        /// The offending vertex index.
        vertex: u32,
        /// Number of vertices.
        n_vertices: usize,
    },
    /// A face uses the same vertex more than once.
    RepeatedVertex {
        /// The face.
        face: u32,
        /// Its three vertex indices.
        vertices: [u32; 3],
    },
    /// Two faces have the same three vertices (in any order).
    DuplicateFace {
        /// The lower-numbered face.
        first: u32,
        /// The higher-numbered face.
        second: u32,
        /// The shared vertices, sorted.
        vertices: [u32; 3],
    },
    /// An edge has more than two incident faces.
    NonManifoldEdge {
        /// The edge as `[min, max]`.
        edge: [u32; 2],
        /// All faces containing the edge, ascending.
        faces: Vec<u32>,
    },
    /// The faces around a vertex form more than one fan.
    NonManifoldVertex {
        /// The vertex.
        vertex: u32,
        /// Number of separate fans.
        fans: u32,
    },
    /// A vertex is used by no face (it has no area and would make the metric singular).
    UnreferencedVertex {
        /// The lowest-numbered unreferenced vertex.
        vertex: u32,
        /// Total number of unreferenced vertices.
        count: usize,
    },
    /// A polygon with other than three corners (only triangle meshes are accepted).
    PolygonFace {
        /// The face (0-based, in file order).
        face: u32,
        /// Its number of corners.
        corners: usize,
    },
    // ---- positions (Surface::validate_positions) ----
    /// The position array has the wrong length.
    PositionCountMismatch {
        /// Number of vertices of the surface.
        expected: usize,
        /// Length of the array that was passed.
        found: usize,
    },
    /// A coordinate is NaN or infinite.
    NonFinitePosition {
        /// The vertex.
        vertex: u32,
    },
    /// A face has zero (or non-finite) area: no normal is defined.
    DegenerateFace {
        /// The face.
        face: u32,
        /// Its area as computed.
        area: f64,
    },
    /// The energy (or its differential) is not finite; `faces` is a pair `(S, T)` whose term
    /// is not finite (coincident barycenters) or, if every term is finite, the pair with the
    /// largest term.
    NonFiniteEnergy {
        /// The pair of faces `(S, T)`.
        faces: [u32; 2],
    },
    // ---- options ----
    /// `Exponents::new` rejected `(alpha, beta)`.
    InvalidExponents {
        /// The exponent `α` that was given.
        alpha: f64,
        /// The exponent `β` that was given.
        beta: f64,
        /// The requirement that is violated.
        requirement: &'static str,
    },
    /// An option has a value outside its range.
    InvalidOption {
        /// Name of the option.
        name: &'static str,
        /// The value, printed.
        value: String,
        /// The requirement that is violated.
        requirement: &'static str,
    },
    /// A feature that a later work package provides was asked for.
    Unsupported {
        /// What was asked for.
        what: &'static str,
    },
    // ---- constraints (WP2) ----
    /// A component's translations are not fixed by any constraint (no barycenter, no pin).
    TranslationNotFixed {
        /// The component.
        component: u32,
    },
    /// A volume constraint or penalty on a component that is not consistently oriented.
    VolumeNeedsOrientation {
        /// The component.
        component: u32,
        /// An edge whose two faces traverse it in the same direction, as `[min, max]`.
        edge: [u32; 2],
    },
    /// A constraint cannot be satisfied (or is violated at the start beyond repair).
    InfeasibleConstraint {
        /// Index of the constraint in the list.
        constraint: usize,
        /// The violation, in the row's unit.
        violation: f64,
    },
    /// A constraint refers to an element (vertex, component, constraint) that does not exist.
    ConstraintOutOfRange {
        /// Index of the constraint (or schedule) in its list.
        constraint: usize,
        /// The offending index.
        index: u32,
    },
    /// A comparison method that needs a soft problem (AQP, L-BFGS, BCQN) was given hard
    /// constraints.
    MethodNeedsSoftProblem {
        /// The method.
        method: &'static str,
    },
    // ---- solves (WP2/WP6) ----
    /// The saddle-point matrix has a zero pivot.
    SingularSaddle {
        /// Index of the pivot.
        pivot: usize,
    },
    /// The iterative solver ended with a relative residual above `1e-1`.
    SolverDiverged {
        /// The relative residual reached.
        residual: f64,
    },
    // ---- obstacles (WP7) ----
    /// An obstacle is malformed.
    ObstacleError {
        /// Index of the obstacle.
        obstacle: usize,
        /// What is wrong.
        message: String,
    },
    // ---- I/O ----
    /// A file could not be read or written.
    Io {
        /// The path.
        path: String,
        /// The system's message.
        message: String,
    },
    /// A file is malformed.
    Parse {
        /// The file, if the text came from one.
        path: Option<String>,
        /// 1-based line number.
        line: usize,
        /// What is wrong.
        message: String,
    },
    // ---- other ----
    /// An error of the generic parts of `ddg-repulsive-curves` (BVH, block cluster tree, ...).
    Curves(RcError),
    /// The operation was cancelled.
    Cancelled,
}

fn write_list(f: &mut fmt::Formatter<'_>, items: &[u32]) -> fmt::Result {
    for (i, x) in items.iter().enumerate() {
        if i > 0 {
            f.write_str(", ")?;
        }
        write!(f, "{x}")?;
    }
    Ok(())
}

impl fmt::Display for RsError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        use RsError::*;
        match self {
            NoFaces => write!(f, "the surface has no faces"),
            TooLarge {
                n_vertices,
                n_faces,
            } => write!(
                f,
                "surface too large for 32-bit ids: {n_vertices} vertices, {n_faces} faces"
            ),
            VertexOutOfRange {
                face,
                corner,
                vertex,
                n_vertices,
            } => write!(
                f,
                "face {face} corner {corner} references vertex {vertex}, but the surface has \
                 {n_vertices} vertices"
            ),
            RepeatedVertex {
                face,
                vertices: [a, b, c],
            } => {
                let r = if a == b || a == c { a } else { b };
                write!(f, "face {face} = ({a}, {b}, {c}) repeats vertex {r}")
            }
            DuplicateFace {
                first,
                second,
                vertices: [a, b, c],
            } => write!(
                f,
                "faces {first} and {second} are duplicates: both have the vertices ({a}, {b}, {c})"
            ),
            NonManifoldEdge {
                edge: [a, b],
                faces,
            } => {
                write!(f, "edge ({a}, {b}) is shared by {} faces: ", faces.len())?;
                write_list(f, faces)
            }
            NonManifoldVertex { vertex, fans } => write!(
                f,
                "vertex {vertex} is non-manifold: its incident faces form {fans} separate fans"
            ),
            UnreferencedVertex { vertex, count } => write!(
                f,
                "vertex {vertex} is not referenced by any face ({count} unreferenced vertices \
                 in total)"
            ),
            PolygonFace { face, corners } => write!(
                f,
                "face {face} has {corners} corners: only triangle meshes are accepted \
                 (triangulate first)"
            ),
            PositionCountMismatch { expected, found } => write!(
                f,
                "expected {expected} positions (one per vertex), got {found}"
            ),
            NonFinitePosition { vertex } => {
                write!(f, "vertex {vertex} has a NaN or infinite coordinate")
            }
            DegenerateFace { face, area } => {
                write!(f, "face {face} is degenerate: its area is {area:e}")
            }
            NonFiniteEnergy { faces: [s, t] } => write!(
                f,
                "the tangent-point energy is not finite: the term of faces ({s}, {t}) \
                 (coincident or nearly coincident barycenters)"
            ),
            InvalidExponents {
                alpha,
                beta,
                requirement,
            } => write!(
                f,
                "invalid exponents (alpha {alpha}, beta {beta}): {requirement}"
            ),
            InvalidOption {
                name,
                value,
                requirement,
            } => write!(f, "invalid option {name} = {value}: {requirement}"),
            Unsupported { what } => write!(f, "not supported yet: {what}"),
            TranslationNotFixed { component } => write!(
                f,
                "the translations of component {component} are not fixed: add a barycenter \
                 constraint or a pin on it"
            ),
            VolumeNeedsOrientation {
                component,
                edge: [a, b],
            } => write!(
                f,
                "a volume needs a consistently oriented component, but the two faces of edge \
                 ({a}, {b}) in component {component} traverse it in the same direction"
            ),
            InfeasibleConstraint {
                constraint,
                violation,
            } => write!(
                f,
                "constraint {constraint} cannot be satisfied (violation {violation:e})"
            ),
            ConstraintOutOfRange { constraint, index } => write!(
                f,
                "constraint {constraint} refers to element {index}, which does not exist"
            ),
            MethodNeedsSoftProblem { method } => write!(
                f,
                "the method {method} needs a problem without hard constraints (use penalties)"
            ),
            SingularSaddle { pivot } => write!(
                f,
                "the saddle-point matrix has a zero pivot at index {pivot}"
            ),
            SolverDiverged { residual } => write!(
                f,
                "the linear solver diverged: relative residual {residual:e} (above 1e-1)"
            ),
            ObstacleError { obstacle, message } => write!(f, "obstacle {obstacle}: {message}"),
            Io { path, message } => write!(f, "{path}: {message}"),
            Parse {
                path: Some(path),
                line,
                message,
            } => write!(f, "{path}:{line}: {message}"),
            Parse {
                path: None,
                line,
                message,
            } => write!(f, "line {line}: {message}"),
            Curves(e) => write!(f, "{e}"),
            Cancelled => write!(f, "cancelled"),
        }
    }
}

impl std::error::Error for RsError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            RsError::Curves(e) => Some(e),
            _ => None,
        }
    }
}

impl From<RcError> for RsError {
    fn from(e: RcError) -> Self {
        RsError::Curves(e)
    }
}

impl From<ddg_core::io::FileError> for RsError {
    fn from(e: ddg_core::io::FileError) -> Self {
        use ddg_core::io::FileError;
        match e {
            FileError::Io { path, message } => RsError::Io { path, message },
            FileError::Parse {
                path,
                line,
                message,
            } => RsError::Parse {
                path,
                line,
                message,
            },
            FileError::InvalidInput { message } => RsError::Io {
                path: String::new(),
                message,
            },
        }
    }
}
